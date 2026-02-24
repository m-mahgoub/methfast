use clap::Parser;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

#[derive(Debug, Clone)]
struct MethInterval {
    start: i32,
    end: i32,
    fraction: f32,
    coverage: i32,
}

#[derive(Debug)]
struct MethRanges {
    by_chrom: HashMap<String, Vec<MethInterval>>,
}

#[derive(Debug)]
struct TargetInterval {
    chrom: String,
    start: i32,
    end: i32,
}

#[derive(Parser, Debug)]
#[command(
    name = "methfast",
    version,
    about = "Extract weighted methylation values for target BED intervals."
)]
struct Cli {
    #[arg(value_name = "METHYLATION_BED")]
    methylation_bed: PathBuf,
    #[arg(value_name = "TARGET_BED")]
    target_bed: PathBuf,

    #[arg(short = 'f', long = "fraction-col", default_value_t = 4)]
    frac_col: usize,
    #[arg(short = 'c', long = "coverage-col", default_value_t = 5)]
    cov_col: usize,
    #[arg(short = 'm', long = "methylated-col", default_value_t = 0)]
    meth_col: usize,
    #[arg(short = 'u', long = "unmethylated-col", default_value_t = 0)]
    unmeth_col: usize,
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,
    #[arg(
        short = 't',
        long = "threads",
        help = "Number of worker threads for processing target intervals"
    )]
    threads: Option<usize>,
}

fn parse_i32_lossy(s: &str) -> i32 {
    s.parse::<i32>().unwrap_or(0)
}

fn parse_f32_lossy(s: &str) -> f32 {
    s.parse::<f32>().unwrap_or(0.0)
}

fn is_gzipped(path: &PathBuf) -> Result<bool, Box<dyn Error>> {
    let mut file = File::open(path)?;
    let mut header = [0_u8; 3];
    let n = file.read(&mut header)?;
    if n < 3 {
        return Ok(false);
    }
    Ok(header == [0x1F, 0x8B, 0x08])
}

fn open_maybe_gz(path: &PathBuf) -> Result<Box<dyn BufRead>, Box<dyn Error>> {
    if is_gzipped(path)? {
        let file = File::open(path)?;
        let decoder = MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        let file = File::open(path)?;
        Ok(Box::new(BufReader::new(file)))
    }
}

fn parse_meth_bed(
    path: &PathBuf,
    frac_col: usize,
    cov_col: usize,
    meth_col: usize,
    unmeth_col: usize,
) -> Result<MethRanges, Box<dyn Error>> {
    let mut by_chrom: HashMap<String, Vec<MethInterval>> = HashMap::new();
    let mut reader = open_maybe_gz(path)?;
    let mut line = String::new();

    let mut prev_chrom = String::new();
    let mut prev_start: i32 = -1;
    let mut prev_end: i32 = -1;
    let mut linenum: usize = 0;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        linenum += 1;

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 4 {
            continue;
        }

        let chrom = fields[0].to_string();
        let start = parse_i32_lossy(fields[1]);
        let end = parse_i32_lossy(fields[2]);

        if prev_start != -1 && chrom == prev_chrom && start < prev_end {
            return Err(format!(
                "Error: Methylation BED file is not sorted. Exiting...\nLine {}: {} {} {}, then {} {} {}",
                linenum, prev_chrom, prev_start, prev_end, chrom, start, end
            )
            .into());
        }

        let field_count = fields.len();
        let (fraction, coverage) = if meth_col > 0
            && meth_col <= field_count
            && unmeth_col > 0
            && unmeth_col <= field_count
        {
            let methylated = parse_i32_lossy(fields[meth_col - 1]);
            let unmethylated = parse_i32_lossy(fields[unmeth_col - 1]);
            let coverage = methylated + unmethylated;
            let fraction = if coverage > 0 {
                methylated as f32 / coverage as f32
            } else {
                0.0
            };
            (fraction, coverage)
        } else if meth_col > 0 && meth_col <= field_count && cov_col > 0 && cov_col <= field_count {
            let methylated = parse_i32_lossy(fields[meth_col - 1]);
            let coverage = parse_i32_lossy(fields[cov_col - 1]);
            let fraction = if coverage > 0 {
                methylated as f32 / coverage as f32
            } else {
                0.0
            };
            (fraction, coverage)
        } else if cov_col > 0 && cov_col <= field_count && frac_col > 0 && frac_col <= field_count {
            let fraction = parse_f32_lossy(fields[frac_col - 1]);
            let coverage = parse_i32_lossy(fields[cov_col - 1]);
            (fraction, coverage)
        } else {
            return Err("Error: invalid column indices".into());
        };

        by_chrom
            .entry(chrom.clone())
            .or_default()
            .push(MethInterval {
                start,
                end,
                fraction,
                coverage,
            });

        prev_chrom = chrom;
        prev_start = start;
        prev_end = end;
    }

    Ok(MethRanges { by_chrom })
}

fn parse_targets(path: &PathBuf) -> Result<Vec<TargetInterval>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut targets = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let mut toks = line.split('\t');
        let Some(chrom) = toks.next() else {
            continue;
        };
        let Some(start_s) = toks.next() else {
            continue;
        };
        let Some(end_s) = toks.next() else {
            continue;
        };

        targets.push(TargetInterval {
            chrom: chrom.to_string(),
            start: parse_i32_lossy(start_s),
            end: parse_i32_lossy(end_s),
        });
    }

    Ok(targets)
}

fn lower_bound_end(intervals: &[MethInterval], start: i32) -> usize {
    let mut lo = 0_usize;
    let mut hi = intervals.len();
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if intervals[mid].end <= start {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

fn compute_target_line(ranges: &MethRanges, target: &TargetInterval) -> String {
    let mut num_positions = 0_usize;
    let mut sum_total_coverage = 0_i32;
    let mut sum_meth_coverage = 0_f32;

    if let Some(intervals) = ranges.by_chrom.get(&target.chrom) {
        let idx = lower_bound_end(intervals, target.start);
        for iv in &intervals[idx..] {
            if iv.start >= target.end {
                break;
            }
            if iv.end > target.start {
                num_positions += 1;
                sum_total_coverage += iv.coverage;
                sum_meth_coverage += iv.fraction * iv.coverage as f32;
            }
        }
    }

    let weighted_fraction = if sum_total_coverage > 0 {
        sum_meth_coverage / sum_total_coverage as f32
    } else {
        0.0
    };

    format!(
        "{}\t{}\t{}\t{}\t{}\t{:.4}",
        target.chrom,
        target.start,
        target.end,
        num_positions,
        sum_total_coverage,
        weighted_fraction
    )
}

fn run(cli: Cli) -> Result<(), Box<dyn Error>> {
    if let Some(threads) = cli.threads {
        if threads > 0 {
            let _ = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global();
        }
    }

    let ranges = parse_meth_bed(
        &cli.methylation_bed,
        cli.frac_col,
        cli.cov_col,
        cli.meth_col,
        cli.unmeth_col,
    )?;
    let targets = parse_targets(&cli.target_bed)?;
    let lines: Vec<String> = targets
        .par_iter()
        .map(|target| compute_target_line(&ranges, target))
        .collect();

    match cli.output {
        Some(path) => {
            let mut out = BufWriter::new(File::create(path)?);
            for line in &lines {
                writeln!(out, "{line}")?;
            }
            out.flush()?;
        }
        None => {
            let stdout = std::io::stdout();
            let mut out = BufWriter::new(stdout.lock());
            for line in &lines {
                writeln!(out, "{line}")?;
            }
            out.flush()?;
        }
    }

    Ok(())
}

fn main() {
    let cli = Cli::parse();
    if let Err(err) = run(cli) {
        eprintln!("{err}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn computes_weighted_fraction_from_intervals() {
        let mut by_chrom: HashMap<String, Vec<MethInterval>> = HashMap::new();
        by_chrom.insert(
            "chr1".to_string(),
            vec![
                MethInterval {
                    start: 10,
                    end: 11,
                    fraction: 1.0,
                    coverage: 5,
                },
                MethInterval {
                    start: 12,
                    end: 13,
                    fraction: 0.5,
                    coverage: 10,
                },
                MethInterval {
                    start: 20,
                    end: 21,
                    fraction: 0.0,
                    coverage: 3,
                },
            ],
        );

        let ranges = MethRanges { by_chrom };
        let target = TargetInterval {
            chrom: "chr1".to_string(),
            start: 9,
            end: 14,
        };
        let line = compute_target_line(&ranges, &target);
        assert_eq!(line, "chr1\t9\t14\t2\t15\t0.6667");
    }

    #[test]
    fn finds_first_candidate_interval_with_binary_search() {
        let intervals = vec![
            MethInterval {
                start: 1,
                end: 2,
                fraction: 0.0,
                coverage: 1,
            },
            MethInterval {
                start: 5,
                end: 6,
                fraction: 0.0,
                coverage: 1,
            },
            MethInterval {
                start: 10,
                end: 11,
                fraction: 0.0,
                coverage: 1,
            },
        ];
        assert_eq!(lower_bound_end(&intervals, 0), 0);
        assert_eq!(lower_bound_end(&intervals, 2), 1);
        assert_eq!(lower_bound_end(&intervals, 6), 2);
        assert_eq!(lower_bound_end(&intervals, 11), 3);
    }
}
