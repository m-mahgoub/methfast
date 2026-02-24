.PHONY: fmt lint test ci install uninstall

fmt:
	cargo fmt --all

lint:
	cargo clippy --all-targets --all-features -- -D warnings

test:
	cargo test --all-targets --all-features

ci: fmt lint test

install:
	cargo install --path . --locked

uninstall:
	cargo uninstall methfast || true
