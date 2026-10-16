.PHONY: install run test test-fast lint format clean

install:
	pip install -e ".[dev]"

run:
	pathogenomics -v run

test:
	pytest

test-fast:
	pytest -m "not slow"

lint:
	ruff check src tests
	ruff format --check src tests

format:
	ruff format src tests
	ruff check --fix src tests

clean:
	rm -rf .pytest_cache .ruff_cache build dist src/*.egg-info
