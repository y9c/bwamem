# Makefile for bwamem project

.PHONY: help install install-dev test test-cov build clean lint format docs check all bwa-lib

# OS detection for sed compatibility
OS := $(shell uname)
ifeq ($(OS), Darwin)
SEDI=sed -i '.bak'
else
SEDI=sed -i
endif

help:  ## Show this help message
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install:  ## Install the package and dependencies
	$(MAKE) bwa-lib
	uv sync
	uv pip install -e . --no-deps

install-dev:  ## Install with development dependencies
	$(MAKE) bwa-lib
	uv sync --extra dev
	uv pip install -e . --no-deps

bwa-lib: bwa/libbwa.a  ## Build BWA static library

bwa/libbwa.a: patches/bwa-makefile-cflags.patch  ## Compile BWA C library with compilation flags
	@echo "Applying patches to BWA source files..."
	@# Check if bwa is a git repo (submodule)
	@if [ -d bwa/.git ] || [ -f bwa/.git ]; then \
		echo "BWA is a git submodule, using git to apply patches..."; \
		cd bwa && (git diff --quiet Makefile || git checkout Makefile) && \
			git apply ../patches/bwa-makefile-cflags.patch; \
	else \
		echo "BWA is not a git repo, using patch command..."; \
		patch -f -p0 -d bwa < patches/bwa-makefile-cflags.patch || true; \
	fi
	cd bwa && $(MAKE) libbwa.a
	@echo "Reverting patches from BWA source files..."
	@# Revert patches to keep bwa directory clean
	@if [ -d bwa/.git ] || [ -f bwa/.git ]; then \
		cd bwa && git checkout Makefile >/dev/null 2>&1 || true; \
	else \
		patch -R -f -p0 -d bwa < patches/bwa-makefile-cflags.patch >/dev/null 2>&1 || true; \
	fi

test:  ## Run tests
	uv run pytest tests/ -v

test-cov:  ## Run tests with coverage
	uv run pytest tests/ --cov=bwamem --cov-report=html --cov-report=term

build:  ## Build the package (wheel and sdist)
	$(MAKE) bwa-lib
	uv build

clean:  ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf bwalib.cpython*.so
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	cd bwa && make clean || true

clean-all: clean  ## Clean everything including BWA artifacts
	rm -f bwa/*.o bwa/*.a bwa/bwa

lint:  ## Run linting
	uv run ruff check bwamem/
	uv run ruff format --check bwamem/

format:  ## Format code
	uv run ruff format bwamem/
	uv run ruff check --fix bwamem/

docs:  ## Build documentation
	@echo "Documentation is in the docs/ directory"
	@echo "README.md contains the main documentation"

check: lint test  ## Run all checks (lint + tests)

all: clean install-dev test build  ## Run full pipeline

rebuild: clean build  ## Clean and rebuild package

.PHONY: rebuild
