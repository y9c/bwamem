# Makefile for bwamem project

.PHONY: help install install-dev test test-cov build clean lint format docs check all bwa-lib

# OS detection for sed compatibility
OS := $(shell uname)
ifeq ($(OS), Darwin)
SEDI=sed -i '.bak'
else
SEDI=sed -i
endif

# Detect if uv is available
UV := $(shell command -v uv 2> /dev/null)

help:  ## Show this help message
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install:  ## Install the package and dependencies
	$(MAKE) bwa-lib
ifdef UV
	uv sync
	uv pip install -e . --no-deps
else
	python3 -m pip install -e .
endif

install-dev:  ## Install with development dependencies
	$(MAKE) bwa-lib
ifdef UV
	uv sync --extra dev
	uv pip install -e . --no-deps
else
	python3 -m pip install -e ".[dev]"
endif

bwa-lib: bwa/libbwa.a  ## Build BWA static library

bwa/libbwa.a:  ## Compile BWA C library with compilation flags
	${SEDI} 's/CFLAGS=.*/CFLAGS=-g\ -Wall\ -Wno-unused-function\ -O2\ -fPIC\ -fno-finite-math-only/' bwa/Makefile
	cd bwa && make libbwa.a

test:  ## Run tests
ifdef UV
	uv run pytest tests/ -v
else
	python3 -m pytest tests/ -v
endif

test-cov:  ## Run tests with coverage
ifdef UV
	uv run pytest tests/ --cov=bwamem --cov-report=html --cov-report=term
else
	python3 -m pytest tests/ --cov=bwamem --cov-report=html --cov-report=term
endif

build:  ## Build the package (wheel and sdist)
	$(MAKE) bwa-lib
ifdef UV
	uv build
else
	python3 -m build
endif

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
ifdef UV
	uv run ruff check bwamem/
	uv run ruff format --check bwamem/
else
	python3 -m ruff check bwamem/ || echo "ruff not installed, skipping lint"
	python3 -m ruff format --check bwamem/ || echo "ruff not installed, skipping format check"
endif

format:  ## Format code
ifdef UV
	uv run ruff format bwamem/
	uv run ruff check --fix bwamem/
else
	python3 -m ruff format bwamem/ || echo "ruff not installed, skipping format"
	python3 -m ruff check --fix bwamem/ || echo "ruff not installed, skipping check"
endif

docs:  ## Build documentation
	@echo "Documentation is in the docs/ directory"
	@echo "README.md contains the main documentation"

check: lint test  ## Run all checks (lint + tests)

all: clean install-dev test build  ## Run full pipeline

rebuild: clean build  ## Clean and rebuild package

.PHONY: rebuild
