# Makefile for Euler Student Guide
# Build LaTeX document and run Python tests

SHELL := /bin/bash

# Python settings
PYTHON := python3
PYTEST := python3 -m pytest
SRC_DIR := src
TEST_DIR := tests

# LaTeX settings
LATEX := pdflatex
LATEX_OPTS := -interaction=nonstopmode -halt-on-error
BIBTEX := bibtex
MAIN_TEX := Euler_Student_Guide.tex
OUTPUT_PDF := $(MAIN_TEX:.tex=.pdf)

# Directories
BUILD_DIR := build
CHAPTERS_DIR := chapters
FIGURES_DIR := figures

.PHONY: all test clean pdf view help install

# Default target
all: pdf test

# Install package in development mode
install:
	pip install -e .

# Run all tests
test:
	$(PYTEST) $(TEST_DIR) -v

# Run tests with coverage
test-cov:
	$(PYTEST) $(TEST_DIR) -v --cov=$(SRC_DIR)/euler --cov-report=term-missing

# Run specific test module
test-module:
	$(PYTEST) $(TEST_DIR)/test_euler.py -v -k $(MODULE)

# Build LaTeX document
pdf: $(OUTPUT_PDF)

$(OUTPUT_PDF): $(MAIN_TEX) $(wildcard $(CHAPTERS_DIR)/*.tex) $(wildcard $(FIGURES_DIR)/*)
	@mkdir -p $(BUILD_DIR)
	@echo "=== Building LaTeX document (pass 1) ==="
	$(LATEX) $(LATEX_OPTS) -output-directory=$(BUILD_DIR) $(MAIN_TEX)
	@echo "=== Building LaTeX document (pass 2) ==="
	$(LATEX) $(LATEX_OPTS) -output-directory=$(BUILD_DIR) $(MAIN_TEX)
	@echo "=== Building LaTeX document (pass 3) ==="
	$(LATEX) $(LATEX_OPTS) -output-directory=$(BUILD_DIR) $(MAIN_TEX)
	@cp $(BUILD_DIR)/$(OUTPUT_PDF) .
	@echo "=== PDF built: $(OUTPUT_PDF) ==="

# Quick build (single pass)
pdf-quick:
	$(LATEX) $(LATEX_OPTS) -output-directory=$(BUILD_DIR) $(MAIN_TEX)
	cp $(BUILD_DIR)/$(OUTPUT_PDF) .

# View PDF
view: pdf
	@if command -v open >/dev/null 2>&1; then open $(OUTPUT_PDF); \
	elif command -v xdg-open >/dev/null 2>&1; then xdg-open $(OUTPUT_PDF); \
	elif command -v evince >/dev/null 2>&1; then evince $(OUTPUT_PDF); \
	else echo "Please open $(OUTPUT_PDF) manually"; fi

# Run Python demos
demo:
	@echo "=== Running Python demos ==="
	$(PYTHON) -c "from euler import *; import math; \
	print('Totient:', [euler_totient(n) for n in range(1, 13)]); \
	print('Basel:', basel_problem()); \
	print('Euler identity:', euler_identity()); \
	print('Partitions:', [partition_function(n) for n in range(1, 11)])"

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR)
	rm -f *.aux *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk
	rm -f *.synctex.gz *.nav *.snm *.vrb
	rm -f $(OUTPUT_PDF)
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete
	rm -rf .pytest_cache .mypy_cache .coverage

# Clean everything including installed package
clean-all: clean
	pip uninstall -y euler-mathematics 2>/dev/null || true

# Check Python code style
lint:
	@echo "Checking Python style..."
	@$(PYTHON) -m flake8 $(SRC_DIR) $(TEST_DIR) --ignore=E501,W503 2>/dev/null || echo "flake8 not installed"
	@$(PYTHON) -m mypy $(SRC_DIR)/euler --ignore-missing-imports 2>/dev/null || echo "mypy not installed"

# Format Python code
format:
	@echo "Formatting Python code..."
	@$(PYTHON) -m black $(SRC_DIR) $(TEST_DIR) 2>/dev/null || echo "black not installed"

# Generate figures
figures:
	@echo "Generating figures..."
	@mkdir -p $(FIGURES_DIR)
	$(PYTHON) scripts/generate_figures.py 2>/dev/null || echo "No figure generation script"

# Help
help:
	@echo "Euler Student Guide - Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all          - Build PDF and run tests (default)"
	@echo "  pdf          - Build LaTeX document (3 passes)"
	@echo "  pdf-quick    - Quick LaTeX build (1 pass)"
	@echo "  test         - Run all Python tests"
	@echo "  test-cov     - Run tests with coverage"
	@echo "  test-module  - Run specific test module (MODULE=name)"
	@echo "  view         - Build and open PDF"
	@echo "  demo         - Run quick Python demos"
	@echo "  install      - Install package in development mode"
	@echo "  clean        - Remove build artifacts"
	@echo "  clean-all    - Clean and uninstall package"
	@echo "  lint         - Check Python code style"
	@echo "  format       - Format Python code with black"
	@echo "  figures      - Generate figures"
	@echo "  help         - Show this help"
	@echo ""
	@echo "Examples:"
	@echo "  make pdf"
	@echo "  make test"
	@echo "  make view"
	@echo "  make test-module MODULE=test_euler"
	@echo "  make demo"