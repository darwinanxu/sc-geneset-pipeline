# sc-gene-set-pipeline

A modular and scalable Python pipeline for functional gene-set scoring and interpretable cell-state annotation in single-cell RNA-seq data.

## Overview

This project is designed as an engineering-focused toolkit for benchmarking and applying gene-set scoring methods in single-cell datasets. Instead of building a one-off notebook analysis, the goal is to provide a reusable pipeline with standardized preprocessing, scoring, evaluation, and reporting.

## Features

- Modular scoring interface
- Support for multiple gene-set scoring methods
- Standardized preprocessing for single-cell data
- Evaluation of method performance and QC confounding
- Script-based and extensible project structure

## Project Structure

```text
sc-gene-set-pipeline/
├── README.md
├── pyproject.toml
├── .gitignore
├── configs/
├── data/
├── notebooks/
├── results/
├── scripts/
├── src/
│   └── sc_gene_set_pipeline/
├── tests/
└── examples/
