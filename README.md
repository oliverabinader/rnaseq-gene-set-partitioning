# RNA-seq Gene Set Partitioning and Downstream Analysis

**Author:** Oliver Abinader  

## Overview

This repository provides workflows for partitioning differential gene expression (DEG) results based on biologically curated gene sets (e.g., Reactome pathways), with a focus on mitochondrial gene classification.

The pipeline enables:

- Separation of DEGs into biologically meaningful groups
- Comparison across multiple datasets or cell lines


## Analyses Included

### 1. Gene Set-Based Partitioning

- Uses curated gene sets (e.g., Reactome pathways)
- Extracts genes belonging to the selected pathways (gene sets of interest)
- Splits DEGs into:
  - Target gene set (e.g., mitochondrial genes)
  - Background / remaining genes

### 2. Multi-Dataset DEG Splitting

- Applies gene set partitioning across multiple datasets
- Supports:
  - Upregulated genes
  - Downregulated genes

### 3. Annotation-Based Gene Classification

- Assigns each gene into:
  - `Mitochondria`
  - `Remaining`
- Adds annotation directly to DEG tables


## Input Requirements

### Differential Expression Data
- Must contain:
  - `GeneSymbol`
  - `log2FoldChange` (optional for filtering)
  - `padj` or `FDR`

### Gene Set File (GMT format)
- Standard `.gmt` file (e.g., MSigDB Reactome)
- OR a simple gene list
