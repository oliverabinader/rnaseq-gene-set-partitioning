# RNA-seq Gene Set Partitioning and Downstream Analysis

**Author:** Oliver Abinader  

## Overview

This repository provides workflows for partitioning differential gene expression (DEG) results based on biologically curated gene sets (e.g., Reactome pathways).

The pipeline enables:

- Separation of DEGs into biologically meaningful groups
- Comparison across multiple datasets or cell lines

---

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
 
Filtering Criteria:

By default, the script applies:
- |log2FoldChange| > 1
- FDR (padj) < 0.05
The thresholds can be modified in the script itself

---

## Input Requirements

### Differential Expression Data
- Must contain:
  - `GeneSymbol`
  - `log2FoldChange` (used for filtering)
  - `padj` or `FDR`

### Gene Set File (GMT format)
- Standard `.gmt` file (e.g., MSigDB Reactome)
- OR a simple gene list

---

## Output Files

Output will be saved as Excel files containing:

- Upregulated genes (gene set vs background)
- Downregulated genes (gene set vs background)
  
Each category is saved as a separate sheet.

---

## Notes

- Input datasets are imported by the user
- Update file paths before running
- GMT files can be downloaded from MSigDB or replaced with user-defined gene sets
