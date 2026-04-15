# ============================================
# Split DEG Results into Gene Set vs Background
# ============================================

# Objective:
# Divide DEGs into:
#   1) Genes belonging to selected pathways (gene sets of interest)
#   2) Remaining background genes

# ============================================
# Load libraries
# ============================================

library(readxl)
library(openxlsx)
library(dplyr)

# ============================================
# User Inputs
# ============================================

# Input DEG files
kuramochi_file <- "/path/to/DE_results.xlsx"

# GMT file
gmt_file <- "/path/to/reactome.v2024.1.Hs.symbols.gmt"

# Output directory
output_dir <- "/path/to/results/"

# Define pathways of interest
reactome_paths <- c(
  "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
  "REACTOME_MITOCHONDRIAL_TRANSLATION_ELONGATION",
  "REACTOME_MITOCHONDRIAL_TRNA_AMINOACYLATION",
  "REACTOME_MITOCHONDRIAL_UNCOUPLING",
  "REACTOME_CELLULAR_RESPONSE_TO_MITOCHONDRIAL_STRESS",
  "REACTOME_MALATE_ASPARTATE_SHUTTLE",
  "REACTOME_CARNITINE_SHUTTLE",
  "REACTOME_TRANSPORT_OF_FATTY_ACIDS",
  "REACTOME_CYTOSOLIC_IRON_SULFUR_CLUSTER_ASSEMBLY",
  "REACTOME_MITOCHONDRIAL_MRNA_MODIFICATION",
  "REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION",
  "REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT",
  "REACTOME_MITOCHONDRIAL_RNA_DEGRADATION"
)

# ============================================
# Load Data
# ============================================

kuramochi_all <- read_excel(kuramochi_file)
#kuramochi_all$GeneSymbol <- trimws(kuramochi_all$GeneSymbol) # We may need to do it to clean spsace in the gene symbol column if present in the dataframe.

# ============================================
# Filter DEGs
# ============================================

# Performing a 2-FC with an FDR<5% on DEGs.
filter_deg <- function(df, fc = 1, fdr = 0.05) {  
  df_up <- df %>% filter(log2FoldChange > fc, padj < fdr)
  df_dn <- df %>% filter(log2FoldChange < -fc, padj < fdr)
  list(up = df_up, dn = df_dn)
}

kur <- filter_deg(kuramochi_all)

# ============================================
# Load GMT Gene Sets
# ============================================

gmt_lines <- readLines(gmt_file)

reactome_list <- lapply(gmt_lines, function(x) strsplit(x, "\t")[[1]])
names(reactome_list) <- sapply(reactome_list, `[`, 1)

reactome_list <- lapply(reactome_list, function(x) x[-c(1,2)])

reactome_subset <- reactome_list[reactome_paths]

all_genes <- unique(unlist(reactome_subset))

# ============================================
# Split Function
# ============================================

split_genes <- function(df, gene_list) {
  list(
    geneset = df %>% filter(GeneSymbol %in% gene_list),
    background = df %>% filter(!GeneSymbol %in% gene_list)
  )
}

kur_up  <- split_genes(kur$up, all_genes)
kur_dn  <- split_genes(kur$dn, all_genes)

# ============================================
# Save Outputs
# ============================================

save_split <- function(up, dn, prefix) {

  wb <- createWorkbook()

  addWorksheet(wb, "Up_Geneset")
  addWorksheet(wb, "Up_Background")
  addWorksheet(wb, "Dn_Geneset")
  addWorksheet(wb, "Dn_Background")

  writeData(wb, "Up_Geneset", up$geneset)
  writeData(wb, "Up_Background", up$background)
  writeData(wb, "Dn_Geneset", dn$geneset)
  writeData(wb, "Dn_Background", dn$background)

  saveWorkbook(wb, file.path(output_dir, paste0(prefix, "_DEG_partitioned_by_geneset.xlsx")), overwrite = TRUE)
}

save_split(kur_up, kur_dn, "Kuramochi")

#Note:
If you just want to split the DE genes into Mito and Background, do the following in the code:
1. Select only Mito pathways 
Mito_pathways <- reactome_paths[str_detect(reactome_paths, "MITOCHONDRIAL|MITOCHONDRION")]
2. When you do this, tweak it a bit: 
reactome_subset <- reactome_list[Mito_pathways]
3. When you read in the DEG data, you can add a new col with either Mitochondria or Remaining. In this analysis, you don't really care about splitting DEGs into up and down lists.
kuramochi_all <- dplyr::kuramochi_all %>% mutate(GeneGroup = ifelse(GeneSymbol %in% all_genes, "Mitochondria", "Background"))
data_to_save <- kuramochi_all %>% rownames_to_column(var = "GeneID")
data_to_save <- data_to_save %>% mutate(PValue = as.numeric(PValue), FDR = as.numeric(FDR)) #just to make these columns numeric before writing to output file.
4. Write the results to an output file
writexl::write_xlsx(
  x = data_to_save, 
  path = "path/Kuramochi_DE_Mitochondria_vs_Remaining.xlsx",
  col_names = TRUE
)
