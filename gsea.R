library(geneset)
library(genekitr)
library(tidyverse)
library(dplyr)
library(readr)
# gskegg <- getKEGG(org = "mmu",category = "pathway")
# gsmsigdb <- getMsigdb(org = "mouse", category = "C5-HPO")
# gswiki <- getWiki(org = "mouse")
# gsenrich <- getEnrichrdb(org = "human", library = "Cancer_Cell_Line_Encyclopedia")
# enrichr_metadata_frame <- enrichr_metadata


# Function to handle multiple gene sets
process_multiple_gene_sets <- function(base_dir, exp_names, gs, type, p_cutoff = 0.05, q_cutoff = 0.05)
{
  for (exp_name in exp_names)
  {
    if (type == "go")
    {
      go_terms_revigo_file <- sprintf("%s/gene_centered_results/go_lists/Revigo_goterms_%s.tsv",
                                base_dir, exp_name)
    }
    else
    {
      go_terms_revigo_file <- ""
    }
      dge_results_file <- sprintf("%s/%s/differential_abundance/tables/differential/condition_%s.deseq2.results_filtered.tsv",
                                base_dir, exp_name, exp_name)
    process_gene_set(dge_results_file, gs, p_cutoff, q_cutoff,exp_name, base_dir, go_terms_revigo_file, type)
  }
}


process_gene_set <- function(dge_results_file, gs, p_cutoff = 0.05, q_cutoff = 0.05, exp_name, base_dir, go_terms_revigo_file, type)
{
  # Read the data
  dge_df <- read_tsv(dge_results_file, show_col_types = FALSE)

  # Subset and rename columns
  dge_df <- subset(dge_df, select = -c(baseMean, lfcSE, pvalue, padj))
  colnames(dge_df) <- c("geneID_symbol", "logfc")
  dge_df <- dge_df %>% arrange(desc(logfc))
  genes_of_interest <- setNames(dge_df$logfc, dge_df$geneID_symbol)
  x <- typeof(genes_of_interest)
  x
  head(genes_of_interest)
  # colnames(genes_of_interest) <- c("logfc", "geneID_symbol")
  gsego <- genGSEA(genelist = genes_of_interest, geneset = gs, p_cutoff = p_cutoff, q_cutoff = q_cutoff)
  if (type == "go")
  {
    gsego <- filter_gsea_by_go_terms(gsego, go_terms_revigo_file)
  }
  # x <- typeof(gsego$gsea_df)
  # y <- typeof(genes_of_interest)
  # plotGSEA(gsego, plot_type = "bar", colour = c("blue", "red"))
  # ggsave(filename = paste0(sprintf("%s/gene_centered_results/go_figures/%s_gsea_plot.png",base_dir,exp_name)),
  #        plot = last_plot(), width = 8, height = 6, dpi = 600)
  # gene_list <- c("Apc","Bub1b","Cdx2","Elavl1","Fpr2","Ptprt","Tlr2","Trp53inp1", "Kras", "Braf", "Akt1", "Src")
  # gsego$gsea_df$geneID_symbol <- gsego$gsea_df$geneID
  plotGSEA(gsego, plot_type = "bar")
  # plotEnrich(gsego$gsea_df, plot_type = "geneheat", fold_change = genes_of_interest)
  ggsave(filename = paste0(sprintf("%s/gene_centered_results/%s_figures/%s_gsea_plot2.png",base_dir,type,exp_name)),
         plot = last_plot())
}


filter_gsea_by_go_terms <- function(gsego, go_terms_file) {
  # Read the GO terms from the external file
  go_terms_to_keep <- read_tsv(go_terms_file, show_col_types = FALSE)$TermID

  # Filter the gsea_df table in the gsego object
  gsego$gsea_df <- gsego$gsea_df %>%
    filter(mmusculus_BP_ID %in% go_terms_to_keep)

  return(gsego)
}

gsgo <- getGO(org = "mouse",ont = "bp")
gskegg <- getKEGG(org = "mmu",category = "pathway")
base_dir <- "/home/direnc/results/tim/rnaseq_mice"
exp_names <- c("nr_r", "nr_nrt", "r_rt", "nrt_rt")
type <- "kegg"
process_multiple_gene_sets(base_dir, exp_names, gskegg, type)