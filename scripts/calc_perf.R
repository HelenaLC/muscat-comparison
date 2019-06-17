suppressPackageStartupMessages({
    library(dplyr)
    library(iCOBRA)
    library(SingleCellExperiment)
})

sce <- readRDS(snakemake@input$sce)
res <- readRDS(snakemake@input$res)$res

gis <- metadata(sce)$gene_info %>% 
    mutate_at("cluster_id", as.character) %>% 
    mutate(is_de = as.integer(!category %in% c("ee", "ep")))
df <- left_join(gis, res, by = c("cluster_id", "gene"))

pval <- data.frame(df$p_val)
colnames(pval) <- snakemake@wildcards$method_id

cd <- COBRAData(pval, 
    padj = df[, c("p_adj.loc", "p_adj.glb")],
    truth = data.frame(is_de = df$is_de))

prf <- suppressMessages(calculate_performance(cd, binary_truth = "is_de"))

saveRDS(prf, snakemake@output$prf)