
#Load packages
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(here)

#Load functions
source(here("R","gene_enrichment_function.R"))

#Load results
input_file <- here("data", "processed", "ranking_and_classification_greenwarld.rds")
entropy_vals <- readRDS(input_file)

#select the highest entropy associated genes
n_top_genes <- 100 #number of top genese
sel_genes <- selectGenes(entropy_vals=entropy_vals,ngenes=n_top_genes)

#go enrichment
go_results <- geneEnrich(sel_genes,entropy_vals)
go_output <- here("data","processed","go_enrichment.rds")
saveRDS(go_results,go_output)

go_plot_dir <- here("results","figures","gene_enrichment","GO")
dir.create(go_plot_dir, recursive=TRUE, showWarnings = FALSE)

go_plots <- vector("list",length(go_results))
names(go_plots) <- names(entropy_vals)

for (i in seq_along(go_results)){
  p <-dotplot(go_results[[i]], showCategory = 10) +
    ggtitle(names(entropy_vals)[i]) +
    theme(axis.text.y = element_text(size = 4))
  go_plots[[i]] <- p
  output_file <- file.path(go_plot_dir, paste0("GO_", names(entropy_vals)[i],".png"))
  ggsave(filename = output_file, plot = p, width = 8, height = 6,dpi = 300)
}  
  
#KEGG pathway enrichment
kegg_results <- keggEnrich(topgenes=sel_genes)
kegg_output <- here("data","processed","KEGG_enrichment.rds")
saveRDS(kegg_results,kegg_output)

kegg_plot_dir <- here("results","figures","gene_enrichment", "KEGG")
dir.create(kegg_plot_dir, recursive = TRUE, showWarnings = FALSE)
kegg_plots <- vector("list",length(kegg_results))
names(kegg_plots) <- names(entropy_vals)

for (i in seq_along(kegg_results)) {
  f <- dotplot(kegg_results[[i]], showCategory = 15) +
    ggtitle(names(entropy_vals)[i]) + theme(axis.text.y = element_text(size=6.5, lineheight=0.8))
kegg_plots[[i]] <- f
output_file <- file.path(kegg_plot_dir,paste0("KEGG_",names(entropy_vals)[i],".png"))
  ggsave(filename = output_file,plot = f, width = 8,height = 6,dpi = 300)
}
