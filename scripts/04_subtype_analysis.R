
#Load packages
library(Seurat)
library(ggplot2)
library(matrixStats)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(here)
library(dplyr)

#Load samples
source(here("R","subtype_overlap_function.R"))
input_file <- here("data","processed","ranking_and_classification.rds")
entropy_vals <- readRDS(input_file)
sample_names <- names(entropy_vals)

#Calculate gene level entropy variability
std <- sdGenes(entropy_vals)

#Genes in the top 10% of entropy variability
genes <- vector("list",length(entropy_vals))
names(genes) <- sample_names

for(i in seq_along(entropy_vals)){
  sample_data <-std[std$Sample==sample_names[i],]
  sample_data <- sample_data[order(-sample_data$SD),] #sort SD in descending order
  threshold <- quantile(sample_data$SD,0.90,na.rm=TRUE)
  genes[[i]] <- sample_data$genes[sample_data$SD > threshold]
}

#Number of high SD genes per sample
for (i in seq_along(entropy_vals)){
  cat(sample_names[i], "high sd genes", length(genes[[i]]), "\n")
}

#Select top 15 genes across all samples
all_genes <- unique(unlist(lapply(genes,head,15)))
sd_gene_matrix <- matrix(NA, nrow=length(all_genes),ncol=length(entropy_vals),dimnames = list(all_genes,sample_names))

for (i in seq_along(entropy_vals)){
  sample_data <- std[std$Sample==sample_names[i],]
  rownames(sample_data) <- sample_data$genes
  matched <- sample_data[all_genes,"SD",drop=TRUE]
  sd_gene_matrix[,i]<-matched
}

#sd heatmap
figure_dir <- here("results","figures","subtype_analysis")
dir.create(figure_dir, recursive = TRUE,showWarnings = FALSE)
png(filename = file.path(figure_dir,"sd_top_genes_heatmap.png"), width=1200,
                         height=1600, res=150)
pheatmap(sd_gene_matrix,cluster_rows = TRUE,cluster_cols = FALSE,show_rownames = TRUE,fontsize_row = 10,fontsize_col=8,cellwidth = 14,cellheight = 10,border_color=NA,main="Standard deviation of Entropy per gene across all samples",color = colorRampPalette(c("white","#4575B4","#D73027"))(100))
dev.off()

#Genes that are shared across all samples
shared_all <- Reduce(intersect,genes)

#Genes that are unique to each sample
unique_genes <- vector("list",length(entropy_vals))
for(i in seq_along(entropy_vals)){
  not_shared <- unlist(genes[names(genes)!=sample_names[i]])
  unique_genes[[i]] <- setdiff(genes[[i]],not_shared)
}

unique_top <- vector("list",length(entropy_vals))
names(unique_top) <- sample_names

for (i in seq_along(entropy_vals)){
  sample_sd <- std[std$Sample == sample_names[i] & std$genes %in% unique_genes[[i]],]
  sample_sd <- sample_sd[order(-sample_sd$SD),]
  unique_top[[i]] <- head(sample_sd$genes,100)
}

#partially shared genes among samples
gene_sample_count <- table(unlist(genes))
partially_shared_genes <- names(gene_sample_count[gene_sample_count>= 2 & 
                                                    gene_sample_count < length(genes)])
#summary statistics
summary_df <- data.frame(
  sample        = sample_names,
  total         = sapply(genes, length),
  unique        = sapply(unique_genes, length),
  shared_all    = sapply(genes, function(x) length(intersect(x,shared_all))),
  partially_shared = sapply(genes,function(x){
    length(intersect(x,partially_shared_genes))
  }),
  stringsAsFactors = FALSE
)
summary_df

#save summary statistics
write.csv(summary_df,file=here("results","tables","subtype_analysis_summary.csv"),
          row.names = FALSE)

#GO enrichment for shared genes among samples
if(length(shared_all)>0){
  shared_go <- enrichGO(gene = shared_all,OrgDb=org.Hs.eg.db, keyType = "SYMBOL",ont = "ALL",
                        pvalueCutoff = 0.05,pAdjustMethod = "BH")
  p_shared_go <- dotplot(shared_go,showCategory=10) + ggtitle("Shared SD genes across samples - GO enrichment") +
    theme(axis.text.y = element_text(size=6.5, lineheight=0.8))
  p_shared_go
  ggsave(filename=file.path(figure_dir,"shared_sd_go.png"),plot=p_shared_go,width=8, height=6,dpi=300)
  
} else {
  message("No genes were shared across all samples")
}

#Go enrichment for unique genes among samples
gbm_samples<- setNames(unique_genes,sample_names)
gbm_samples <- gbm_samples[sapply(gbm_samples,length)>0]

if(length(gbm_samples)>0){
  u.go.results <- compareCluster(geneClusters  = gbm_samples , fun="enrichGO",
                                 OrgDb=org.Hs.eg.db, keyType = "SYMBOL",ont = "BP",
                                 pvalueCutoff = 0.05,pAdjustMethod = "BH")
  
  p <- dotplot(u.go.results, showCategory = 10) + ggtitle("Unique SD genes among samples - GO enrichment") +
    theme(axis.text.y = element_text(size=3 , lineheight=0.8),
          axis.text.x = element_text(angle=45, hjust=1, size=9))
  p
  ggsave(filename=file.path(figure_dir,"unique_sd_go.png"),plot=p,width = 10,height=12,units="in",dpi=300)
  
} else {
  message("No sample-specific genes were identified")
}

length(shared_all)
head(shared_all)

#KEGG Pathway enrichment for shared genes among samples
if(length(shared_all)>0){
  mapped_genes <- bitr(shared_all, fromType="SYMBOL",toType="ENTREZID", OrgDb=org.Hs.eg.db)
  
  #remove duplicated and missing mappings
  mapped_genes <- mapped_genes[!is.na(mapped_genes$ENTREZID),]
  mapped <- mapped_genes[!duplicated(mapped_genes$ENTREZID),]
  shared_entrez <- unique(mapped_genes$ENTREZID)
  cat("Shared genes:", length(shared_all), "\n",
    "Mapped to Entrez IDs:", length(shared_entrez), "\n")
  
  #KEGG enrichment
  shared_kegg<- enrichKEGG(shared_entrez,organism = "hsa",pvalueCutoff=0.05, pAdjustMethod = "BH")
  
  if(!is.null(shared_kegg) && nrow(as.data.frame(shared_kegg))>0){
    p_shared_kegg <- dotplot(shared_kegg,showCategory=10) + ggtitle("Shared SD genes across samples - Pathway enrichment")
    p_shared_kegg
    ggsave(filename=file.path(figure_dir,"shared_sd_kegg.png"),plot=p_shared_kegg,width=8,height=6,dpi=300)
    
  } else {
    message("No significantly enriched KEGG pathways were found for shared genes.")
  }
  
}else{
  message("No shared genes available")
}

#KEGG Pathway enrichment for unique genes among samples
u.entrez <- lapply(unique_genes,function(gene){
  if (length(gene)==0) {
    return (NULL)
  }
  #covert SYMBOL to ENTREZID
  mapped <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID", OrgDb=org.Hs.eg.db)
  
  #Remove missing and duplicate IDS
  mapped <- mapped[!is.na(mapped$ENTREZID), ]
  mapped <- mapped[!duplicated(mapped$ENTREZID), ]
  unique(mapped$ENTREZID)
})

gbm_samp <- setNames(u.entrez, sample_names)
gbm_samp<- gbm_samp[sapply(gbm_samp,length)>0]
cat("Samples with KEGG-compatible unique genes:",length(gbm_samp),"\n")

if(length(gbm_samp)>0){
  u.results <- compareCluster(geneClusters  = gbm_samp , fun="enrichKEGG",
                              organism="hsa", pvalueCutoff = 0.05,pAdjustMethod = "BH")
  
  #check if any significant pathway is found
  u.results_df <- as.data.frame((u.results))
  
  if (nrow(u.results_df)>0){
    p_unique_kegg <- dotplot(u.results, showCategory = 10) + ggtitle("Unique SD genes across samples - KEGG Pathway enrichment") +
      theme(axis.text.y = element_text(size=6.5, lineheight=0.8),axis.text.x = element_text(angle=45, hjust=1, size=9))
    p_unique_kegg
    ggsave(filename = file.path(figure_dir,"unique_sd_kegg.png"),plot=p_unique_kegg,width=19,height=12,dpi=300)
  } else {
    message("No significantly enriched KEGG pathways were found for the unique genes.")
  }
  
} else {
  
  message("No samples contained genes that could be mapped to KEGG.")
}

saveRDS(list(sd_genes = std, high_sd_genes = genes, shared_genes = shared_all, 
             unique_genes = unique_genes, unique_top_genes = unique_top,
             partially_shared_genes = partially_shared_genes, summary = summary_df, 
             sd_gene_matrix = sd_gene_matrix ), file = here( "data", "processed", "subtype_analysis.rds" ) )

# Sensitivity analysis excluding mitochondrial and ribosomal genes
shared_mt <- shared_all[grepl("^MT-", shared_all)]
shared_ribo <- shared_all[grepl("^RPS|^RPL", shared_all)]
data.frame(category = c("All shared genes","Mitochondrial genes","Ribosomal genes"),
           n = c(length(shared_all),length(shared_mt),length(shared_ribo)),
           percentage = c(100, 100 * length(shared_mt) / length(shared_all), 100 * length(shared_ribo) / length(shared_all)))

shared_filtered <- shared_all[!grepl("^MT-", shared_all) &!grepl("^RPS|^RPL", shared_all)]
cat("Shared genes before filtering:", length(shared_all), "\n",
  "Mitochondrial genes removed:", length(shared_mt), "\n",
  "Ribosomal genes removed:", length(shared_ribo), "\n",
  "Shared genes retained:", length(shared_filtered), "\n")

#Filtered GO enrichment

shared_filtered_go <- enrichGO(gene = shared_filtered,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",
                               ont = "ALL",pvalueCutoff = 0.05,pAdjustMethod = "BH")

if (!is.null(shared_filtered_go) && nrow(as.data.frame(shared_filtered_go)) > 0) {
  p_shared_filtered_go <- dotplot(shared_filtered_go,showCategory = 10) +
    ggtitle("Shared high-variability genes excluding mitochondrial and ribosomal genes") +
    theme(axis.text.y = element_text(size = 6.5))
  p_shared_filtered_go
  ggsave(filename = file.path(figure_dir,"shared_sd_go_filtered.png"),
    plot = p_shared_filtered_go,width = 8, height = 6, dpi = 300)
  
} else {
  message("No significant GO enrichment after excluding mitochondrial and ribosomal genes.")
}

#Filtered KEGG enrichment

mapped_filtered <- bitr(shared_filtered,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
mapped_filtered <- mapped_filtered[!is.na(mapped_filtered$ENTREZID),]
mapped_filtered <- mapped_filtered[!duplicated(mapped_filtered$ENTREZID),]
shared_filtered_entrez <- unique(mapped_filtered$ENTREZID)
cat("Filtered shared genes:", length(shared_filtered), "\n","Mapped Entrez IDs:", length(shared_filtered_entrez), "\n")

# KEGG enrichment after excluding mitochondrial and ribosomal genes

shared_filtered_kegg <- enrichKEGG(gene = shared_filtered_entrez,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")
if (!is.null(shared_filtered_kegg) && nrow(as.data.frame(shared_filtered_kegg)) > 0) {
  p_shared_filtered_kegg <- dotplot(shared_filtered_kegg,showCategory = 10) +
    ggtitle("Shared high-variability genes excluding mitochondrial and ribosomal genes - KEGG") +
    theme(axis.text.y = element_text(size = 7))
  p_shared_filtered_kegg
  ggsave(filename = file.path(figure_dir,"shared_sd_kegg_filtered.png"),
         plot = p_shared_filtered_kegg, width = 8, height = 6, dpi = 300)
} else {
  message("No significant KEGG enrichment after excluding mitochondrial and ribosomal genes.")
  
}

#filtered genes summary
shared_filter_summary <- data.frame(category = c("All shared high SD genes", "Mitochondrial genes","Ribosomal genes","filtered shared genes"),
                                    n = c(length(shared_all),length(shared_mt),length(shared_ribo),length(shared_filtered)))
shared_filter_summary$percentage <- 100 * shared_filter_summary$n / length(shared_all)
shared_filter_summary

length(setdiff(shared_all, shared_filtered))
length(intersect(shared_mt, shared_ribo))
head(shared_filtered, 50) #filtered genes

shared_filtered_mt <- shared_filtered[grepl("MT-",shared_filtered)]
shared_filtered_ribo <- shared_filtered[grepl("^RPS|^RPL",shared_filtered)]
length(shared_filtered_mt)
length(shared_filtered_ribo)

shared_filtered_housekeeping <- shared_filtered[grepl("^MT-|^RPS|^HIST",shared_filtered)]
length(shared_filtered_housekeeping)

as.data.frame(shared_filtered_go) |> 
  dplyr::select(ID,Description, GeneRatio,BgRatio,p.adjust) |>
  head (20)

as.data.frame(shared_filtered_kegg) |> 
  dplyr::select(ID,Description,GeneRatio,BgRatio,p.adjust) |>
  head(20)

#Filter mitochondrial and ribosomal genes from sample specific genes
unique_filtered <- lapply(unique_genes, function(gene_set) {
  gene_set[!grepl("^MT-|^RPS|^RPL", gene_set)]
})
names(unique_filtered) <- sample_names
unique_filtered_summary <- data.frame(sample=sample_names, original_unique=sapply(unique_genes,length),
                                      filtered_unique <- sapply(unique_filtered,length))
unique_filtered_summary$removed <- unique_filtered_summary$original_unique - unique_filtered_summary$filtered_unique
unique_filtered_summary
unique_filtered_summary[order(unique_filtered_summary$filtered_unique, decreasing = TRUE),] #samples with enough genes
sum(sapply(unique_filtered, length) > 0)
lapply(unique_filtered, head, 20)

# GO enrichment of filtered unique genes
unique_go_input <- unique_filtered[sapply(unique_filtered, length) >= 10]
unique_go_input <- lapply(unique_go_input, unique)
unique_go_results <- compareCluster(geneClusters = unique_go_input, fun = "enrichGO",OrgDb = org.Hs.eg.db,
                                    keyType = "SYMBOL",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH")

# Check whether enrichment was detected
unique_go_df <- as.data.frame(unique_go_results)
nrow(unique_go_df)
head(unique_go_df[, c("Cluster","Description","GeneRatio","BgRatio","p.adjust")],20)

unique_go_df |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::count(Cluster, name = "n_significant") |>
  dplyr::arrange(desc(n_significant))

unique_go_plot <- dotplot(unique_go_results,showCategory = 10) +
  ggtitle("Sample-specific high-SD genes: GO Biological Process") +
  theme(axis.text.y = element_text(size = 7),axis.text.x = element_text(angle = 45,hjust = 1),
        plot.title = element_text(face = "bold",hjust = 0.5))
unique_go_plot
ggsave(filename = file.path(figure_dir,"unique_sd_go_filtered.png"),
       plot = unique_go_plot,width = 12, height = 10,dpi = 300)
