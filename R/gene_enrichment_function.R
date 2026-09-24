
# Extract top genes that contribute to different entropic groups

selectGenes <- function(entropy_vals,ngenes){
  
  topGenes <- vector("list",length(entropy_vals))
  for (i in seq_along(entropy_vals)){
    gbm <- entropy_vals[[i]]
    group <- sort(unique(gbm$group))
    gs <- GetAssayData(gbm, assay="Entropy", layer = "counts")
    
    genes <- vector("list",length(group))
    names(genes) <- group
    for (g in group){
      spots <- rownames(gbm@meta.data[gbm@meta.data$group == g,])
      gs.matrix <- gs[,spots,drop=FALSE]
      mean_geneEnt <- rowMeans(gs.matrix)
      mean_geneEnt <-sort(mean_geneEnt[mean_geneEnt != 0],decreasing=TRUE) 
      genes[[g]] <- names(head(mean_geneEnt,ngenes))
      
    }
    
    topGenes[[i]] <- genes
  }
  return(topGenes)
}

## Gene Ontology enrichment 

geneEnrich <- function(topgenes,entropy_vals){
  results <- vector("list",length(topgenes))
  
  for(i in seq_along(topgenes)){
    cat("Go enrichment in process for ",names(entropy_vals)[i], "\n")
    gs <- GetAssayData(entropy_vals[[i]],assay="Entropy",layer="counts")
    bg_genes <- rownames(gs)
    sample <- topgenes[[i]]
    groups <- names(sample)
    
    results[[i]] <- compareCluster(geneClusters  = sample , fun="enrichGO",
                                   OrgDb=org.Hs.eg.db, keyType = "SYMBOL",ont = "BP",
                                   pvalueCutoff = 0.05,pAdjustMethod = "BH",
                                   universe=bg_genes)
  }
  
  return(results)
}

## KEGG Pathway enrichment

keggEnrich <- function(topgenes){
  results <- vector("list",length(topgenes))
  
  for(i in seq_along(topgenes)){
    sample <- topgenes[[i]]
    entrez <- lapply(sample,function(genes){
      res <- bitr(genes, fromType="SYMBOL",toType="ENTREZID",
                  OrgDb=org.Hs.eg.db)$ENTREZID
    })
    results[[i]] <- compareCluster(geneClusters = entrez,
                                   fun = "enrichKEGG",organism="hsa",pvalueCutoff=0.05)
  }
  return(results)
}