
#Calculate gene level variability across spots per sample

sdGenes <- function(entropy_vals){
  sd_list <- vector("list",length(entropy_vals))
  for (i in seq_along(entropy_vals)){
    gs <- GetAssayData(entropy_vals[[i]],assay="Entropy",layer="counts")
    sd_list[[i]] <- data.frame(Sample=names(entropy_vals)[i],genes = rownames(gs), 
                               SD = matrixStats::rowSds(as.matrix(gs)),
                               stringsAsFactors = FALSE )
    
  }
  sd_genes <- dplyr::bind_rows(sd_list)
  return(sd_genes)
}

#correlation calculation
calculate_signature_correlations <- function(entropy_vals) {
  
  sample_names <- names(entropy_vals)
  results <- lapply(sample_names,function(sample_name) {
    gbm <- entropy_vals[[sample_name]]
    data.frame(sample = sample_name,GSC = cor(gbm$mean_entropy,gbm$GSC1,use = "complete.obs",method="pearson"),
               MES = cor(gbm$mean_entropy,gbm$SubType1, use = "complete.obs",method="pearson"),
               PN = cor(gbm$mean_entropy,gbm$SubType2,use = "complete.obs",method="pearson"),
               CL = cor(gbm$mean_entropy,gbm$SubType3, use = "complete.obs",method="pearson"),
               NL = cor(gbm$mean_entropy,gbm$SubType4, use = "complete.obs",method="pearson"))
              }
  )
  dplyr::bind_rows(results)
}

#fishers exact test 
calculate_fisher_overlap <- function(entropy_vals,quantile_cutoff) {
  
  #validate quantile cutoff
  if(quantile_cutoff <= 0 || quantile_cutoff >= 1){
    stop ("quantile_cutoff must be between 0 and 1.")
  }
  
  sample_names <- names(entropy_vals)
  results <- lapply(sample_names, function(sample_name) {
    gbm <- entropy_vals[[sample_name]]
    entropy <- gbm$mean_entropy
    gsc <- gbm$GSC1
    
    # Remove spots with missing values
    valid <- complete.cases(entropy, gsc)
    entropy <- entropy[valid]
    gsc <- gsc[valid]
    
    #test for enough observations for meaningful test
    if(length(entropy)<2){
      return(data.frame(sample=sample_name,observed=NA,expected=NA,ratio = NA,pvalue=NA))
    }
    
    # Calculate sample-specific thresholds
    ent_threshold <- quantile(entropy, quantile_cutoff, na.rm = TRUE)
    gsc_threshold <- quantile(gsc, quantile_cutoff, na.rm = TRUE)
    
    # Define high/low groups
    high_entropy <- entropy > ent_threshold
    high_gsc <- gsc > gsc_threshold
    
    # Contingency table
    contingency_matrix <- table(high_entropy,high_gsc)
    
    # Fisher's exact test
    fisher_result <- fisher.test(contingency_matrix)
    
    # Observed overlap
    observed <- sum( high_entropy & high_gsc)
    
    # Expected overlap under independence
    total <- length(high_entropy)
    expected <- total * mean(high_entropy) * mean(high_gsc)
    ratio <- if(expected>0){
      observed/expected
    } else {
      NA
    }
    
    data.frame(sample = sample_name, observed = observed, expected = expected,
               ratio = ratio, pvalue = fisher_result$p.value)
  })
  
  dplyr::bind_rows(results)

}