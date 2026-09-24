
#Shannon Entropy Calculation

calculate_entropy <- function(gbm){
  gs_matrix <- GetAssayData(gbm,assay="SCT", layer="counts") #raw counts gene-spot matrix
  total_count <- colSums(gs_matrix) #total number of counts per spot
  total_count[total_count==0] <- NA #avoid division by zero
  p <- sweep(gs_matrix,2,total_count,FUN="/") #pdf across spots per gene
  gs_ent <- -(p*log2(p)) #calculates the entropy
  gs_ent[!is.finite(gs_ent)] <- 0 #replaces NaN or Inf given by 0log0 to 0
  return (gs_ent)
}

#Entropy Ranking

#Calculate global entropy values
  gbm_calc <- function(ent_val){
  global_min_ent <- Inf
  global_max_ent <- -Inf
  
  for (i in seq_along(ent_val)){
    gbm_ent <- ent_val[[i]] 
    gs.entropy <- GetAssayData(gbm_ent, assay="Entropy",layer="counts")  #Gene x Spot entropy matrix
    global_min_ent <- min(global_min_ent, min(gs.entropy, na.rm = TRUE))
    global_max_ent <- max(global_max_ent, max(gs.entropy, na.rm = TRUE))
  }
  gbm_global <- c(min_entropy = global_min_ent, max_entropy = global_max_ent)
  return(gbm_global)
}

#Scale entropy values
gbm_scale <- function(ent_val,gbm_global){ 
  gbm.data <- vector("list",length(ent_val))
  for (i in seq_along(ent_val)){
    gbm <- GetAssayData(ent_val[[i]], assay = "Entropy", layer = "counts")
    ent_scaled <- (gbm - gbm_global[1]) / (gbm_global[2] - gbm_global[1]) 
    gbm.data[[i]] <- ent_scaled
    ent_val[[i]][["Scaled_Entropy"]] <- CreateAssayObject(counts = as(ent_scaled,"CsparseMatrix"))
    DefaultAssay(ent_val[[i]]) <- "Scaled_Entropy"
  }
  names(gbm.data) <- names(ent_val)
  return(list(entropy_data = gbm.data, entropy_vals=ent_val))
}

#Calculate entropy statisitcs
entropy_stats <- function(scale_data, entropy_vals){ 
  entropy_df <- vector("list",length(scale_data))
  
  for (i in seq_along(scale_data)){
    mean_entropy <- colMeans(scale_data[[i]], na.rm = TRUE)
    sd_entropy <- matrixStats::colSds(as.matrix(scale_data[[i]]),na.rm=TRUE)
    entropy_df[[i]] <- data.frame(sample = names(scale_data)[i], spot = colnames(scale_data[[i]]), 
                                  mean_entropy=mean_entropy,sd_entropy = sd_entropy)
  }
  
  entropy_data <- dplyr::bind_rows(entropy_df)
  return (entropy_data)
}

#Add mean and sd entropy to metadata
add_to_metadata <- function(entropy_data,entropy_vals){
  for (i in seq_along(entropy_vals)){
    sample_data <- entropy_data[entropy_data$sample==names(entropy_vals)[i],]
    mean_ent <- sample_data$mean_entropy
    names(mean_ent) <- sample_data$spot
    std <- sample_data$sd_entropy
    names(std) <- sample_data$spot
    
    entropy_vals[[i]] <- AddMetaData(entropy_vals[[i]],metadata = mean_ent, col.name = "mean_entropy")
    entropy_vals[[i]] <- AddMetaData(entropy_vals[[i]],metadata = std, col.name = "sd_entropy")
  }
  return(entropy_vals)
}

#detect outliers
detect_outliers <- function(entropy_data){
  samples <- unique(entropy_data$sample)
  outlier_data <- lapply(samples, function(sample_name){
    sample_data <- entropy_data[entropy_data$sample == sample_name,]
    q1 <- quantile(sample_data$score,0.25)
    q3 <- quantile(sample_data$score,0.75)
    iqr_val <- q3 - q1
    low_bound <- q1 - 1.5 * iqr_val
    up_bound <- q3 + 1.5 * iqr_val
    sample_data[sample_data$score < low_bound| sample_data$score > up_bound,]
  })
  do.call(rbind,outlier_data)
}

#classification using GMM
gmm_classification <- function(entropy_data){
  split_data <- split(entropy_data,entropy_data$sample)
  gmm_models <- list()
  classified_data <- list()
  
  for (sample_name in names(split_data)){
    sample_data <- split_data[[sample_name]]
    fit <- Mclust(cbind(sample_data$mean_entropy,sample_data$sd_entropy))
    sample_data$group <- fit$classification 
    classified_data[[sample_name]] <- sample_data
    gmm_models[[sample_name]] <- fit
  } 
  return (list(data=do.call(rbind, classified_data), models=gmm_models))
}