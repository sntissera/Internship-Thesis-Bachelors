
#Packages

library(here)
library(Seurat)

## Load the functions
source(here("R","entropy_calculation_function.R"))

##Input

input_file <- here("data","raw","Ren")
file_paths <- list.dirs(input_file,full.names=TRUE,recursive = FALSE)
sample_names <- basename(file_paths)

message("Number of samples: ",length(sample_names))
message("Samples detected: ",paste(sample_names,collapse = ", "))
file_name <- "filtered_feature_bc_matrix.h5"

#Process samples

entropy_vals <- vector("list",length(sample_names))
names(entropy_vals) <- sample_names

for (i in seq_along(file_paths)) {
  
  cat("Sample",i,":",sample_names[i],"\n")
  
  #Load ST data
  gbm <- Load10X_Spatial(data.dir = file_paths[i],filename = file_name)
  
  #Remove zero UMI
  gbm <- subset(gbm, subset=nCount_Spatial > 0)
  #Normalization
  gbm <- SCTransform(gbm,assay ="Spatial", verbose=FALSE)
  
  #Calculate transcriptional entropy
  gbm_ent <- calculate_entropy(gbm)
  
  #Store entropy as a new assay
  gbm[["Entropy"]] <- CreateAssayObject(counts = as(gbm_ent,"CsparseMatrix"))
  
  #Store processed sample
  entropy_vals[[sample_names[i]]] <- gbm
  
}

#Save results

output_file <- here("data","processed","entropy_vals_ren.rds")
saveRDS(entropy_vals,file=output_file)


