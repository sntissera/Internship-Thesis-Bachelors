
# Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(matrixStats)
library(mclust)
library(here)

# Load the functions
source(here("R","entropy_calculation_function.R"))

# Load the samples
input_file <- here("data","processed","entropy_vals_ren.rds")
entropy_vals <- readRDS(input_file)
sample_names <- names(entropy_vals)

#Generate global statistics
gbm_global <- gbm_calc(entropy_vals)
scale_vals <- gbm_scale(entropy_vals,gbm_global)
entropy_data <- entropy_stats(scale_vals$entropy_data,scale_vals$entropy_vals)
entropy_vals <- add_to_metadata(entropy_data,scale_vals$entropy_vals)

#Calculate entropy score 
for (i in seq_along(entropy_vals)){
  score <- entropy_vals[[i]]$mean_entropy + entropy_vals[[i]]$sd_entropy
  entropy_vals[[i]] <- AddMetaData(entropy_vals[[i]],metadata = score,
                                          col.name = "Score")
  sample_index <- entropy_data$sample == sample_names[i]
  entropy_data$score[sample_index] <- score
}

#Boxplot 
entropy_boxplot <- ggplot(data=entropy_data, mapping = aes(x=sample,y=mean_entropy,fill=sample)) + geom_boxplot() +  
  stat_summary(fun = "mean",geom = "point", shape = 1, size = 2, color = "black") +
  labs(title = "Mean Entropy distribution across glioblastoma samples", x = "Sample", y = "Mean Entropy") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
entropy_boxplot

ggsave(filename = here("results", "plots", "Ren", "mean_entropy_boxplot.png"),plot = entropy_boxplot,
       width = 8,height = 6,dpi = 300) #save the boxplot

#Boxplot statistics

#median entropy per sample
sample_median <- sapply(sample_names,function(sample_name){
  median(entropy_data$mean_entropy[entropy_data$sample==sample_name])
})

sample_median <- sort(sample_median, decreasing=TRUE) #Provide the sample medians in descending order
sample_median

highest_entropy <- names(which.max(sample_median)) #sample with highest median entropy
highest_entropy

lowest_entropy <- names(which.min(sample_median)) #sample_with lowest median entropy
lowest_entropy

#entropy variability between samples
iqr_summary <- data.frame(sample=sample_names,IQR=sapply(sample_names,function(sample_name){
  IQR(entropy_data$mean_entropy[entropy_data$sample==sample_name])
}))

iqr_summary <- iqr_summary[order(iqr_summary$IQR,decreasing = TRUE),] #IQR in descending order
rownames(iqr_summary) <- NULL
iqr_summary

#Box plot (mean + sd)
score_boxplot <- ggplot(data=entropy_data, mapping = aes(x=sample,y=score,fill=sample)) + geom_boxplot() +        
  stat_summary(fun = "mean",geom = "point", shape = 1, size = 2, color = "black") +
  labs(title = "Mean + Standard deviation of Entropy distribution across glioblastoma samples", x = "Sample", y = "Mean + SD Entropy") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
score_boxplot

ggsave(filename = here("results", "plots", "Ren", "entropy_score_boxplot.png"),plot = score_boxplot,
       width = 8,height = 6,dpi = 300)

score_median <- sapply(names(entropy_vals),function(sample_name){
  median(entropy_data$score[entropy_data$sample==sample_name],na.rm=TRUE)
})
sort(score_median,decreasing =TRUE) #Median entropy in decreasing order
names(which.max(score_median)) #highest median entropy sample
names(which.min(score_median)) #lowest median entropy sample

#Mean entropy Violin plot
entropy_violin <- ggplot(data=entropy_data, mapping = aes(x=sample,y=mean_entropy,fill=sample)) + geom_violin() +         
  stat_summary(fun = "mean",geom = "point", shape = 1, size = 2, color = "black") +
  labs(title = "Mean Entropy distribution across glioblastoma samples", x = "Sample", y = "Mean Entropy") + 
  theme(axis.text.x = element_text(angle=45,hjust=1))
entropy_violin

ggsave(filename = here("results", "plots", "Ren", "mean_entropy_violinplot.png"),plot = entropy_violin,
       width = 8,height = 6,dpi = 300)

#Score Violin plot
score_violin <- ggplot(data=entropy_data, mapping = aes(x=sample,score,fill=sample)) + geom_violin() +
  stat_summary(fun = "mean",geom = "point", shape = 1, size = 2, color = "black") +
  labs(title = "Mean + Standard deviation distribution across glioblastoma samples", x = "Sample", y = "Mean + SD Entropy") + 
  theme(axis.text.x = element_text(angle=45,hjust=1))
score_violin

ggsave(filename = here("results", "plots","Ren", "entropy_score_violinplot.png"),plot = score_violin,
       width = 8,height = 6,dpi = 300)

#detect outliers
outlier_data <- detect_outliers(entropy_data)
entropy_data$new_spot <- paste0(entropy_data$spot,"_",entropy_data$sample) #new spot lables to distinguish between similar barcodes across samples
outlier_data$new_spot <- paste0(outlier_data$spot,"_",outlier_data$sample)
entropy_data$is_outlier <- entropy_data$new_spot %in% outlier_data$new_spot #Add to entropy data if outlier or not

#Add outliers into metadata and count them
for(i in seq_along(entropy_vals)){
  outlier_status <- entropy_data[entropy_data$sample == sample_names[i],]$is_outlier
  names(outlier_status) <- entropy_data[entropy_data$sample == sample_names[i],]$spot
  entropy_vals[[i]] <- AddMetaData(entropy_vals[[i]],metadata=outlier_status,col.name = "is_outlier")
  outlier_number <- sum(entropy_vals[[i]]$is_outlier,na.rm=TRUE)
  out <-paste0(sample_names[i],":",outlier_number)
  print(out)
}

#gmm classification
entropy_clean <- entropy_data[entropy_data$is_outlier==FALSE,]
set.seed(123)
gmm <- gmm_classification(entropy_clean)
entropy_clean <- gmm$data
gmm_models <- gmm$models

#classification statistics
gmm_stat <- lapply(names(gmm_models),function(sample_name){
  fit <- gmm_models[[sample_name]]
  means <- fit$parameters$mean
  data.frame(sample=sample_name,component = seq_len(fit$G),mean_entropy = means[1, ],
             sd_entropy = means[2, ],component_size = as.numeric(table(fit$classification)))
}) |> 
  dplyr::bind_rows()
gmm_stat

#Mean entropy across outliers per sample
m.persample <- data.frame(sample=names(entropy_vals))
for (i in seq_along(entropy_vals)){
  entropy_dt <- entropy_data[entropy_data$sample == names(entropy_vals)[i],]
  outlier_spots <- entropy_dt[entropy_dt$is_outlier == TRUE,]
  m.persample$m_oe[i] <- mean(outlier_spots$mean_entropy) #m_oe -> mean outlier entropy per sample
  m.persample$sd_oe[i] <- mean(outlier_spots$sd_entropy)
  m.persample$score[i] <- mean(outlier_spots$score)
}
m.persample

#Add the regions into metadata
for(i in seq_along(entropy_vals)){
  category <- entropy_clean[entropy_clean$sample == sample_names[i],]$group
  names(category) <- entropy_clean[entropy_clean$sample == sample_names[i],]$spot
  entropy_vals[[i]] <- AddMetaData(entropy_vals[[i]],metadata = category, col.name = "group")
  entropy_vals[[i]]$group[entropy_vals[[i]]$is_outlier == TRUE] <- "outlier" #add as outlier to the spots without a region
}

#save classified groups
output_file <- here("data","processed","ranking_and_classification_ren.rds")
saveRDS(entropy_vals,output_file)

#Visualise discrete entropy regions
for (i in seq_along(entropy_vals)){
  groups <- sort(unique(entropy_vals[[i]]$group))
  n_groups <- length(groups)
  colours <- setNames(colorRampPalette(c("royalblue", "lightblue", "lightyellow", "orange", "red", "darkred"))(n_groups),groups)
  spatial_plot <- print((SpatialDimPlot(entropy_vals[[i]], group.by = "group", 
                        cols = colours)) +
          labs(title = sample_names[i]))
  ggsave(filename = here("results","figures","classification","Ren",paste0(sample_names[i],".png")),
         plot = spatial_plot,width = 8,height = 6,dpi = 300)
}

