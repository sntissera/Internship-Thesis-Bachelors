
#Load the packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(here)

#Load fucntions
source(here("R","subtype_overlap_function.R"))

#GBM subtype markers
subtypes <- list(MES = c("YKL40","CD44","STAT3","CHI3L1","NF1","PTEN","MET","VIM"), 
                 PN =  c("PDGFRA","OLIG2","NKX2-2","DCX","SOX","ASCL1","DLL3","TCF4"),
                 CL =c("EGFR","NES","NOTCH3","JAG1","LFNG","SMO","GAS1","GLI2"),
                 NL = c("NEFL","SYT1","GABRA1","SLC12A5"))
#MES: mesenchyma, PN: Proneural, CL: Classical, NL:neural 

#Glioma stem like(GSC) markers
gsc <- list(GSC=c("CD133","CD15","CD44","CD70","S100A4","CD49f","L1CAM","A2B5","SOX2","NANOG","OLIG2","MYC","MSI1","NES","OCT4","ALDH1A3"))

#analysis parameters
high_expression_quantile <- 0.75

#Load input
input_file <- here("data","processed","ranking_and_classification.rds")
entropy_vals <- readRDS(input_file)
sample_names <- names(entropy_vals)

#output directories
figure_dir <- here("results","figures","expression_patterns")
table_dir <- here("results","tables","expression_patterns")
dir.create(table_dir,recursive=TRUE,showWarnings = FALSE)
dir.create(figure_dir,recursive=TRUE,showWarnings = FALSE)

#GSC and GBM subtype module scores
marker_availability <- vector("list",length(entropy_vals))
names(marker_availability) <- sample_names
for (i in seq_along(entropy_vals)){
  DefaultAssay(entropy_vals[[i]]) <- "SCT"
  available_genes <- rownames(entropy_vals[[i]][["SCT"]])
  signature_list <- c(gsc,subtypes)
  availability <- lapply(signature_list,function(markers){
    present <- intersect(markers,available_genes)
    data.frame(total_markers = length(markers), available_markers = length(present),
               missing_markers = length(setdiff(markers,available_genes)),
               stringsAsFactors = FALSE)
  })
  availability <- do.call(rbind,availability)
  availability$signature <- names(signature_list)
  availability$sample <- sample_names[i]
  marker_availability[[i]] <- availability
  
  entropy_vals[[i]] <- AddModuleScore(entropy_vals[[i]],features = gsc,name="GSC")
  entropy_vals[[i]] <- AddModuleScore(entropy_vals[[i]],features = subtypes,name="SubType")
 
#Plots 
  image_theme <- theme(plot.title = element_text(size = 14,face = "bold",hjust = 0.5,margin = margin(b = 15)),
                       plot.margin = margin(20, 20, 20, 20), panel.background = element_blank(), plot.background = element_blank(),
                       legend.position = "right",legend.direction = "vertical",legend.title = element_text(size = 10,face = "bold",margin = margin(b = 8)),
                       legend.text = element_text(size = 9),legend.key.height = grid::unit(0.8, "cm"),legend.key.width = grid::unit(0.4, "cm"))
  
  
  
  p1 <- SpatialFeaturePlot(entropy_vals[[i]],features="mean_entropy")+
    labs(title= paste0(names(entropy_vals)[i]," Entropy")) + image_theme
  
  p2 <- SpatialFeaturePlot(entropy_vals[[i]],features="GSC1") + 
    labs(title= "Glioma Stem Cells") + image_theme
  
  p3 <- SpatialFeaturePlot(entropy_vals[[i]],features="SubType1") + 
    labs(title = "Mesenchyma") + image_theme
  
  p4 <- SpatialFeaturePlot(entropy_vals[[i]],features="SubType2") + 
    labs(title = "Proneural") + image_theme
  
  p5 <- SpatialFeaturePlot(entropy_vals[[i]],features="SubType3") + 
    labs(title = "Classical") + image_theme
  
  p6 <- SpatialFeaturePlot(entropy_vals[[i]],features="SubType4") + 
    labs(title = "Neural") + image_theme
  
  pair1 <- p1 + p2 + plot_layout(ncol = 2)
  pair2 <- p3 + p4 + plot_layout(ncol = 2)
  pair3 <- p5 + p6 + plot_layout(ncol = 2)
  
  ggsave(filename = file.path(figure_dir, paste0(names(entropy_vals)[i]), "_entropy_gsc.png"), plot = pair1, width = 10, height = 4.5)
  ggsave(filename = file.path(figure_dir, paste0(names(entropy_vals)[i]), "_subtypes_1_2.png"), plot = pair2, width = 10, height = 4.5)
  ggsave(filename = file.path(figure_dir, paste0(names(entropy_vals)[i]), "_subtypes_3_4.png"), plot = pair3, width = 10, height = 4.5)
}

sample_figure_dir <- file.path(figure_dir,sample_names[i])
dir.create(sample_figure_dir,recursive=TRUE,showWarnings=FALSE)
ggsave (filename=file.path(sample_figure_dir,"entropy_gsc.png"),plot=pair1,width=10,height=4.5,dpi=300)
ggsave(filename=file.path(sample_figure_dir,"subtypes_1_2.png"),plot=pair2, width=10,height=4.5,dpi=300)
ggsave(filename=file.path(sample_figure_dir,"subtypes_3_4.png"),plot=pair3, width=10,height=4.5,dpi=300)

module_scores <- here("data","processed","module_scores.rds")
saveRDS(entropy_vals,module_scores)

marker_availability_df <- dplyr::bind_rows(marker_availability)
write.csv(marker_availability_df,file=file.path(table_dir,"marker_availability.csv"))

#correlation analysis
correlation_results <- calculate_signature_correlations(entropy_vals)
correlation_results
write.csv(correlation_results,file=file.path(table_dir,"entropy_correlations.csv"))

# Fisher exact test to test for entropy GSC overlap
fisher_results <- calculate_fisher_overlap(entropy_vals,quantile_cutoff = high_expression_quantile)
fisher_results
write.csv(fisher_results,file=file.path(table_dir,"fisher_exact_test.csv"))

