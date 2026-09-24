
#Load packages
library(GEOquery)
library(here)

#Create directory to store downloaded dataset
raw_dir <- here("data","raw")
if(!dir.exists(raw_dir)){
  dir.create(raw_dir,recursive=TRUE)
}

#Download supplementary files given a GEO accession number
getGEOSuppFiles("GSE237183",makeDirectory = TRUE,baseDir = raw_dir)
