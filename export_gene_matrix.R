gc()
rm(list = ls())

set.seed(411)
#####----------------------------------------------------------------------#####
##### PACKAGE INSTALLATION
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/UMP_oral_cancer/official/ump_oral_cancer_project"

source(file.path(path.to.main.src, "import_libraries.R"))
source(file.path(path.to.main.src, "helper_functions.R"))

library(DESeq2)

new.packages <- c("umap", "corrr", "factoextra")
for (pkgs in new.packages){
  if (pkgs %in% installed.packages() == FALSE){
    install.packages(pkgs)
  }  else {
    lapply(c(pkgs), 
           require, 
           character.only = TRUE)
  }
} 
library(umap)
if ("ConsensusClusterPlus" %in% installed.packages() == FALSE){
  BiocManager::install("ConsensusClusterPlus")
  BiocManager::install("sva")
}

library("ConsensusClusterPlus")
library("sva")

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
data.version <- "20240630"
code.version <- "v17"
output.version <- "focus_v17_20240703"

#####----------------------------------------------------------------------#####
##### PATHS
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/GSHD_HN01/outdir"
PROJECT <- "UMP_Oral_cancer"

path.to.main.input <- "/media/hieunguyen/GSHD_HN01/raw_data/UMP_Oral_cancer/input"
path.to.main.output <- file.path(outdir, PROJECT, output.version)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output", data.version, code.version)

count.matrix.be <- readRDS(file.path(path.to.01.output, "count_matrix_all_samples_all_genes_vst.rds"))

write.csv(count.matrix.be %>% as.data.frame() %>% rownames_to_column("Gene"), file.path(path.to.01.output, "all_gene_expression_matrix.csv"))