gc()
rm(list = ls())

set.seed(411)
#####----------------------------------------------------------------------#####
##### PACKAGE INSTALLATION
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/UMP_oral_cancer/official"

source(file.path(path.to.main.src, "import_libraries.R"))
source(file.path(path.to.main.src, "helper_functions.R"))
source(file.path(path.to.main.src, "config.R"))

library(DESeq2)
library(ggpubr)

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
output.version <- "focus_v17_20240630"

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
up.count.thres <- config.params[[code.version]][["up.count.thres"]]
down.count.thres <- config.params[[code.version]][["down.count.thres"]]
input.k <- config.params[[code.version]][["input.k"]]
top_variable_genes <- config.params[[code.version]][["top_variable_genes"]]
batch.effect.sw <- config.params[[code.version]][["batch.effect"]]
var.method <- config.params[[code.version]][["var.method"]]

#####----------------------------------------------------------------------#####
##### PATHS
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/HNSD_mini/data/outdir"
PROJECT <- "UMP_oral_cancer"

path.to.main.input <- "/media/hieunguyen/HNSD_mini/data/UMP_Oral_cancer/input"
path.to.main.output <- file.path(outdir, PROJECT, output.version)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output", data.version, code.version)
path.to.04.output <- file.path(path.to.main.output, "04_output", data.version, code.version)
dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)

traindf <- readxl::read_excel(file.path(path.to.main.input, data.version, "Data for model 30.06.24.xlsx"), sheet = "RNAseq") %>%
  subset(select = -c(Group, `No.`))
testdf <- readxl::read_excel(file.path(path.to.main.input, data.version, "Data for model 30.06.24.xlsx"), sheet = "No RNAseq")
colnames(testdf) <- c("SampleID", colnames(testdf)[2:length(colnames(testdf))])

# convert male and female to 1 and 2
traindf <- traindf %>% rowwise() %>%
  mutate(Gender = ifelse(Gender == "Male", 1, 2)) 
traindf[["Pathological diagnosis"]] <- unlist(lapply(
  traindf[["Pathological diagnosis"]], function(x){
    return(str_replace(x, "grade ", "") %>% as.numeric())
  }))

testdf <- testdf %>% rowwise() %>%
  mutate(Gender = ifelse(Gender == "Male", 1, 2))
testdf[["Pathological diagnosis"]] <- unlist(lapply(
  testdf[["Pathological diagnosis"]], function(x){
    return(str_replace(x, "grade ", "") %>% as.numeric())
  }))

write.csv(traindf, file.path(path.to.04.output, "traindf.csv"))
write.csv(testdf, file.path(path.to.04.output, "testdf.csv"))


