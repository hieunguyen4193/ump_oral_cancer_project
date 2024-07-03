# gc()
# rm(list = ls())

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
output.version <- "focus_v17_20240703"

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
path.to.05.output <- file.path(path.to.main.output, "05_output", data.version, code.version)
path.to.06.output <- file.path(path.to.main.output, "06_output", data.version, code.version)
dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)

library(heatmaply)

##### Training model results
cv.scoredf <- readxl::read_excel(file.path(path.to.05.output, "all_CV_scores_final.xlsx")) %>% subset(select = -c(`...1`))
cv.scoredf <- colMeans(cv.scoredf) %>% as.data.frame() %>%
  rownames_to_column("Model")
colnames(cv.scoredf) <- c("Model", "Avg_Accuracy_10FoldCV")
p <- cv.scoredf %>% ggplot(aes(x = Model, y = Avg_Accuracy_10FoldCV, fill = Model)) + geom_bar(stat = "identity") + theme_pubr()
ggsave(plot = p, filename = "compare_model.svg", path = file.path(path.to.01.output), device = "svg", width = 14, height = 10, dpi = 300)

##### plot heat map
path.to.deseq.res <- file.path(path.to.main.output, "02_output", data.version, code.version, "merged.cluster13", "DESEQ2_RNAseq_data_one_vs_one.obj.rds")
dds <- readRDS(path.to.deseq.res)
resdf <- results(dds$cluster1_vs_cluster2) %>% 
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  subset(padj <= 0.05) %>%
  rowwise() %>%
  mutate(abs.log2FoldChange = abs(log2FoldChange)) %>%
  subset(abs.log2FoldChange >= 1)

top100.up <- subset(resdf, resdf$log2FoldChange > 0) %>% arrange(desc(abs.log2FoldChange)) %>% head(25)
top100.down <- subset(resdf, resdf$log2FoldChange < 0) %>% arrange(desc(abs.log2FoldChange)) %>% head(25)
top.genedf <- rbind(top100.up, top100.down)

umapdf <- read.csv(file.path(path.to.01.output, "umap_RNAseq.csv"))
sample.order <- umapdf[order(umapdf$merged.cluster13),]$SampleID

full.count.matrix <- counts(dds$cluster1_vs_cluster2, normalized = TRUE)
dds.count <- counts(dds$cluster1_vs_cluster2, normalized = TRUE)[top.genedf$Gene, sample.order]
count.matrix.cluster1 <- full.count.matrix[, subset(colData(dds$cluster1_vs_cluster2), colData(dds$cluster1_vs_cluster2)$merged.cluster13 == 1) %>% row.names()]
count.matrix.cluster2 <- full.count.matrix[, subset(colData(dds$cluster1_vs_cluster2), colData(dds$cluster1_vs_cluster2)$merged.cluster13 == 2) %>% row.names()]
write.csv(full.count.matrix, file.path(path.to.01.output, "full_count_matrix_merged_cluster13.csv"))
write.csv(count.matrix.cluster1, file.path(path.to.01.output, "count_matrix_merged_clsuter13_cluster1.csv"))
write.csv(count.matrix.cluster2, file.path(path.to.01.output, "count_matrix_merged_clsuter13_cluster2.csv"))

dds.count.scaled <- (dds.count - rowMeans(dds.count))/rowSds(dds.count)
dds.count.scaled.pivot <- dds.count.scaled %>% 
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(!Gene, names_to = "Sample", values_to = "exprs") %>%
  rowwise() %>%
  mutate(cluster = subset(umapdf, umapdf$SampleID == Sample)$merged.cluster13)

p <- ggheatmap(dds.count, scale = "row", dend = "none",
               row_dend_left = FALSE, row_text_angle = 0,
               column_text_angle = 90,
               scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                 low = "blue",
                 high = "red"))

ggsave(plot = p, filename = "heatmap_DEG.svg", path = file.path(path.to.01.output), device = "svg", width = 20, height = 15, dpi = 300)