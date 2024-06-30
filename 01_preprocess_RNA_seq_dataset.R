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
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

writexl::write_xlsx(as.data.frame(config.params[[code.version]]), file.path(path.to.01.output, "config.params.logs.xlsx"))

#####----------------------------------------------------------------------#####
##### PREPROCESSING KALLISTO INPUT OBJECT
#####----------------------------------------------------------------------#####
kallisto.obj <- readRDS(file.path(path.to.main.input, data.version, "txi.kallisto.RDS"))

all.samples <- colnames(kallisto.obj$counts)

meta.data <- read.csv(file.path(path.to.main.input, data.version, "Metadata.txt"), sep = "\t")
meta.data <- subset(meta.data, meta.data$SampleID %in% all.samples)

cluster.score <- read.csv(file.path(path.to.main.input, data.version, "cluster_score.csv"), sep = ";")

deseq.obj <- DESeqDataSetFromTximport(txi = kallisto.obj, colData = meta.data, design = ~SampleID)
deseq.obj <- estimateSizeFactors(deseq.obj)

deseq.obj <- deseq.obj[rowSums(counts(deseq.obj)) > down.count.thres, ]
vst_mat <- vst(deseq.obj, blind = TRUE)
count.matrix <- assay(vst_mat)

saveRDS(count.matrix, file.path(path.to.01.output, "count_matrix.raw.rds"))
##### batch effect correction
batchdf <- data.frame(sampleid = colnames(count.matrix)) %>%
  rowwise() %>%
  mutate(batch = subset(meta.data, meta.data$SampleID == sampleid)$seq_date)
if (batch.effect.sw == TRUE){
  print("perform batch effect correction with limma")
  count.matrix.be <- limma::removeBatchEffect(x = count.matrix, batch = meta.data$seq_date)    
} else if (batch.effect.sw == FALSE){
  count.matrix.be <- count.matrix
}

saveRDS(count.matrix.be, file.path(path.to.01.output, "count_matrix_all_samples_all_genes_vst.rds"))

##### Select top variable genes only
if (var.method == "mad"){
  gene.mad <- apply(count.matrix.be, 1, mad) %>% sort()
  topVar.genes <- tail(gene.mad, top_variable_genes) %>% names()
} else if (var.method == "none"){
  topVar.genes <- rowSds(count.matrix.be) %>% sort() %>% tail(top_variable_genes) %>% names()    
}

saveRDS(topVar.genes, file.path(path.to.01.output, "selected_topVar_genes.rds"))
count.matrix.be.raw <- count.matrix.be
count.matrix.be <- count.matrix.be[topVar.genes, ]

##### run consensus clustering
title <- file.path(path.to.01.output, "test")
results <- ConsensusClusterPlus(count.matrix.be, 
                                maxK= input.k, 
                                reps=100, 
                                pItem=1,
                                pFeature=1,
                                title=title, 
                                clusterAlg="km",
                                distance="euclidean",
                                seed=411,plot="png")

icl <- calcICL(results)
icldf <- icl$itemConsensus %>% as.data.frame()

return_consensus_cluster <- function(input.k, sample.id, icldf){
  
  tmpdf <- subset(icldf, icldf$k == input.k & icldf$item == sample.id)  
  return(subset(tmpdf, tmpdf$itemConsensus == max(tmpdf$itemConsensus))$cluster)
}

meta.data <- meta.data %>% rowwise() %>% 
  mutate(consensus.cluster = return_consensus_cluster(input.k, SampleID, icldf))

output.metadata <- subset(meta.data, select = c(SampleID, consensus.cluster))
colnames(output.metadata) <- c("SampleID", "RNA.consensus.cluster")

write.csv(output.metadata, file.path(path.to.01.output, "metadata_with_consensus_RNA_clusters.csv"))

##### run UMAP to validate the consensus clustering
umap.res <- umap(count.matrix.be  %>% t())
umapdf <- umap.res$layout %>% as.data.frame() %>%
  rownames_to_column("SampleID")
umapdf <- merge(umapdf, subset(meta.data, select = c(SampleID, consensus.cluster)), by.x = "SampleID", by.y = "SampleID")
colnames(umapdf) <- c("SampleID", "V1", "V2", "RNA.consensus.cluster")
umapdf$RNA.consensus.cluster <- factor(umapdf$RNA.consensus.cluster, levels = unique(umapdf$RNA.consensus.cluster))
umap.plot.RNA <- umapdf %>% ggplot(aes(x = V1, y = V2, color = RNA.consensus.cluster)) + geom_point(size = 2)

ggsave(plot = umap.plot.RNA, filename = sprintf("UMAP_RNAseq_data.%s.svg", code.version), path = path.to.01.output, device = "svg", width = 14, height = 10)

##### Apply k-mean clustering to have 2 clusters only

km.res <- kmeans(count.matrix.be %>% t(), 2, nstart = 1)
umapdf <- umapdf %>% rowwise() %>%
  mutate(kmean.cluster = km.res$cluster[[SampleID]])
umapdf$kmean.cluster <- factor(umapdf$kmean.cluster, levels = c(1, 2))
umap.plot.kmean <- umapdf %>% ggplot(aes(x = V1, y = V2, color = kmean.cluster)) + geom_point(size = 2)

ggsave(plot = umap.plot.kmean, filename = sprintf("UMAP_RNAseq_data_KMean_2_clusters.%s.svg", code.version), path = path.to.01.output, device = "svg", width = 14, height = 10)

cluster.scoredf <- read.csv(file.path(path.to.main.input, data.version, "cluster_score.csv"), sep = ";")[, c("SampleID", "cluster2kmean", "cluster3kmean")]
colnames(cluster.scoredf) <- c("SampleID", "kmean.2clusters.DrNam", "kmean.3clusters.DrNam")
umapdf <- merge(umapdf, cluster.scoredf, by.x = "SampleID", by.y = "SampleID")



umapdf <- umapdf %>% rowwise() %>%
  mutate(merged.cluster12 = ifelse(RNA.consensus.cluster %in% c(1,2), 1, 2)) %>%
  mutate(merged.cluster13 = ifelse(RNA.consensus.cluster %in% c(1,3), 1, 2)) %>%
  mutate(merged.cluster23 = ifelse(RNA.consensus.cluster %in% c(2,3), 1, 2))

sample.order <- umapdf[order(umapdf$merged.cluster13),]$SampleID

write.csv(umapdf, file.path(path.to.01.output, "umap_RNAseq.csv"))

write.table(count.matrix %>% as.data.frame() %>% rownames_to_column("Gene") , file.path(path.to.01.output, "final_input_to_CIBERSORTx.tsv"), sep = "\t", row.names = FALSE)

write.table(count.matrix[topVar.genes, ] %>% as.data.frame() %>% rownames_to_column("Gene") , file.path(path.to.01.output, "final_input_to_CIBERSORTx.topVarGenes.tsv"), sep = "\t", row.names = FALSE)

#####----------------------------------------------------------------------#####
##### UPDATE 15.05.2024
#####----------------------------------------------------------------------#####

##### Plot cluster - marker genes
marker.genes <- c("RESF1",
                  "SLAIN2",
                  "HNRNPD",
                  "ADNP",
                  "SLK",
                  "WAC",
                  "MGAT1",
                  "BAG1",
                  "MARCKS",
                  "OGFR",
                  "GIGYF1",
                  "FKBP8", 
                  "MKI67", 
                  "LAMC1", 
                  "LAMC2", 
                  "MYCL", 
                  "VEGFA", 
                  "CDH3", 
                  "CDH1", 
                  "SCD1", 
                  "CTNNB1")
for (input.gene in intersect(marker.genes, rownames(count.matrix.be.raw))){
  dir.create(file.path(path.to.01.output, "marker_genes_boxplot"), showWarnings = FALSE, recursive = TRUE)
  tmpdf <- count.matrix.be.raw[input.gene, ] %>% as.data.frame()
  colnames(tmpdf) <- "exprs"
  tmpdf <- tmpdf %>% rownames_to_column("Sample") %>%
    rowwise() %>%
    mutate(cluster.k3 = subset(umapdf, umapdf$SampleID == Sample)$RNA.consensus.cluster) %>%
    mutate(cluster.merge = subset(umapdf, umapdf$SampleID == Sample)$merged.cluster13)
  tmpdf$cluster.k3 <- factor(tmpdf$cluster.k3, levels = c(1,2,3))
  tmpdf$cluster.merge <- factor(tmpdf$cluster.merge, levels = c(1,2))
  tmpdf %>% ggplot(aes(x = cluster.k3, y = exprs, fill = cluster.k3)) + geom_boxplot()
  p <- tmpdf %>% ggplot(aes(x = cluster.k3, y = exprs, fill = cluster.k3)) + geom_boxplot() + geom_jitter(width = 0.02) + theme_pubr()
  ggsave(plot = p, path = file.path(path.to.01.output, "marker_genes_boxplot"), filename = sprintf("Gene_%s.3clusters.svg", input.gene), device = "svg", width = 10, height = 10, dpi = 300)
  p <- tmpdf %>% ggplot(aes(x = cluster.merge, y = exprs, fill = cluster.merge)) + geom_boxplot() + geom_jitter(width = 0.02) + theme_pubr()
  ggsave(plot = p, path = file.path(path.to.01.output, "marker_genes_boxplot"), filename = sprintf("Gene_%s.merged_2_clusters.svg", input.gene), device = "svg", width = 10, height = 10, dpi = 300)
}

##### CIBERSORTx results
path.to.cibersort.output <- "/media/hieunguyen/HNSD_mini/data/outdir/UMP_Oral_cancer/focus_v17_20240425/CIBERSORTx"
all.files <- Sys.glob(file.path(path.to.cibersort.output, "*"))

for (selected.cluster.label in c("RNA.consensus.cluster", "merged.cluster13")){
  dir.create(file.path(path.to.01.output, selected.cluster.label), showWarnings = FALSE, recursive = TRUE)
  
  for (file in all.files){
    output.name <- str_replace(basename(file), ".csv", ".svg")
    tmpdf <- read.csv(file) %>% subset(select = -c(P.value, Correlation, RMSE))
    tmpdf <- tmpdf %>% rowwise() %>%
      mutate(cluster = subset(umapdf, umapdf$SampleID == Mixture)[[selected.cluster.label]])
    tmpdf$cluster <- factor(tmpdf$cluster, levels = seq(1, length(unique(umapdf[[selected.cluster.label]]))))
    tmpdf <- tmpdf %>% column_to_rownames("Mixture") 
    p <- tmpdf[, c("Malignant", "cluster")] %>% ggplot(aes(x = cluster, y = Malignant, fill = cluster)) + geom_boxplot() + theme_pubr()
    ggsave(plot = p, path = file.path(path.to.01.output, selected.cluster.label), filename = output.name, device = "svg", width = 10, height = 10, dpi = 300)
  }
}