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
library(umap)
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

if ("ConsensusClusterPlus" %in% installed.packages() == FALSE){
  BiocManager::install("ConsensusClusterPlus")
  BiocManager::install("sva")
}

library("ConsensusClusterPlus")

# BiocManager::install("msigdbr", update = FALSE)
# BiocManager::install("org.Hs.eg.db", update = FALSE)
# remove.packages("clusterProfiler")
# remove.packages("DOSE")
# remove.packages("GOSemSim")
# path.to.install.dir <- "/media/hieunguyen/HNHD01/offline_pkgs/clusterProfiler"
# BiocManager::install("fgsea", update = FALSE)
# install.packages("tidygraph")
# BiocManager::install("enrichplot", update = FALSE)
# install.packages("ggplotify")
# for (pkg in c("HDO.db_0.99.1.tar.gz",
#               "yulab.utils_0.1.4.tar.gz",
#               "GO.db_3.18.0.tar.gz",
#               "GOSemSim_2.28.1.tar.gz",
#               "DOSE_3.28.2.tar.gz",
#               "gson_0.1.0.tar.gz",
#               "downloader_0.4.tar.gz",
#               "clusterProfiler_4.10.1.tar.gz")){
#     install.packages(file.path(path.to.install.dir, pkg), type = "source", repos = NULL)
# }

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
data.version <- "240319"
code.version <- "v17"
output.version <- "focus_v17_20240425"

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
up.count.thres <- config.params[[code.version]][["up.count.thres"]]
down.count.thres <- config.params[[code.version]][["down.count.thres"]]
input.k <- config.params[[code.version]][["input.k"]]
top_variable_genes <- config.params[[code.version]][["top_variable_genes"]]

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
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)


cluster.scoredf <- read.csv(file.path(path.to.main.input, data.version, "cluster_score_old.csv"), sep = ",")
convert.colname.clusterscore <- list(
  "Độ.sâu.xâm.lấn..DOI." = "feature1",
  "Độ.sừng.hóa.tế.bào" = "feature2",
  "Dị.dạng.nhân" = "feature3",
  "Kiểu.xâm.lấn" = "feature4",
  "Thấm.nhập.lympho.tương.bào" = "feature5",
  "Độ.ác.tính.mô.học.theo.Bryne" = "feature6",
  "Phân.độ.mô.học.theo.WHO" = "feature7",
  "Xâm.nhiễm.thần.kinh..PNI." = "feature8",
  "Kiểu.hình.xâm.nhiễm.xấu.nhất..WPOI."  = "feature9" ,
  "Xâm.nhiễm.lympho..LHR." = "feature10",
  "Độ.ác.tính.mô.học.theo.Brandwwein.Gensler" = "feature11",
  "Chẩn.đoán.vị.trí" = "feature12",
  "Kích.thước.bướu.lâm.sàng..cT." = "feature13",
  "Di.căn.hạch.lâm.sàng..cN." = "feature14",
  "Di.căn.xa" = "feature15",
  "Giai.đoạn.lâm.sàng..TNM." = "feature16",
  "Hút.thuốc" = "feature17",
  "Uống.rượu" = "feature18",
  "Chẩn.đoán.giải.phẫu.bệnh" = "feature19",
  "Di.căn.hạch.giải.phẫu.bệnh..pN." = "feature20",
  "Kiểu.tăng.trưởng...Xem.đại.thể." = "feature21"
)

remove.feature.names <- c(
  "Giai.đoạn.lâm.sàng..TNM.",
  "Hút.thuốc",
  "Uống.rượu"
)

remove.features <- convert.colname.clusterscore[remove.feature.names] %>% unlist()

colnames(cluster.scoredf) <- unlist(lapply(colnames(cluster.scoredf), function(x){
  if (x %in% names(convert.colname.clusterscore)){
    return(convert.colname.clusterscore[[x]])
  } else {
    return(x)
  }
}))
cluster.scoredf <- subset(cluster.scoredf, select = -c(X))
featuredf <- cluster.scoredf[, c(3, seq(34, 39), 41, 43, seq(44, 52), seq(54, 57))]

##### manually convert feature from "character" to numerical values
featuredf$feature13 <- to_vec( for(item in featuredf$feature13) as.numeric(str_replace(item, "T", "")))
featuredf$feature14 <- to_vec( for(item in featuredf$feature14) as.numeric(str_replace(item, "N", "")))
featuredf$feature19 <- to_vec( for(item in featuredf$feature19) as.numeric(str_replace(item, "grad ", "")))
convert.feature16 <- seq(1, length(unique(featuredf$feature16)))
names(convert.feature16) <- unique(featuredf$feature16)
featuredf <- featuredf %>% rowwise() %>%
  mutate(feature15 = ifelse(feature15 == "M0", 1, 2)) %>%
  mutate(feature16 = convert.feature16[[feature16]])
featuredf[is.na(featuredf)] <- 0

meta.data <- read.csv(file.path(path.to.01.output, "metadata_with_consensus_RNA_clusters.csv")) %>%
  subset(select = -c(X))
colnames(meta.data) <- c("SampleID", "RNA.consensus.cluster")

featuredf.raw <- featuredf
write.csv(featuredf, file.path(path.to.05.output, "featuredf.raw.csv"))

write.csv(featuredf[, setdiff(colnames(featuredf), c("feature16"))], file.path(path.to.05.output, "featuredf.remove_TNM.csv"))

gpb.features <- unlist(lapply(seq(1, 13), function(x){return(sprintf("feature%s", x))}))
featuredf13 <- featuredf[, c(gpb.features, c("SampleID"))]
write.csv(featuredf13, file.path(path.to.05.output, "featuredf.GPB_only.csv"))

featuredf <- featuredf[, setdiff(colnames(featuredf), remove.features)]
write.csv(featuredf, file.path(path.to.05.output, "featuredf.filtered.csv"))

#####----------------------------------------------------------------------#####
##### UMAP on score input features
#####----------------------------------------------------------------------#####

umap.input.mat <- featuredf %>% column_to_rownames("SampleID") 
umap.res <- umap(umap.input.mat)
umapdf <- umap.res$layout %>% as.data.frame() %>% 
  rownames_to_column("SampleID")

pca.output <- prcomp(t(umap.input.mat), scale. = FALSE)
pcadf <- pca.output$rotation[, c("PC1", "PC2")] %>% as.data.frame()

##### assign clusters based on UMAP of featuredf
title <- file.path(path.to.05.output, "test")
results <- ConsensusClusterPlus(featuredf %>% column_to_rownames("SampleID") %>% as.matrix() %>% t(), 
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

featuredf <- featuredf %>% rowwise() %>%
  mutate(score.consensus.cluster = return_consensus_cluster(input.k, SampleID, icldf))

umapdf$score.consensus.cluster <- featuredf$score.consensus.cluster
umapdf$score.consensus.cluster <- factor(umapdf$score.consensus.cluster, levels = seq(1, input.k))
umapdf$score.consensus.cluster <- factor(umapdf$score.consensus.cluster, levels = seq(1, input.k))
umap.plot.score <- umapdf %>% ggplot(aes(x = V1, y = V2, color = score.consensus.cluster)) + geom_point(size = 2)
umap.plot.score.RNA.overlay <- umapdf %>% ggplot(aes(x = V1, y = V2, color = RNA.consensus.cluster)) + geom_point(size = 2)

write.csv(umapdf, file.path(path.to.05.output, "umap_from_score.csv"))

