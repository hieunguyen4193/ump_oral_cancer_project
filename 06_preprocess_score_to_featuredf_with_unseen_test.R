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
#   install.packages(file.path(path.to.install.dir, pkg), type = "source", repos = NULL)
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
path.to.06.output <- file.path(path.to.main.output, "06_output", data.version, code.version)
dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)

# cluster.scoredf <- readxl::read_excel(file.path(path.to.main.input, data.version, "Supplementary document S1.xlsx"), sheet = "Samples with RNAseq")

cluster.scoredf <- readxl::read_excel(file.path(path.to.main.input, data.version, "Supplementary document S1 28.05.24.xlsx"), sheet = "Samples with RNAseq")
cluster.scoredf.raw <- cluster.scoredf
cluster.scoredf$Gender <- to_vec( for(item in cluster.scoredf$Gender) if (item == "Male") 1 else 0)
cluster.scoredf$`Smoking tobacco` <- to_vec( for(item in cluster.scoredf$`Smoking tobacco`) if (item == "NA") 0 else item)
cluster.scoredf$`Alcohol drinking` <- to_vec( for(item in cluster.scoredf$`Alcohol drinking`) if (item == "NA") 0 else item)
cluster.scoredf.raw <- cluster.scoredf
cluster.test.scoredf <- readxl::read_excel(file.path(path.to.main.input, data.version, "Supplementary document S1.xlsx"), sheet = "Samples with no RNAseq")
cluster.test.scoredf$Gender <- to_vec( for(item in cluster.test.scoredf$Gender) if (item == "Male") 1 else 0)
cluster.test.scoredf$`Smoking tobacco` <- to_vec( for(item in cluster.test.scoredf$`Smoking tobacco`) if (item == "NA") 0 else item)
cluster.test.scoredf$`Alcohol drinking` <- to_vec( for(item in cluster.test.scoredf$`Alcohol drinking`) if (item == "NA") 0 else item)

all.cols <- c("Gender", "Primary tumor site",
                  "Clinical tumor size (cT)", "Clinical node mestatasis (cN)", "Distant metastasis (M)",
                  "TNM classification", "Smoking tobacco", "Alcohol drinking",
                  "Pathological diagnosis", "Pathological lymph node metastasis (pN)", "Growth patterns",
                  "Depth of invasion (DOI)", "Degree of keratinization", "Nuclear polymorphism",
                  "Pattern of invasion", "Lymphoplasmacytic \r\ninfiltration", "Bryne scores",
                  "WHO system", "Perineural invasion (PNI)", "Worst pattern of invasion (WPOI)",
                  "Lymphocytic host response (LHR)", "BrandweinGensler risk level")
all.features <- to_vec( for (item in seq(1, length(all.cols))) sprintf("feature%s", item))
names(all.features) <- all.cols

colnames(cluster.scoredf) <- unlist(lapply(colnames(cluster.scoredf), function(x){
  if (x %in% names(all.features)){
    return(all.features[[x]])
  } else {
    return(x)
  }
}))

colnames(cluster.test.scoredf) <- unlist(lapply(colnames(cluster.test.scoredf), function(x){
  if (x %in% names(all.features)){
    return(all.features[[x]])
  } else {
    return(x)
  }
}))

cluster.test.scoredf <- cluster.test.scoredf[, c("No.", all.features)]
featuredf <- cluster.scoredf[, c("BlockID", all.features)]

convert.feature.to.num <- function(featuredf){
  featuredf$feature3 <- to_vec( for(item in featuredf$feature3) as.numeric(str_replace(item, "T", "")))
  featuredf$feature4 <- to_vec( for(item in featuredf$feature4) as.numeric(str_replace(item, "N", "")))
  featuredf$feature9 <- to_vec( for(item in featuredf$feature9) as.numeric(str_replace(item, "grad ", "")))
  convert.feature6 <- seq(1, length(unique(featuredf$feature6)))
  names(convert.feature6) <- unique(featuredf$feature6)
  featuredf <- featuredf %>% rowwise() %>%
    mutate(feature5 = ifelse(feature5 == "M0", 1, 2)) %>%
    mutate(feature6 = convert.feature6[[feature6]])
  featuredf$feature1 <- as.numeric(featuredf$feature1)
  featuredf$feature7 <- as.numeric(featuredf$feature7)
  featuredf$feature8 <- as.numeric(featuredf$feature8)
  return(featuredf)
}

featuredf <- convert.feature.to.num(featuredf)
test.featuredf <- convert.feature.to.num(cluster.test.scoredf)

matchdf <- read.csv(file.path(path.to.main.input, data.version, "cluster_score_old.csv"), sep = ",")[, c("BlockID", "SampleID")]
featuredf <- merge(featuredf, matchdf, by.x = "BlockID", by.y = "BlockID")[, c("SampleID", all.features)]

featuredf.old <- read.csv(file.path(path.to.05.output, "featuredf.raw.csv")) %>% subset(select = -c(X))

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

convert.colname.clusterscore.engvn <- list(
  "Độ.sâu.xâm.lấn..DOI." = "Depth of invasion (DOI)",
  "Độ.sừng.hóa.tế.bào" = "Degree of keratinization",
  "Dị.dạng.nhân" = "Nuclear polymorphism",
  "Kiểu.xâm.lấn" = "Pattern of invasion",
  "Thấm.nhập.lympho.tương.bào" = "Lymphoplasmacytic \r\ninfiltration",
  "Độ.ác.tính.mô.học.theo.Bryne" = "Bryne scores",
  "Phân.độ.mô.học.theo.WHO" = "WHO system",
  "Xâm.nhiễm.thần.kinh..PNI." = "Perineural invasion (PNI)" ,
  "Kiểu.hình.xâm.nhiễm.xấu.nhất..WPOI."  = "Worst pattern of invasion (WPOI)" ,
  "Xâm.nhiễm.lympho..LHR." = "Lymphocytic host response (LHR)",
  "Độ.ác.tính.mô.học.theo.Brandwwein.Gensler" = "BrandweinGensler risk level",
  "Chẩn.đoán.vị.trí" = "Primary tumor site" ,
  "Kích.thước.bướu.lâm.sàng..cT." = "Clinical tumor size (cT)",
  "Di.căn.hạch.lâm.sàng..cN." = "Clinical node mestatasis (cN)",
  "Di.căn.xa" = "Distant metastasis (M)",
  "Giai.đoạn.lâm.sàng..TNM." = "TNM classification",
  "Hút.thuốc" = "Smoking tobacco",
  "Uống.rượu" = "Alcohol drinking",
  "Chẩn.đoán.giải.phẫu.bệnh"= "Pathological diagnosis",
  "Di.căn.hạch.giải.phẫu.bệnh..pN." =  "Pathological lymph node metastasis (pN)",
  "Kiểu.tăng.trưởng...Xem.đại.thể." = "Growth patterns"
)

namedf <- data.frame(vn.name = names(convert.colname.clusterscore))
namedf$old.feature <- unlist(convert.colname.clusterscore)
namedf <- namedf %>% rowwise() %>%
  mutate(eng.name = convert.colname.clusterscore.engvn[[vn.name]]) %>%
  mutate(new.feature = all.features[[eng.name]])

featuredf <- featuredf %>% column_to_rownames("SampleID")
featuredf.old <- featuredf.old%>% column_to_rownames("SampleID")

sample.order <- rownames(featuredf)
feature.order <- to_vec( for (item in colnames(featuredf.old))  subset(namedf, namedf$old.feature == item)$new.feature)
write.csv(featuredf[, feature.order] %>% rownames_to_column("SampleID"), file.path(path.to.06.output, "featuredf.final.csv"))
write.csv(featuredf[, c(feature.order, "feature1")] %>% rownames_to_column("SampleID"), file.path(path.to.06.output, "featuredf.final_with_Gender.csv"))
write.csv(test.featuredf[,  c("No.", feature.order)] , file.path(path.to.06.output, "test_featuredf.final.csv"))
write.csv(test.featuredf[, c("No.", feature.order, "feature1")], file.path(path.to.06.output, "test_featuredf.final_with_Gender.csv"))

# for (f in namedf$old.feature){
#   tmp.old <- featuredf.old[sample.order, f]
#   tmp.new <- featuredf[sample.order, subset(namedf, namedf$old.feature == f)$new.feature]
#   if (sum(tmp.old - tmp.new) !=0){
#     print(f)
#   }
# }

writexl::write_xlsx(namedf, file.path(path.to.06.output, "convert_feature_nam_df.xlsx"))
