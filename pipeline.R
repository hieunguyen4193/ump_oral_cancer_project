path.to.main.src <- "/media/hieunguyen/HNSD01/src/UMP_oral_cancer/official/ump_oral_cancer_project"

source(file.path(path.to.main.src, "install_cluster_profiler.R"))

source(file.path(path.to.main.src, "01_preprocess_RNA_seq_dataset.R"))
source(file.path(path.to.main.src, "02_DEG_pathway_analysis_2_clusters.R"))
source(file.path(path.to.main.src, "03_DEG_pathway_analysis_3_clusters.R"))
source(file.path(path.to.main.src, "04_preprocess_score_features.R"))
source(file.path(path.to.main.src, "06_plot.R"))

