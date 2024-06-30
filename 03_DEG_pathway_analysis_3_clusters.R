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
library("sva")

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
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

all.names.3clusters <- c("kmean.3clusters.DrNam",
                         "RNA.consensus.cluster")

for (cluster.name in all.names.3clusters){
  path.to.main.input <- "/media/hieunguyen/HNSD_mini/data/UMP_Oral_cancer/input"
  path.to.main.output <- file.path(outdir, PROJECT, output.version)
  dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)
  
  path.to.01.output <- file.path(path.to.main.output, "01_output", data.version, code.version)
  path.to.03.output <- file.path(path.to.main.output, "03_output", data.version, code.version, cluster.name)
  dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.to.03.output, "DESEQ_results"), showWarnings = FALSE, recursive = TRUE)
  
  dir.create(file.path(path.to.03.output, "pathway_results_ORA"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.to.03.output, "pathway_plots_ORA"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.to.03.output, "pathway_objs_ORA"), showWarnings = FALSE, recursive = TRUE)
  
  dir.create(file.path(path.to.03.output, "pathway_results_GSEA"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.to.03.output, "pathway_plots_GSEA"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.to.03.output, "pathway_objs_GSEA"), showWarnings = FALSE, recursive = TRUE)
  
  dir.create(file.path(path.to.03.output, "count_matrix"), showWarnings = FALSE, recursive = TRUE)
  
  writexl::write_xlsx(as.data.frame(config.params[[code.version]]), file.path(path.to.03.output, "config.params.logs.xlsx"))
  
  umapdf.rna <- read.csv(file.path(path.to.01.output, "umap_RNAseq.csv"))
  
  #####----------------------------------------------------------------------#####
  ##### DIFFERENTIAL GENE EXPRESSION ANALYSIS
  #####----------------------------------------------------------------------#####
  kallisto.obj <- readRDS(file.path(path.to.main.input, data.version, "txi.kallisto.RDS"))
  all.samples <- colnames(kallisto.obj$counts)
  meta.data <- read.csv(file.path(path.to.main.input, data.version, "Metadata.txt"), sep = "\t")
  meta.data <- subset(meta.data, meta.data$SampleID %in% all.samples)
  meta.data <- merge(meta.data, umapdf.rna, by.x = "SampleID", by.y = "SampleID")
  
  #####----------------------------------------------------------------------#####
  ##### SUBSET KALLISTO OBJECT
  #####----------------------------------------------------------------------#####
  run_deseq <- function(cluster1, cluster2, kallisto.obj, input.metadata, cluster.name){
    filtered.metadata <- subset(input.metadata, input.metadata[[cluster.name]] %in% c(cluster1, cluster2))
    subset.kallisto.obj <- kallisto.obj
    for (name in names(subset.kallisto.obj)){
      if (name != "countsFromAbundance"){
        tmp <- subset.kallisto.obj[[name]]
        subset.kallisto.obj[[name]] <- tmp[, filtered.metadata$SampleID]
      }
    }
    filtered.metadata$condition <- factor(filtered.metadata[[cluster.name]], levels = c(cluster1, cluster2))
    subset.kallisto.obj[["countsFromAbundance"]] <- kallisto.obj$countsFromAbundance
    deseq.obj <- DESeqDataSetFromTximport(txi = subset.kallisto.obj, colData = filtered.metadata, design = ~condition)
    deseq.obj <- estimateSizeFactors(deseq.obj)
    deseq.obj <- deseq.obj[rowSums(counts(deseq.obj)) > down.count.thres, ]
    dds <- DESeq(deseq.obj)
    return(dds)  
  }
  
  dds.obj <- hash()
  
  dds.obj[["cluster1_vs_cluster2"]] <- run_deseq(cluster1 = 1, cluster2 = 2, kallisto.obj = kallisto.obj, input.metadata = meta.data, cluster.name = cluster.name)
  dds.obj[["cluster1_vs_cluster3"]] <- run_deseq(cluster1 = 1, cluster2 = 3, kallisto.obj = kallisto.obj, input.metadata = meta.data, cluster.name = cluster.name)
  dds.obj[["cluster2_vs_cluster3"]] <- run_deseq(cluster1 = 2, cluster2 = 3, kallisto.obj = kallisto.obj, input.metadata = meta.data, cluster.name = cluster.name)
  
  saveRDS(dds.obj, file.path(path.to.03.output, "DESEQ2_RNAseq_data_one_vs_one.obj.rds"))
  
  ##### save normalized matrix
  for (i in names(dds.obj)){
    save.count.matrix <- counts(dds.obj[[i]], normalized = TRUE)
    phenotypedf <- dds.obj[[i]]@colData %>% as.data.frame() %>% subset(select = c(SampleID, RNA.consensus.cluster))
    write.table(save.count.matrix, file.path(path.to.03.output, "count_matrix", sprintf("count_matrix_%s.csv", i)), sep = "\t")
    write.table(phenotypedf, file.path(path.to.03.output, "count_matrix", sprintf("phenotype_%s.csv", i)), sep = "\t")
  }
  
  all.diff.genes <- hash()
  all.raw.diff.genes <- hash()
  
  for (i in names(dds.obj)){
    tmpdf <- dds.obj[[i]] %>% results() %>% as.data.frame() %>% rownames_to_column("Gene")
    all.raw.diff.genes[[i]] <- tmpdf
    tmpdf <- tmpdf %>% subset(padj <= 0.05) %>%
      rowwise() %>%
      mutate(abs.logFC = abs(log2FoldChange)) %>%
      arrange(desc(abs.logFC))
    writexl::write_xlsx(tmpdf, file.path(path.to.03.output, "DESEQ_results", sprintf("DESEQ_result_%s.xlsx", i)))
    all.diff.genes[[i]] <- tmpdf
  }
  
  #####----------------------------------------------------------------------#####
  ##### PATHWAY ANALYSIS HELPER FUNCTIONS
  #####----------------------------------------------------------------------#####
  generate_all_pathway_analysis_ORA <- function(cluster1, cluster2){
    raw.gene.list <- all.raw.diff.genes[[sprintf("cluster%s_vs_cluster%s", cluster1, cluster2)]]
    raw.gene.list <- bitr(raw.gene.list$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    sig.gene.list.up <- all.diff.genes[[sprintf("cluster%s_vs_cluster%s", cluster1, cluster2)]] %>% subset(log2FoldChange > 0) %>% subset(select = c(Gene))
    sig.gene.list.down <- all.diff.genes[[sprintf("cluster%s_vs_cluster%s", cluster1, cluster2)]] %>% subset(log2FoldChange <= 0) %>% subset(select = c(Gene))
    
    if (nrow(all.diff.genes[[sprintf("cluster%s_vs_cluster%s", cluster1, cluster2)]] ) > 1){
      sig.gene.list.up <- bitr(sig.gene.list.up$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      sig.gene.list.down <- bitr(sig.gene.list.down$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      all.enrich.res <- hash()
      all.enrich.df <- hash()
      
      ##### GO - up
      enrich.go <- enrichGO(
        gene = sig.gene.list.up$ENTREZID,
        OrgDb = org.Hs.eg.db,
        universe = raw.gene.list$ENTREZID,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE) 
      
      all.enrich.res[["GO_up"]] <- enrich.go
      enrich.go <- enrich.go %>% as.data.frame()
      enrich.go <- subset(enrich.go, select = c(ID, Description, GeneRatio, BgRatio, p.adjust, Count)) %>%
        arrange(desc(Count))
      all.enrich.df[["GO_up"]] <- enrich.go
      
      ##### GO - down
      enrich.go <- enrichGO(
        gene = sig.gene.list.down$ENTREZID,
        OrgDb = org.Hs.eg.db,
        universe = raw.gene.list$ENTREZID,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE) 
      
      all.enrich.res[["GO_down"]] <- enrich.go
      enrich.go <- enrich.go %>% as.data.frame()
      enrich.go <- subset(enrich.go, select = c(ID, Description, GeneRatio, BgRatio, p.adjust, Count)) %>%
        arrange(desc(Count))
      all.enrich.df[["GO_down"]] <- enrich.go
      
      ##### KEGG - up
      enrich.KEGG <- enrichKEGG(gene = sig.gene.list.up$ENTREZID,
                                organism     = 'hsa',
                                pvalueCutoff = 0.05)
      
      all.enrich.res[["KEGG_up"]] <- enrich.KEGG
      enrich.KEGG <- enrich.KEGG %>% as.data.frame()
      enrich.KEGG <- subset(enrich.KEGG, select = c(ID, category, subcategory, Description, GeneRatio, BgRatio, p.adjust, Count)) %>%
        arrange(desc(Count))
      all.enrich.df[["KEGG_up"]] <- enrich.KEGG
      
      ##### KEGG - down
      enrich.KEGG <- enrichKEGG(gene = sig.gene.list.down$ENTREZID,
                                organism     = 'hsa',
                                pvalueCutoff = 0.05)
      
      all.enrich.res[["KEGG_down"]] <- enrich.KEGG
      enrich.KEGG <- enrich.KEGG %>% as.data.frame()
      enrich.KEGG <- subset(enrich.KEGG, select = c(ID, category, subcategory, Description, GeneRatio, BgRatio, p.adjust, Count)) %>%
        arrange(desc(Count))
      all.enrich.df[["KEGG_down"]] <- enrich.KEGG
      
      ##### msigdbr enrichment UP
      all.cats <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
      for (selectedC in all.cats){
        m_t2g <- msigdbr(species = "Homo sapiens", category = selectedC) %>% 
          dplyr::select(gs_name, entrez_gene, gs_description)
        em <- enricher(sig.gene.list.up$ENTREZID, TERM2GENE=m_t2g)
        all.enrich.res[[sprintf("%s_up", selectedC)]] <- em
        em <- em %>% as.data.frame()
        if ("geneID" %in% colnames(em)){
          em <- subset(em, select = -c(geneID)) %>%
            arrange(desc(Count))          
        }
        all.enrich.df[[sprintf("%s_up", selectedC)]] <- em
      }
      
      ##### msigdbr enrichment DOWN
      all.cats <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
      for (selectedC in all.cats){
        m_t2g <- msigdbr(species = "Homo sapiens", category = selectedC) %>% 
          dplyr::select(gs_name, entrez_gene, gs_description)
        em <- enricher(sig.gene.list.down$ENTREZID, TERM2GENE=m_t2g)
        all.enrich.res[[sprintf("%s_down", selectedC)]] <- em
        em <- em %>% as.data.frame()
        if ("geneID" %in% colnames(em)){
          em <- subset(em, select = -c(geneID)) %>%
            arrange(desc(Count))          
        }
        all.enrich.df[[sprintf("%s_down", selectedC)]] <- em
      }
      
      return(list(res = all.enrich.res, dfs = all.enrich.df))
    } else {
      return(NULL)
    }
  }
  
  generate_results_from_enricher_ORA <- function(cluster1, cluster2){
    output <- generate_all_pathway_analysis_ORA(cluster1 = cluster1, cluster2 = cluster2)
    saveRDS(output, file.path(path.to.03.output, "pathway_objs_ORA", sprintf("pathway_obj_cluster%s_vs_cluster%s.rds", cluster1, cluster2)))
    for (i in names(output$dfs)){
      tmp.table <- output$dfs[[i]]
      if (nrow(tmp.table) != 0){
        print(sprintf("Saving results from pathway %s", i))
        writexl::write_xlsx(tmp.table, file.path(path.to.03.output, "pathway_results_ORA", sprintf("enriched_pathway_%s.xlsx", i)))
      }
    }
    
    for (i in names(output$res)){
      if (nrow(output$res[[i]]) != 0){
        print(sprintf("Saving results from pathway %s", i))
        p <- dotplot(output$res[[i]], showCategory=15) + ggtitle(sprintf("Top 20 pathways in gene set %s", i)) +
          theme(axis.text.x = element_text(size = 22),
                axis.text.y = element_text(size = 22),
                axis.title.x = element_text(size = 22),
                axis.title.y = element_text(size = 22),
                title = element_text(size = 22))    
        ggsave(plot = p, filename = sprintf("dotplot_pval_geneRatio_pathway_%s.svg", i), 
               path = file.path(path.to.03.output, "pathway_plots_ORA"), 
               device = "svg", 
               width = 290, 
               height = 290, units = "px")
      } else {
        p <- ggplot() + ggtitle("No pathway to show")
      }
      print(p)
    }
  }
  
  
  generate_all_pathway_analysis_GSEA <- function(cluster1, cluster2){
    all.enrich.res <- hash()
    all.enrich.df <- hash()
    
    raw.gene.list <- all.raw.diff.genes[[sprintf("cluster%s_vs_cluster%s", cluster1, cluster2)]]
    convertdf <- bitr(raw.gene.list$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    raw.gene.list <- merge(raw.gene.list, convertdf, by.x = "Gene", by.y = "SYMBOL") 
    raw.gene.list <- raw.gene.list %>% arrange(desc(log2FoldChange))
    gene.list <- raw.gene.list$log2FoldChange 
    names(gene.list) <- raw.gene.list$ENTREZID
    
    ##### GO
    enrich.go <- gseGO(geneList = gene.list,
                       OrgDb = org.Hs.eg.db,
                       minGSSize = 100,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       verbose = FALSE)
    
    all.enrich.res[["GO"]] <- enrich.go
    enrich.go <- enrich.go %>% as.data.frame() %>% rowwise() %>% mutate(absNES = abs(NES)) %>% arrange(desc(absNES))
    # enrich.go <- subset(enrich.go, select = -c(core_enrichment, leading_edge))
    all.enrich.df[["GO"]] <- enrich.go
    
    ##### KEGG
    enrich.KEGG <- gseKEGG(geneList     = gene.list,
                           organism     = 'hsa',
                           minGSSize    = 120,
                           pvalueCutoff = 0.05,
                           verbose      = FALSE)
    
    all.enrich.res[["KEGG"]] <- enrich.KEGG
    enrich.KEGG <- enrich.KEGG %>% as.data.frame() %>% rowwise() %>% mutate(absNES = abs(NES)) %>% arrange(desc(absNES))
    # enrich.KEGG <- subset(enrich.KEGG, select = -c(core_enrichment, leading_edge))
    
    all.enrich.df[["KEGG"]] <- enrich.KEGG
    
    ##### msigdbr enrichment
    all.cats <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
    for (selectedC in all.cats){
      m_t2g <- msigdbr(species = "Homo sapiens", category = selectedC) %>% 
        dplyr::select(gs_name, entrez_gene, gs_description)
      em <- GSEA(gene.list, TERM2GENE=m_t2g)
      all.enrich.res[[selectedC]] <- em
      em <- em %>% as.data.frame() %>% rowwise() %>% mutate(absNES = abs(NES)) %>% arrange(desc(absNES))
      if ("geneID" %in% colnames(em)){
        # em <- subset(em, select = -c(core_enrichment, leading_edge)) 
      }
      
      all.enrich.df[[selectedC]] <- em
    }
    return(list(res = all.enrich.res, dfs = all.enrich.df))
  }
  
  generate_results_from_enricher_GSEA <- function(cluster1, cluster2){
    output <- generate_all_pathway_analysis_GSEA(cluster1 = cluster1, cluster2 = cluster2)
    saveRDS(output, file.path(path.to.03.output, "pathway_objs_GSEA", sprintf("pathway_obj_cluster%s_vs_cluster%s.rds", cluster1, cluster2)))
    for (i in names(output$dfs)){
      tmp.table <- output$dfs[[i]]
      if (nrow(tmp.table) != 0){
        print(sprintf("Saving results from pathway %s", i))
        writexl::write_xlsx(tmp.table, file.path(path.to.03.output, "pathway_results_GSEA", sprintf("enriched_pathway_%s.xlsx", i)))
      }
    }
    
    for (i in names(output$res)){
      if (nrow(output$res[[i]]) != 0){
        print(sprintf("Saving results from pathway %s", i))
        p <- ridgeplot(output$res[[i]], showCategory=15, orderBy = "NES") + ggtitle(sprintf("Top 20 pathways in gene set %s", i))    
        ggsave(plot = p, filename = sprintf("ridge_plot_pathway_%s.svg", i), path = file.path(path.to.03.output, "pathway_plots_GSEA"), device = "svg", width = 290, height = 290, units = "px")
      } else {
        p <- ggplot() + ggtitle("No pathway to show")
      }
      print(p)
    }
  }
  
  #####----------------------------------------------------------------------#####
  ##### MAIN RUN
  #####----------------------------------------------------------------------#####
  generate_results_from_enricher_ORA(cluster1 = 1, cluster2 = 2)
  generate_results_from_enricher_ORA(cluster1 = 1, cluster2 = 3)
  generate_results_from_enricher_ORA(cluster1 = 2, cluster2 = 3)
  
  generate_results_from_enricher_GSEA(cluster1 = 1, cluster2 = 2)
  generate_results_from_enricher_GSEA(cluster1 = 1, cluster2 = 3)
  generate_results_from_enricher_GSEA(cluster1 = 2, cluster2 = 3)
}  
