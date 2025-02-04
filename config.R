library(hash)
config.params <- hash()

all.down.count.thres <- c(100, 1000, 2000)
input.k <- 3
up.count.thres <- 7000
all.top_variable_genes <- c(500, 1000, 2000)
all.var.method <- c("none", "mad")

count <- 1

for (down.count.thres in all.down.count.thres){
  for (top_variable_genes in all.top_variable_genes){
    for (var.method in all.var.method){
      for (batch.effect in c(TRUE, FALSE)){
        config.params[[sprintf("v%s", count)]] <-  list(
          down.count.thres = down.count.thres,
          up.count.thres = up.count.thres,
          input.k = input.k,
          top_variable_genes = top_variable_genes,
          batch.effect = batch.effect,
          var.method = var.method
        )
        count <- count + 1
      }
    }
  }
}

library(tidyverse)
library(dplyr)

configdf <- data.frame(code.version = names(config.params))

configdf <- configdf %>% rowwise() %>% 
  mutate(down.count.thres = config.params[[code.version]]$down.count.thres) %>%
  mutate(top_variable_genes = config.params[[code.version]]$top_variable_genes) %>%
  mutate(var.method = config.params[[code.version]]$var.method) %>%
  mutate(batch.effect = config.params[[code.version]]$batch.effect)
configdf$input.k <- input.k
configdf$up.count.thres <- up.count.thres

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/UMP_Oral_cancer"
PROJECT <- "UMP_oral_cancer"
output.version <- "output_20240411"
path.to.main.input <- "/media/hieunguyen/GSHD_HN01/raw_data/UMP_Oral_cancer/input"
path.to.main.output <- file.path(outdir, PROJECT, output.version)

write.csv(configdf, file.path(path.to.main.output, "config.csv"))