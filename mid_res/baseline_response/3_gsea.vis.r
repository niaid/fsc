# R version 4.0.5
# visualization of gene set enrichment results. 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))
# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")

# res text gsea curated 
d = data.table::fread(here('mid_res/baseline_response/dataV3/g0jaccard.curated.txt')) %>% 
  filter(include ==1) %>% 
  mutate(signal = paste(celltype, pathway, sep = '~')) 

# gsea res raw 
g0 = readRDS(file = here('mid_res/baseline_response/dataV3/g0.rds'))

# filter to the gene sets from curated results. 
g0.sub = list()
for (i in 1:length(g0)) {
  g0.sub[[i]] =
    g0[[i]] %>% 
    mutate(signal = paste(celltype, pathway, sep = '~')) %>% 
    filter(signal %in% d$signal) %>% 
    mutate(celltype = str_replace_all(celltype,pattern = '_',replacement = ' '))
}
names(g0.sub) = names(g0)
# plot
p = PlotFgsea(gsea_result_list = g0.sub, padj_filter = 0.01)
ggsave(p,filename = paste0(figpath, 'gsea.g0sub.baseline.pdf'), width = 9.5, height = 6)
# save object test curated fgsea formatted results. 
saveRDS(g0.sub, file = paste0(datapath, 'g0.sub.rds'))



