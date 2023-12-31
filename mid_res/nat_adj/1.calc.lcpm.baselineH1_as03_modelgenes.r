# R version 4.0.5 
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(scglmmr))
source(here('functions/MattPMutils.r'))

# set save paths 
figpath = here("mid_res/nat_adj/figures/V4/"); dir.create(figpath)
datapath = here("mid_res/nat_adj/generated_data//V4/"); dir.create(datapath)

# define high responders 
high.responders = c("205","207","209","212","215","234","237","245","250","256")

# read pb data, subset to day 0 non adj, subset out day 0 metadata. 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))
pb = lapply(pb, function(x) x = x[ ,1:40])
cnames = gsub("~.*","",colnames(pb[[1]]))
pb = lapply(pb, function(x){
  x %>% as.data.frame() %>% setNames(nm = cnames) %>% as.matrix() 
})
d0 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV3/d0.rds'))
d0d = lapply(pb, function(x){ x = x[ , rownames(d0)]})

# make a list of genes indxed by celltype for genes to fit from H5 model 
av_tidy = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/git_ignore/av_tidy.rds'))
genes.test = lapply(av_tidy , function(x) unique(x$gene))

# check names are the same of indexed genes and average data
stopifnot(isTRUE(all.equal(names(d0d), names(genes.test))))

# calcualte log CPM of the baseline pseudobulk data
av = list()
for (i in 1:length(d0d)) {
  dge =  edgeR::DGEList( d0d[[i]] ) 
  dge = dge[genes.test[[i]], ]
  av[[i]]  = edgeR::cpm(dge, log = TRUE)
}
names(av) = names(d0d)
saveRDS(av,file = paste0(datapath,'av.rds'))

# tidy aggregated data 
av0 = list()
for (i in 1:length(pb)) {
  ct = names(av)[i]
  gs = rownames(av[[i]])
  av0[[i]] = GetTidySummary(
    av.exprs.list = av,
    celltype.index = i,
    genes.use = gs) %>%
    mutate(response = if_else(str_sub(sample, 1, 3) %in%  high.responders, 'High', "Low")) %>%
    mutate(response = factor(response, levels = c('Low', 'High')))
}
names(av0) = names(pb)
saveRDS(av0, file = paste0(datapath,'av0.rds'))
