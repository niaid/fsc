# average gene distributions - part 1 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
#suppressMessages(library(scglmmr))
source(here('functions/MattPMutils.r'))
source('functions/scglmmr.functions.R')

datapath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gene_dist/");
dir.create(datapath)

# calculate logcpm 
pb = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/pb12.rds'))
samplemd = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/samplemd12.rds'))
time.group = factor(samplemd$time.group)
av = list()
for (i in 1:length(pb)) {
  dge =  edgeR::DGEList( pb[[i]] ) 
  gtable = edgeR::filterByExpr(y = dge$counts, min.count = 3,  design = time.group)
  dge = dge[gtable, ]
  av[[i]]  = edgeR::cpm(dge, log = TRUE)
}
names(av) = names(pb)

# create tidy data for gene visualization and calculation of average signature scores 
# get tidy summary data 
av_tidy = list()
for (i in 1:length(av)) {
  ct = names(av)[i]
  gs = rownames(av[[i]])
  av_tidy[[i]] = GetTidySummary(av.exprs.list = av, 
                                celltype.index = i,
                                genes.use = gs)  %>% 
    mutate(cohort = if_else( str_sub(sample, 1,2) == "H5", "H5N1", "H1N1")) %>% 
    mutate(group = paste(str_sub(sample, -2,-1), cohort))
  av_tidy[[i]]$group= factor(av_tidy[[i]]$group, 
                             levels = 
                               c("d0 H1N1" ,"d1 H1N1" ,"d0 H5N1","d1 H5N1"))
  av_tidy[[i]]$group = plyr::revalue(av_tidy[[i]]$group, 
                                     c("d0 H1N1" = "d0 No-AS03",
                                       "d1 H1N1" = "d1 No-AS03",
                                       "d0 H5N1" = "d0 AS03",
                                       "d1 H5N1" = "d1 AS03"))
}
names(av_tidy) = names(av)
saveRDS(av_tidy, file = paste0(datapath,'av_tidy.rds'))


