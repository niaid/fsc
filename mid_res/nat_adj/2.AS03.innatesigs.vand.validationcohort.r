# R version 4.0.5 
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(scglmmr))
suppressMessages(library(magrittr))
source(here('functions/MattPMutils.r'))

# set save paths 
figpath = here("mid_res/nat_adj/figures/V4/")
datapath = here("mid_res/nat_adj/generated_data//V4/")

# define adjuvant signatures 
# leading edge combined signature 
# upregulated genes only 
# mono and mDC specific and combined genes
li.up = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li.up.rds'))
li2.up = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li2.up.rds'))
li.full = unlist(li.up,use.names = FALSE, recursive = TRUE) %>%  unique 
mono.genes  = li.up$CD14_Mono %>%  unlist() %>% unique()
mdc.genes = li2.up$mDC %>%  unlist() %>% unique()
as03.sig = list('AS03_signature' = li.full, "AS03_Monocyte" = mono.genes, 'AS03_mDC' = mdc.genes)

# vand validation for the highlighted signatures
vand.fit = readRDS(file = here('mid_res/vand/generated_data/fit1.rds'))
vand.rank = ExtractResult(model.fit.list = vand.fit,what = 'lmer.z.ranks', coefficient.number = 1, coef.name = 'delta')
gvand = FgseaList(rank.list.celltype = vand.rank, pathways = as03.sig)

# leading edge from each cell tpe 
dc.va.le = gvand$DNC %>%  filter(pathway == 'AS03_mDC') %$% leadingEdge %>%  unlist()
saveRDS(dc.va.le,file = paste0(datapath, 'dc.va.le.rds'))

mono.va.le = gvand$MNC %>%  filter(pathway == 'AS03_Monocyte') %$% leadingEdge %>%  unlist() 
saveRDS(mono.va.le, file = paste0(datapath, 'mono.va.le.rds'))


#mono
enrline = list(geom_line(color = "deepskyblue3", size = 2 ))
p = fgsea::plotEnrichment(pathway = as03.sig$AS03_Monocyte, stats = vand.rank$MNC) + enrline
# mDC
p = fgsea::plotEnrichment(pathway = as03.sig$AS03_mDC, stats = vand.rank$DNC) + enrline
ggsave(p, filename = paste0(figpath, 'dc.vand.enr.pdf'), width = 4, height = 2.5)
