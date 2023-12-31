suppressMessages(library(here))
suppressMessages(library(tidyverse))
source(here('functions/scglmmr.functions.R'))
suppressMessages(library(magrittr))
source(here('functions/MattPMutils.r'))

# set save paths 
figpath = here("mid_res/nat_adj/figures/V4/")
datapath = here("mid_res/nat_adj/generated_data//V4/")

# load mono and mDC AS03 specific signatures 
mono.genes = readRDS(file = here('mid_res/combined_contrast/generated_data/mo.noifn.rds'))
mdc.genes = readRDS(file = here('mid_res/combined_contrast/generated_data/dc.noifn.rds'))

# extract leading edge into signatures 
mono.genes = mono.genes$leadingEdge %>% unlist() %>% unique()
mdc.genes = mdc.genes$leadingEdge %>% unlist() %>% unique()
li.full = c(mono.genes, mdc.genes) %>% unique()
as03.sig = list('AS03_signature' = li.full, "AS03_Monocyte" = mono.genes, 'AS03_mDC' = mdc.genes)
as03.sig.list = as03.sig
saveRDS(as03.sig.list, file = paste0(datapath, 'as03.sig.list.rds'))

# validation cohort for the highlighted signatures
vand.fit = readRDS(file = here('mid_res/vand/generated_data/fit1.rds'))
vand.rank = ExtractResult(model.fit.list = vand.fit,what = 'lmer.z.ranks', coefficient.number = 1, coef.name = 'delta')

# Enrichment of AS03 signatures without ifn signatures 
gvand = FgseaList(rank.list.celltype = vand.rank, pathways = as03.sig)
saveRDS(gvand,file = paste0(datapath,'gvand.rds'))

gvand$MNC
# pathway         pval         padj  log2err        ES      NES size                                  leadingEdge celltype
# 1: AS03_signature 1.910373e-34 5.731119e-34 1.529705 0.5714156 3.144152  263 SERPINA1,SECTM1,HCK,PLSCR1,LILRA1,FCGR1B,...      MNC
# 2:  AS03_Monocyte 4.563863e-26 6.845794e-26 1.326716 0.5979691 3.088277  170 SERPINA1,SECTM1,HCK,PLSCR1,LILRA1,FCGR2A,...      MNC
# 3:       AS03_mDC 6.087321e-21 6.087321e-21 1.186651 0.5916694 2.984150  145 SERPINA1,SECTM1,HCK,LILRA1,FCGR1B,FCGR2A,...      MNC

gvand$DNC
#       pathway         pval         padj   log2err        ES      NES size                                 leadingEdge celltype
# 1: AS03_signature 4.154660e-17 1.246398e-16 1.0672100 0.4932650 2.402933  257 PSME2,FPR1,SLC31A2,PSMB10,PSME1,ALDH1A1,...      DNC
# 2:       AS03_mDC 1.761296e-12 2.641944e-12 0.9101197 0.5311486 2.420236  146 PSME2,FPR1,SLC31A2,PSMB10,PSME1,ALDH1A1,...      DNC
# 3:  AS03_Monocyte 6.171525e-10 6.171525e-10 0.8012156 0.4774659 2.201396  163   FPR1,SLC31A2,DPYD,IFIH1,LILRB3,PLSCR1,...      DNC

# leading edge from each innate cell types - AS03 signature 
dc.as03.sig.validated = gvand$DNC %>%  filter(pathway == 'AS03_mDC') %$% leadingEdge %>%  unlist()
mono.as03.sig.validated = gvand$MNC %>%  filter(pathway == 'AS03_Monocyte') %$% leadingEdge %>%  unlist() 
saveRDS(dc.as03.sig.validated, file = paste0(datapath,'dc.as03.sig.validated.rds'))
saveRDS(mono.as03.sig.validated, file = paste0(datapath,'mono.as03.sig.validated.rds'))

data.table::fwrite(list(dc.as03.sig.validated),file = paste0(datapath,'dc.as03.sig.validated.txt'),sep = '\t')
data.table::fwrite(list(mono.as03.sig.validated),file = paste0(datapath,'mono.as03.sig.validated.txt'),sep = '\t')

# plot enrichment distributions 
enrline = list(geom_line(color = "deepskyblue3", size = 2 ))
#mono
p = fgsea::plotEnrichment(pathway = as03.sig$AS03_Monocyte, stats = vand.rank$MNC) + enrline
ggsave(p, filename = paste0(figpath, 'mono.vand.enr.2.pdf'), width = 5, height = 3)
# mDC
p = fgsea::plotEnrichment(pathway = as03.sig$AS03_mDC, stats = vand.rank$DNC) + enrline
ggsave(p, filename = paste0(figpath, 'dc.vand.enr.2.pdf'), width = 5, height = 3)








