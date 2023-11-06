# H5 vs H1 day 1 DE contrast pseudobulk model Enrichment
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(variancePartition))
suppressMessages(library(magrittr))
#suppressMessages(library(scglmmr))
source(here('functions/scglmmr.functions.R'))

# specify output directories and init parallel opts
datapath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/")
dir.create(datapath)
figpath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/figuresV4/gsea/")
dir.create(figpath, recursive = TRUE)

# parallel options 
BiocParallel::register(BiocParallel::SnowParam(4))
pparam = BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# load combined modules -- rm the RP gene outliers 
cmod = readRDS(file = here('signature_curation/combined_sig_sub.rds'))

# load contrast fit results 
fit12e = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/fit12e.rds'))
toprank = ExtractResult(
  model.fit.list = fit12e,
  what = 'lmer.z.ranks',
  coefficient.number = 1,
  coef.name = 'delta'
)

# gsea on combined modules
gc = FgseaList(rank.list.celltype = toprank,pathways = cmod, BPPARAM = pparam)
saveRDS(gc,file = paste0(datapath, 'gc.rds'))
gc = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/gc.rds'))

# index leading edge genes 
li = LeadingEdgeIndexed(gsea.result.list = gc,padj.threshold = 0.05)
li = Filter(li,f =  length)
saveRDS(li,file = paste0(datapath, 'li.rds'))

# jaccard enrichment
gsub = lapply(gc, function(x) x %>%  filter(padj < 0.05))
gsub = Filter(gsub, f = nrow)
ji = EnrichmentJaccard(gsealist = gsub, 
                       indexedgenes = li, 
                       saveplot = FALSE, 
                       figpath = figpath, 
                       returnJaccardMtx = TRUE)
saveRDS(ji$sortedgsea,file = paste0(datapath, 'sortedgsea.rds'))

# subet of leading edge genes upregulated in contrast model only. 
g.up = lapply(gc, function(x)
  x %>%  filter(padj < 0.1 & NES > 0)
  )
li.up = LeadingEdgeIndexed(gsea.result.list = g.up, padj.threshold = 0.1)
li.up = Filter(li.up, f = length)
saveRDS(li.up,file = paste0(datapath, 'li.up.rds'))

# add top genes by effect 
res = ExtractResult(model.fit.list = fit12e, coefficient.number = 1, coef.name = 'delta')
topgene = lapply(res, function(x)
  x %>% 
    filter(logFC > 0.25 & P.Value < 0.03) %$% 
    gene
  )
saveRDS(topgene, file = paste0(datapath, 'topgene.rds'))


