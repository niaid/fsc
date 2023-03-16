# R version 4.0.5
# H1N1 Unadjuvanted group baseline gene set enrichment analysis
# gsea based on the ranks of the contrast high vs low responder pre vaccinaiton
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))

# set save paths 
datapath = here("mid_res/baseline_response/dataV3/"); dir.create(datapath)
figpath = here("mid_res/baseline_response/figuresV3/"); dir.create(figpath)

# parallel options for FseaList
BiocParallel::register(BiocParallel::SnowParam(4))
pparam = BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# load combined pathways 
mods = readRDS(file = here('signature_curation/combined_signatures.rds'))

# load baseline contrast, rank genes run gsea 
cont0 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV3/cont0.rds'))
r0 = ExtractResult(model.fit.list = cont0, what = 'gene.t.ranks',coefficient.number = 1, coef.name = 'adjmfc')

# run fgea on each cell type
g0 = FgseaList(rank.list.celltype = r0, pathways = mods, BPPARAM = pparam)
saveRDS(object = g0, file = paste0(datapath, 'g0.rds'))
