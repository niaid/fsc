# gene set nenrichment signals in vand cohort 
# R 4.0.5 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(variancePartition))
suppressMessages(library(scglmmr))

###### save paths  
datapath = here("mid_res/vand/generated_data/")

# parallel opts
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# load signatures 
li.up = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li.up.rds'))
li2.up = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li2.up.rds'))
li = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li.rds'))
#li.full = unlist(li.up,use.names = FALSE, recursive = TRUE) %>%  unique 
mono.combined  = li.up$CD14_Mono %>%  unlist() %>% unique()
mdc.combined = li2.up$mDC %>% unlist() %>% unique()
nb.combined = li.up$BC_Naive %>% unlist() %>% unique()
t.combined = list(
  'Tcell.combined' = c(
    li.up$CD4_CD161_Mem_Tcell,
    li.up$CD4_CD25_Tcell,
    li.up$CD4_Efct_Mem_Tcell,
    li.up$CD4Naive_Tcell,
    li.up$CD8_Mem_Tcell,
    li.up$CD8_Naive_Tcell
  ) %>%
    unlist() %>%
    unique()
)
  

# add additional b cell signatures from apoptosis hypothesis 
sig.test = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/bsig_data/sig.test.rds'))

# add combined signals 
li.up$CD14_Mono$combined.signature = mono.combined
li2.up$mDC$combined.signature = mdc.combined
li$BC_Naive$combined.signature = nb.combined
li$BC_Naive = c(li$BC_Naive, sig.test)

# load vand fits and extract ranks 
fit1 = readRDS(file = here('mid_res/vand/generated_data/fit1.rds'))
vand.rank = ExtractResult(model.fit.list = fit1,
                          what = 'lmer.z.ranks', 
                          coefficient.number = 1, 
                          coef.name = 'delta')

# CD14 monocyte test in total sorted monocyte
mv = FgseaList(
  rank.list.celltype = list('MNC' = vand.rank$MNC),
  pathways = li.up$CD14_Mono,
  BPPARAM = pparam
)

# mDC test in sorted DC
dcv = FgseaList(
  rank.list.celltype = list('DNC' = vand.rank$DNC),
  pathways = li2.up$mDC,
  BPPARAM = pparam
)

# naive BC test in sorted total B
bcv = FgseaList(
  rank.list.celltype = list('BCL' = vand.rank$BCL),
  pathways = li$BC_Naive,
  BPPARAM = pparam
)


# T cell combined in sorted T celsl 
tcv = FgseaList(
  rank.list.celltype = list('TCL' = vand.rank$TCL),
  pathways = c(
    li.up$CD4_CD161_Mem_Tcell,
    li.up$CD4_CD25_Tcell,
    li.up$CD4_Efct_Mem_Tcell,
    li.up$CD4Naive_Tcell,
    li.up$CD8_Mem_Tcell,
    li.up$CD8_Naive_Tcell,
    t.combined
  ),
  BPPARAM = pparam
)

# save objects
saveRDS(object = mv,file = paste0(datapath, 'mv.rds'))
saveRDS(object = dcv,file = paste0(datapath, 'dcv.rds'))
saveRDS(object = bcv,file = paste0(datapath, 'bcv.rds'))
saveRDS(object = tcv,file = paste0(datapath, 'tcv.rds'))