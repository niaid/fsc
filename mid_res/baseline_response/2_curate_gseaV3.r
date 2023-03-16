# R version 4.0.5
# Curate H1N1 Unadjuvanted group baseline gene set enrichment analysis
# calculate pairwise jaccard index and reduce enrichments to major signals with low mutual information
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))

# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")

# read baseline enrichemnt results 
mrm = readRDS(file = here('signature_curation/module_rmlist.rds'))
g0 = readRDS(file = here('mid_res/baseline_response/dataV3/g0.rds'))
g0 = lapply(g0, function(x) x %>% filter(!pathway %in% mrm))
filtered_g0 = lapply(g0, function(x) x %>% filter(padj < 0.05))

# compute jaccard index of leadingedge genes within celltype  
li = LeadingEdgeIndexed(gsea.result.list = filtered_g0, padj.threshold = 0.05)

jres = EnrichmentJaccard(gsealist = filtered_g0, indexedgenes = li, 
                         saveplot = FALSE,
                         figpath = figpath,
                         returnJaccardMtx = TRUE, 
                         fontsize_row = 7.5, fontsize_col = 7.5)
d = jres$sortedgsea %>% 
  mutate(leadingEdge = map_chr(leadingEdge, toString)) %>% 
  select(celltype, av_jaccard,everything())
write_delim(d,file = paste0(datapath, 'g0jaccard.csv'),delim = ',')

# save the jaccard matrices 
jmats = jres$jaccard_matrix_list
saveRDS(jmats ,file = paste0(datapath, 'jmats.rds'))

sessionInfo()