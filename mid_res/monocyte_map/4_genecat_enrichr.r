# R 4.0.5
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
library(enrichR)
figpath = here("mid_res/monocyte_map/figures/"); dir.create(figpath, recursive = TRUE)
datapath = here("mid_res/monocyte_map/generated_data/"); dir.create(datapath, recursive = TRUE)


# monocyte leading edge d1 
g1c = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/g1c.rds'))
g1c = lapply(g1c, function(x) x %>%  filter(NES > 0))
mono.le = LeadingEdgeIndexed(gsea.result.list = g1c,padj.threshold = 0.05)


# branch ependent genes 
de.branch = readRDS(file = here('mid_res/monocyte_map/data/de_branch.rds'))
de.branch.sub = de.branch %>% filter(qval < 0.05)
branch.genes = as.character(de.branch.sub$gene_short_name)


# enrichr
dbs <- c("GO_Molecular_Function_2015",
         "GO_Cellular_Component_2015",
         "GO_Biological_Process_2015")

cat2.ifn =  c('IFITM2', 'PTPN1', 'EIF4E2', 'IFITM3', 'HLA-C')
cat1.ifn = intersect(mono.le$CD14_Mono$`reactome interferon signaling`, branch.genes)
cat1.ifn = setdiff(cat1.ifn, cat2.ifn)


# mtor 
cat2.mtor  =c('CTSC', 'PFKL', 'ACTR3', 'CITED2', 'PGK1', 'INSIG1', 'CHST2')
cat1.mtor = intersect(mono.le$CD14_Mono$`HALLMARK MTORC1 signaling`, branch.genes)
cat1.mtor = setdiff(cat1.mtor, cat2.mtor)

# Cat 1 Mtor 
# https://maayanlab.cloud/Enrichr/enrich?dataset=eb3028063e4833deb6212e24235dffee

# Cat 2 mtor 
# https://maayanlab.cloud/Enrichr/enrich?dataset=ace0e94a6d3c1fa0037be65c0f677aa2

# Cat 1 IFN 
# https://maayanlab.cloud/Enrichr/enrich?dataset=537c52fe1a6621033587e9048bf20e98

# Cat 2 ifn 
# https://maayanlab.cloud/Enrichr/enrich?dataset=e3ef3178aeb695f05f90ae855b80da21



