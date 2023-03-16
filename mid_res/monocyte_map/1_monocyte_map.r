# R 3.5 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
suppressMessages(library(here))
#source("de_workflow-master/downstream_analysis_functions.r")
source("functions/analysis_functions.R")
btm = readRDS("signature_curation/BTM_li.rds")
datapath = here('mid_res/monocyte_map/generated_data/'); dir.create(datapath)

# select celltype on which to to run pseudotime analysis 
celltype_use = c("CD14_Mono", "CD16_Mono")

# subset time cohort 
sub = ReadCohort(joint_object_dir = "data/h1h5_annotated_with_meta.rds", cohort = "H1N1") %>% 
  SetAllIdent(id = "time_cohort") %>% 
  SubsetData(ident.use  = "d1", subset.raw = TRUE) %>% 
  SetAllIdent(id = "celltype_joint") %>% 
  SubsetData(ident.use = celltype_use, subset.raw = TRUE)

# addd proteins as meta data 
prot_dat = as.data.frame(t(sub@assay$CITE@data))
sub = AddMetaData(sub, metadata = prot_dat)

#### use DDR tree to calculate trajectory 
library(monocle)
sm = monocle::importCDS(sub)
sm = BiocGenerics::estimateSizeFactors(sm)
sm <- detectGenes(sm, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(sm),num_cells_expressed >= 15))

## Select genes 
time1_genes = differentialGeneTest(sm[expressed_genes,],fullModelFormulaStr = "~timepoint",cores = 4)
rpgene =  grep(pattern = "RPL|RPS|MT-|RP11", x = expressed_genes, value = TRUE)
t1genes = time1_genes %>% 
  rownames_to_column("gene") %>% 
  filter(qval < 0.15) %>% 
  arrange(qval) %>% 
  filter(!gene %in% rpgene)
t1gene = t1genes %$% gene

# set ordering filter and reduce dimensions by the ddr tree algorithm 
msm = setOrderingFilter(sm, ordering_genes = t1gene)
sm = reduceDimension(sm, max_components = 2, 
                     reduction_method = "DDRTree", 
                     residualModelFormulaStr = "~ sampleid")
sm = orderCells(sm, reverse = TRUE)
saveRDS(sm, file = paste0(datapath, "sm_cd14_cd16_d1_monocle_object.rds"))
