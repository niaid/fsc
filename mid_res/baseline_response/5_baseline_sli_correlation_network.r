# baseline nettwork correlation 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
source('functions/MattPMutils.r')

# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")

################################
# Part I correct intracellular correlations for overlapping gene content
################################

# pairwise jaccard index of all leading edge genes from modules 
library(scglmmr)
g0.sub = readRDS(file = here("mid_res/baseline_response/dataV3/g0.sub.rds"))
li.index = scglmmr::LeadingEdgeIndexed(gsea.result.list = g0.sub,padj.threshold = 0.05)

# remove lists with no enrichments indexed by cell type 
g0.sub = base::Filter(g0.sub, f = nrow)
li.index = base::Filter(li.index, f = length)

# reorder enrichment to match leading edge index
g0.sub = g0.sub[names(li.index)]
stopifnot(isTRUE(all.equal(names(g0.sub), names(li.index))))

ss = li.index[c('CD14_Mono', 'CD16_Mono', 'mDC', 'MAIT_Like')]
ss = lapply(ss, function(x) x %>%  unlist() %>% unname() %>% unique())


VennDiagram::venn.diagram(ss, 
                          # color 
                          fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"), 
                          # text
                          cat.cex = 0.4,cat.default.pos = "inner",cat.fontfamily = "Helvetica",
                          # file 
                          imagetype="png", height = 1080, width = 1080, resolution = 300,compression = "lzw",
                          filename = paste0(figpath,"gpvenn.png")
)
go = VennDiagram::calculate.overlap(ss)
go <- unlist(
  lapply(
    1:length(ss), function(j) {
           combn(names(ss), j, simplify = FALSE)
           }),
  recursive = FALSE
  )
names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
str(combs)


# get jaccard index matrix of intercellular leading edge genes 
# from different enrichments 
jmat = EnrichmentJaccard(gsealist = g0.sub, 
                         indexedgenes = li.index, 
                         returnJaccardMtx = TRUE)
jmat = jmat$jaccard_matrix_list

# load of baseline module expression across donors only of 
# leading edge genes from baseline enrichments 
ds = readRDS(here('mid_res/baseline_response/dataV3/ds.rds'))
data.table::fwrite(ds,file = paste0(here('git_ignore/ds.csv')),sep = ',')

# created shorter names 
new.names = data.table::fread(
  file = here('mid_res/baseline_response/dataV3/baseline.module.name.shortened.txt'),
  sep = '\t'
  )
new.names$cname2 = paste(new.names$celltype, new.names$shortname, sep = ' :: ')
colnames(ds) = new.names$cname2

# split module expression by cell types to subtract jacard similariry index 
# from intracellular spearman correlation coefficient caluclated below 
# also add shorter names names 
ds2 = ds %>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column('cname2') %>% 
  full_join(new.names, by = "cname2") %>% 
  select(-c('module', 'cname', 'shortname')) %>% 
  select(celltype, everything())
# split
ds2 = split(ds2,f = ds2$celltype)

# calculate Spearman correlation for intracellular correlations 
ds.cor = lapply(ds2, function(x) { 
  mtr =  x %>%  
    select(-c('celltype')) %>% 
    remove_rownames() %>% 
    column_to_rownames('cname2') %>% 
    t() 
  return(Hmisc::rcorr(mtr, type = 'spearman')$r)
  })

# subset to celltypes with multiple enrichments 
ds.cor = ds.cor[names(jmat)]

# append the names of the jaccard matrix with cell type so that format is 
# celltype :: module name to match correlation matrix
for (i in 1:length(jmat)) {
  mt = jmat[[i]]
  colnames(mt) = rownames(mt) = plyr::mapvalues(
    rownames(mt), 
    from = new.names$module, 
    to = new.names$shortname
    )
  colnames(mt) = rownames(mt) = paste(names(jmat)[i], colnames(mt),sep = ' :: ')
  jmat[[i]] = mt
}

# calculate shared latent informaiton of intracellular correlations
# subtract jaccard similariry from the spearman correlation coefficient 
# diagonal corrects to 0 
sli = list()
for (i in 1:length(jmat)) {
  stopifnot(isTRUE(all.equal(rownames(ds.cor[[i]]), rownames(jmat[[i]]))) )
  sli[[i]] = ds.cor[[i]] - jmat[[i]]
}

# now calculate the full spearman correlation matrix without subtracting jmat 
# replace the matrix values from intracellular correlations with the SLI values 
# which are stored in sli list by celltype 
spearmanmat = Hmisc::rcorr(ds, type = 'spearman')
saveRDS(spearmanmat, file = paste0(datapath, 'spearmanmat.rds'))
spearmanmat = readRDS(file = here('mid_res/baseline_response/dataV3/spearmanmat.rds'))
rhomat = spearmanmat$r
mat = rhomat
for (i in 1:length(sli)) {
  # ged index of cols and rows 
  row.replace = which(rownames(mat) %in% rownames(sli[[i]]))
  col.replace = which(colnames(mat) %in% colnames(sli[[i]]))
  
  # we need to be replacing values along the square diagonal of the matrix 
  # confirm these are the same 
  stopifnot(isTRUE(all.equal(row.replace, col.replace)))
  
  # check structure 
  stopifnot(isTRUE(all.equal(
    # full spearman matrix subset by rows of celltype i 
    mat[row.replace, col.replace],
    # original spearman correlation matrix for celltype i 
    ds.cor[[i]]
    )))
  
  # replace iteratively 
  mat[row.replace,col.replace] = sli[[i]]
  
  # this should be going down with each iteration 
  print(sum(mat))

}

# save SLI corrected matrix 
saveRDS(mat,file = paste0(datapath,'mat.rds'))
mat = readRDS(file = here('mid_res/baseline_response/dataV3/mat.rds'))

################################
# Part II visualization
################################

# plot pre and post adjustment correlation map 
# spearmancorrelation coefficient rhomat
cu = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(11)
range <- max(abs(rhomat))
# without clustering
dev.off()
pdf(file = paste0(figpath,'preadjusted.cormat.baseline.pdf'),width = 6, height = 5)
pheatmap::pheatmap(rhomat, color = cu, 
                   cluster_rows = FALSE, cluster_cols = FALSE, 
                   border_color = NA,
                   breaks = seq(-range, range, length.out = 11),
                   fontsize_row = 5, fontsize_col = 0.01)
dev.off()

# plot sli corrected matrix 
pdf(file = paste0(figpath,'post.SLIadjusted.cormat.baseline.pdf'),width = 6, height = 5)
pheatmap::pheatmap(mat, color = cu, 
                   cluster_rows = FALSE, cluster_cols = FALSE, 
                   border_color = NA,
                   breaks = seq(-range, range, length.out = 11),
                   fontsize_row = 5, fontsize_col = 0.01)
dev.off()

# figure  -- full network visualization as matrix (non pruned)
# cluster the square matrix 
#diag(mat)[diag(mat) > 0] = 0 
range <- max(abs(mat))
cu3 = BuenColors::jdb_palette('solar_flare', type = 'continuous') %>%  as.vector
cu3 = cu3[seq(from = 0 , to = 1000,length.out = 20)]
rownames(mat) = str_replace_all(string = rownames(mat),pattern = '_',replacement = ' ')
pdf(file = paste0(figpath,'post.clustered.SLIadjusted.cormat.baseline.pdf'),width = 8, height = 6.5)
pheatmap::pheatmap(mat, 
                   color = cu3, 
                   #cluster_rows = FALSE, cluster_cols = FALSE, 
                   border_color = NA,
                   treeheight_row = 10, 
                   treeheight_col = 20,
                   breaks = seq(-range, range, length.out = 20),
                   fontsize_row = 5.5,
                   fontsize_col = 0.01)
dev.off()


# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scglmmr_0.1.0   magrittr_2.0.1  here_1.0.1      forcats_0.5.1   stringr_1.4.0   dplyr_1.0.4     purrr_0.3.4    
# [8] readr_1.4.0     tidyr_1.1.2     tibble_3.0.6    ggplot2_3.3.3   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] lme4_1.1-26                 tidyselect_1.1.0            RSQLite_2.2.7               AnnotationDbi_1.52.0       
# [5] htmlwidgets_1.5.3           grid_4.0.5                  BiocParallel_1.24.1         scatterpie_0.1.7           
# [9] munsell_0.5.0               codetools_0.2-18            statmod_1.4.35              withr_2.4.3                
# [13] colorspace_2.0-0            GOSemSim_2.16.1             Biobase_2.50.0              knitr_1.39                 
# [17] rstudioapi_0.13             stats4_4.0.5                ggsignif_0.6.0              DOSE_3.16.0                
# [21] Rdpack_2.1.1                MatrixGenerics_1.2.1        emmeans_1.5.4               GenomeInfoDbData_1.2.4     
# [25] polyclip_1.10-0             bit64_4.0.5                 farver_2.0.3                pheatmap_1.0.12            
# [29] rprojroot_2.0.2             downloader_0.4              coda_0.19-4                 vctrs_0.4.1                
# [33] generics_0.1.2              TH.data_1.0-10              xfun_0.30                   doParallel_1.0.16          
# [37] R6_2.5.0                    GenomeInfoDb_1.26.7         graphlayouts_0.7.2          locfit_1.5-9.4             
# [41] pals_1.7                    bitops_1.0-6                cachem_1.0.4                fgsea_1.16.0               
# [45] DelayedArray_0.16.3         assertthat_0.2.1            scales_1.1.1                multcomp_1.4-16            
# [49] ggraph_2.0.5                nnet_7.3-15                 enrichplot_1.10.2           gtable_0.3.0               
# [53] egg_0.4.5                   tidygraph_1.2.0             sandwich_3.0-0              rlang_1.0.2                
# [57] slanter_0.2-0               splines_4.0.5               rstatix_0.7.0               dichromat_2.0-0            
# [61] broom_0.7.5                 checkmate_2.0.0             BiocManager_1.30.10         reshape2_1.4.4             
# [65] abind_1.4-5                 modelr_0.1.8                backports_1.2.1             qvalue_2.22.0              
# [69] Hmisc_4.5-0                 clusterProfiler_3.18.1      tools_4.0.5                 ellipsis_0.3.2             
# [73] gplots_3.1.1                RColorBrewer_1.1-2          BiocGenerics_0.36.1         Rcpp_1.0.6                 
# [77] plyr_1.8.6                  progress_1.2.2              base64enc_0.1-3             zlibbioc_1.36.0            
# [81] RCurl_1.98-1.3              prettyunits_1.1.1           ggpubr_0.4.0                rpart_4.1-15               
# [85] viridis_0.5.1               cowplot_1.1.1               S4Vectors_0.28.1            zoo_1.8-8                  
# [89] SummarizedExperiment_1.20.0 haven_2.3.1                 ggrepel_0.9.1               cluster_2.1.2              
# [93] fs_1.5.0                    variancePartition_1.25.6    data.table_1.14.0           DO.db_2.9                  
# [97] openxlsx_4.2.3              reprex_1.0.0                mvtnorm_1.1-1               packrat_0.7.0              
# [101] matrixStats_0.58.0          hms_1.0.0                   GSVA_1.38.2                 xtable_1.8-4               
# [105] pbkrtest_0.5-0.1            RhpcBLASctl_0.21-247.1      XML_3.99-0.6                rio_0.5.16                 
# [109] jpeg_0.1-8.1                BuenColors_0.5.6            readxl_1.3.1                IRanges_2.24.1             
# [113] gridExtra_2.3               compiler_4.0.5              maps_3.4.0                  KernSmooth_2.23-18         
# [117] crayon_1.4.1                shadowtext_0.0.9            minqa_1.2.4                 htmltools_0.5.2            
# [121] ggfun_0.0.4                 Formula_1.2-4               lubridate_1.7.9.2           DBI_1.1.1                  
# [125] tweenr_1.0.2                dbplyr_2.1.0                MASS_7.3-53.1               boot_1.3-27                
# [129] Matrix_1.3-2                car_3.0-10                  cli_3.3.0                   rbibutils_2.0              
# [133] parallel_4.0.5              igraph_1.2.6                GenomicRanges_1.42.0        pkgconfig_2.0.3            
# [137] rvcheck_0.1.8               foreign_0.8-81              foreach_1.5.1               xml2_1.3.2                 
# [141] annotate_1.68.0             XVector_0.30.0              GeneOverlap_1.26.0          estimability_1.3           
# [145] rvest_0.3.6                 digest_0.6.27               graph_1.68.0                cellranger_1.1.0           
# [149] fastmatch_1.1-0             htmlTable_2.1.0             edgeR_3.32.1                GSEABase_1.52.1            
# [153] curl_4.3                    gtools_3.8.2                nloptr_1.2.2.2              nlme_3.1-152               
# [157] lifecycle_1.0.0             jsonlite_1.7.2              aod_1.3.1                   carData_3.0-4              
# [161] mapproj_1.2.8               viridisLite_0.3.0           limma_3.46.0                pillar_1.4.7               
# [165] lattice_0.20-41             fastmap_1.1.0               httr_1.4.2                  survival_3.2-10            
# [169] GO.db_3.12.1                glue_1.6.2                  zip_2.1.1                   iterators_1.0.13           
# [173] png_0.1-7                   bit_4.0.4                   ggforce_0.3.3               stringi_1.5.3              
# [177] blob_1.2.1                  org.Hs.eg.db_3.12.0         latticeExtra_0.6-29         caTools_1.18.1             
# [181] memoise_2.0.0 