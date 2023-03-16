# R version 4.0.5 
# H5 vs H1 day 1 DE contrast pseudobulk model 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(variancePartition))
suppressMessages(library(scglmmr))
source("functions/analysis_functions.R")

# make output directories 
datapath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/"); 
dir.create(datapath, recursive = TRUE)
figpath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/figuresV4/"); 
dir.create(figpath, recursive = TRUE)

# parallel options for dream lme4 fits 
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# read metadata to exract sample names in analysis 
samplemd = readRDS(file = here('data/samplemd.rds')) 
d1sx = samplemd %>% filter(time_cohort =='d1') %$% sample

# read processed pseudobulk data 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))

# remove cell type string from sample names and subset to day 1 cohort
cnames = gsub("~.*","",colnames(pb[[1]]))
pb = lapply(pb, function(x){ 
  x %>% as.data.frame() %>% 
    setNames(nm = cnames) %>% 
    select(all_of(d1sx)) %>% 
    as.matrix() 
  })

# sample metadata for contrast model 
samplemd =readRDS(file = here('data/samplemd.rds')) 
samplemd = 
  samplemd %>% 
  filter(! time_cohort == 'd7') %>% 
  mutate(scaledage = (age - mean(age)) / sd(age)) %>% 
  rename('subjectid' = sampleid) %>% 
  rename('group' = adjmfc.group) %>%  
  mutate(group = ifelse(group %in% c('high', 'low'), "NOAS03", "AS03")) %>% 
  mutate(time.group = paste(timepoint, group, sep = "_"))
# get rid of space
samplemd$time.group = str_replace_all(
  string = samplemd$time.group,
  pattern = ' ',
  replacement = ''
)
# re-level time.group into ordered combined factor 
samplemd$time.group = factor(
  samplemd$time.group,
  levels = c('d0_AS03', 'd1_AS03', 'd0_NOAS03', 'd1_NOAS03')
  )
# format metadata 
samplemd = samplemd %>% 
  remove_rownames() %>% 
  column_to_rownames('sample')

# designmat 
met = samplemd[ ,c('gender', 'scaledage', 'time.group')]
mat = model.matrix( ~ 0 + time.group +  gender + scaledage, data = met)

################################
# specify random intercept model formula with combined factor 
f1 <- ~ 0 + time.group + gender + scaledage + (1|subjectid) 

# specify contrast matrix to test the fold change difference 
# based on levels of time.group this should be cmat = c(-1, 1, 1, -1, 0, 0)
L2 = makeContrastsDream(formula = f1, data = samplemd,
  contrasts = c(delta = "(time.groupd1_AS03 - time.groupd0_AS03) - (time.groupd1_NOAS03 - time.groupd0_NOAS03)")
  )
plotContrasts(L2) + ggsave(filename = paste0(figpath,'contrastmodel.pdf'), width = 7, height = 4)

# fit model on each subset 
# init store 
fit1 = v1 = list()
for (i in 1:length(pb)) {
  
  # init data 
  meta = samplemd
  form = f1 
  contrast_matrix = L2
  counts = pb[[i]]
  
  # dge list 
  d = edgeR::DGEList(counts = counts, samples = meta)
  
  # filter cell type specific lowly expressed genes and calc norm factors 
  gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, design = mat)
  print(names(pb)[i]);print(table(gtable))
  d = d[gtable, keep.lib.sizes=FALSE]
  d = edgeR::calcNormFactors(object = d)
  
  
  # get voom observation level weights 
  v = voomWithDreamWeights(counts = d,
                           formula = form,
                           data = meta,
                           BPPARAM = pparam,
                           plot = TRUE, save.plot = TRUE)
  
  # fit contrast mixed model 
  fitmm = dream(exprObj = v, 
                formula = form,
                data = meta,
                L = contrast_matrix,
                BPPARAM = pparam, 
                useWeights = TRUE, REML = TRUE)
  
  # save results 
  v1[[i]] = v
  fit1[[i]] = fitmm
}
names(v1) = names(fit1) = names(pb)

# save day 1 contrast fit 
saveRDS(object = fit1, file = paste0(datapath, 'fit1.rds'))
saveRDS(object = v1, file = paste0(datapath, 'v1.rds'))


# apply eBayes to model fit 
fit1 = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/fit1.rds'))
fit1e  = lapply(fit1, variancePartition::eBayes) 
saveRDS(object = fit1e, file = paste0(datapath, 'fit1e.rds'))


# save model fitting data 
saveRDS(object = samplemd, file = paste0(datapath, 'samplemd.rds'))
saveRDS(object = pb, file = paste0(datapath, 'pb.rds'))
saveRDS(object = L2, file = paste0(datapath, 'L2.rds'))

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] viridis_0.5.1            viridisLite_0.3.0        scglmmr_0.1.0            variancePartition_1.25.6 BiocParallel_1.24.1      limma_3.46.0            
# [7] magrittr_2.0.1           SeuratObject_4.0.0       Seurat_4.0.1             here_1.0.1               forcats_0.5.1            stringr_1.4.0           
# [13] dplyr_1.0.4              purrr_0.3.4              readr_1.4.0              tidyr_1.1.2              tibble_3.0.6             ggplot2_3.3.3           
# [19] tidyverse_1.3.0         
# 
# loaded via a namespace (and not attached):
# [1] estimability_1.3            scattermore_0.7             coda_0.19-4                 bit64_4.0.5                 irlba_2.3.3                
# [6] multcomp_1.4-16             DelayedArray_0.16.3         data.table_1.14.0           rpart_4.1-15                RCurl_1.98-1.3             
# [11] doParallel_1.0.16           generics_0.1.0              snow_0.4-3                  BiocGenerics_0.36.1         RhpcBLASctl_0.21-247.1     
# [16] cowplot_1.1.1               TH.data_1.0-10              RSQLite_2.2.7               shadowtext_0.0.9            RANN_2.6.1                 
# [21] future_1.21.0               bit_4.0.4                   enrichplot_1.10.2           spatstat.data_2.1-0         xml2_1.3.2                 
# [26] lubridate_1.7.9.2           httpuv_1.5.5                SummarizedExperiment_1.20.0 assertthat_0.2.1            hms_1.0.0                  
# [31] promises_1.2.0.1            progress_1.2.2              caTools_1.18.1              dbplyr_2.1.0                readxl_1.3.1               
# [36] igraph_1.2.6                DBI_1.1.1                   htmlwidgets_1.5.3           spatstat.geom_2.0-1         stats4_4.0.5               
# [41] ellipsis_0.3.1              ggpubr_0.4.0                backports_1.2.1             annotate_1.68.0             aod_1.3.1                  
# [46] deldir_0.2-10               MatrixGenerics_1.2.1        vctrs_0.3.6                 Biobase_2.50.0              ROCR_1.0-11                
# [51] abind_1.4-5                 cachem_1.0.4                withr_2.4.1                 ggforce_0.3.3               emmeans_1.5.4              
# [56] sctransform_0.3.2           prettyunits_1.1.1           goftest_1.2-2               cluster_2.1.2               DOSE_3.16.0                
# [61] lazyeval_0.2.2              crayon_1.4.1                labeling_0.4.2              edgeR_3.32.1                pkgconfig_2.0.3            
# [66] tweenr_1.0.2                GenomeInfoDb_1.26.7         nlme_3.1-152                rlang_0.4.10                globals_0.14.0             
# [71] lifecycle_1.0.0             miniUI_0.1.1.1              sandwich_3.0-0              downloader_0.4              modelr_0.1.8               
# [76] cellranger_1.1.0            rprojroot_2.0.2             polyclip_1.10-0             GSVA_1.38.2                 matrixStats_0.58.0         
# [81] lmtest_0.9-38               graph_1.68.0                Matrix_1.3-2                carData_3.0-4               boot_1.3-27                
# [86] zoo_1.8-8                   reprex_1.0.0                pheatmap_1.0.12             ggridges_0.5.3              png_0.1-7                  
# [91] bitops_1.0-6                KernSmooth_2.23-18          blob_1.2.1                  qvalue_2.22.0               parallelly_1.23.0          
# [96] rstatix_0.7.0               S4Vectors_0.28.1            ggsignif_0.6.0              scales_1.1.1                memoise_2.0.0              
# [101] GSEABase_1.52.1             plyr_1.8.6                  ica_1.0-2                   gplots_3.1.1                zlibbioc_1.36.0            
# [106] compiler_4.0.5              scatterpie_0.1.7            RColorBrewer_1.1-2          lme4_1.1-26                 fitdistrplus_1.1-3         
# [111] cli_2.5.0                   XVector_0.30.0              lmerTest_3.1-3              listenv_0.8.0               patchwork_1.1.1            
# [116] pbapply_1.4-3               MASS_7.3-53.1               mgcv_1.8-34                 tidyselect_1.1.0            stringi_1.5.3              
# [121] GOSemSim_2.16.1             locfit_1.5-9.4              ggrepel_0.9.1               grid_4.0.5                  fastmatch_1.1-0            
# [126] tools_4.0.5                 rio_0.5.16                  future.apply_1.7.0          parallel_4.0.5              rstudioapi_0.13            
# [131] foreign_0.8-81              foreach_1.5.1               gridExtra_2.3               farver_2.0.3                Rtsne_0.15                 
# [136] ggraph_2.0.5                digest_0.6.27               rvcheck_0.1.8               BiocManager_1.30.10         shiny_1.6.0                
# [141] Rcpp_1.0.6                  car_3.0-10                  GenomicRanges_1.42.0        broom_0.7.5                 egg_0.4.5                  
# [146] later_1.1.0.1               RcppAnnoy_0.0.18            org.Hs.eg.db_3.12.0         httr_1.4.2                  AnnotationDbi_1.52.0       
# [151] Rdpack_2.1.1                colorspace_2.0-0            rvest_0.3.6                 XML_3.99-0.6                fs_1.5.0                   
# [156] tensor_1.5                  reticulate_1.18             IRanges_2.24.1              splines_4.0.5               uwot_0.1.10                
# [161] statmod_1.4.35              spatstat.utils_2.1-0        graphlayouts_0.7.2          plotly_4.9.3                xtable_1.8-4               
# [166] jsonlite_1.7.2              nloptr_1.2.2.2              tidygraph_1.2.0             ggfun_0.0.4                 R6_2.5.0                   
# [171] pillar_1.4.7                htmltools_0.5.1.1           mime_0.10                   glue_1.4.2                  fastmap_1.1.0              
# [176] minqa_1.2.4                 clusterProfiler_3.18.1      codetools_0.2-18            fgsea_1.16.0                mvtnorm_1.1-1              
# [181] lattice_0.20-41             spatstat.sparse_2.0-0       numDeriv_2016.8-1.1         pbkrtest_0.5-0.1            curl_4.3                   
# [186] leiden_0.3.7                gtools_3.8.2                zip_2.1.1                   openxlsx_4.2.3              GO.db_3.12.1               
# [191] survival_3.2-10             munsell_0.5.0               DO.db_2.9                   GenomeInfoDbData_1.2.4      iterators_1.0.13           
# [196] haven_2.3.1                 reshape2_1.4.4              gtable_0.3.0                rbibutils_2.0               spatstat.core_2.0-0
