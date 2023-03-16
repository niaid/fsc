# Early kinetics of baseline states 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))
# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")


# load single cell data 
h1 = readRDS(file = here('data/h1h5_annotated_with_meta.rds'))
md = h1@meta.data %>% 
  filter(cohort == 'H1N1') %>% 
  filter(time_cohort == 'd1') %>% 
  mutate(group_id = factor(adjmfc.time, levels = c("d0 low", "d1 low", "d0 high", "d1 high"))) %>% 
  mutate(subjectid = factor(sampleid)) %>% 
  mutate(sex = factor(gender)) %>% 
  mutate(scaled.age = as.numeric(scale(age))) %>% 
  mutate(celltype = celltype_joint) 
  
# add a covariate for number of cells per sample  
ncell_met = md %>% group_by(sample) %>% summarize(ncells = n())
md$ncell = plyr::mapvalues(x = md$sample, from = ncell_met$sample, to = ncell_met$ncells)
md$ncell = as.numeric(md$ncell)
md$log10ncell = log10(md$ncell) 

# subset normalized RNA data 
norm.rna = h1@data[ ,rownames(md)]


# get leading edge genes from cur. baseline mods 
g0.sub = readRDS(file = here('mid_res/baseline_response/dataV3/g0.sub.rds'))
li.g0 = LeadingEdgeIndexed(gsea.result.list = g0.sub, padj.threshold = 0.05)
li.g0 = base::Filter(length, li.g0)

# metadata by cell type 
cts = names(li.g0)
md = md %>% filter( celltype %in% cts )
ct.md = split( md, f = md$celltype )

# module score for each cell type of specific baseline enriched leading edge genes.
mod_scores = list()
for (i in 1:length(ct.md)) {
  
  # init data for subset i 
  rna = norm.rna[ ,rownames(ct.md[[i]])]
  mod.list = li.g0[[i]]
  
  # calculate single cell score for baseline-enriched module 
  mod_scores[[i]] = WeightedCellModuleScore(
    gene_matrix = rna, 
    module_list = mod.list, 
    threshold = 0, 
    cellwise_scaling = TRUE, 
    return_weighted = FALSE 
    )
  
  # add a "null" score of Gaussian noise as a reference 
  mod_scores[[i]]$null = rnorm(n = nrow(mod_scores[[i]]), mean = 0, sd = 1)
}

# specify save paths for marginal means plots.
plot_savepath1 = paste0(figpath, "/marginalmeans.m1/"); dir.create(plot_savepath1)
plot_savepath2 = paste0(figpath, "/marginalmeans.m2/"); dir.create(plot_savepath2)


# specify the 2 models 
f1 = 'modulescore ~ 0 + group_id + (1|subjectid)'
f2 = 'modulescore ~ 0 + group_id + log10ncell + scaled.age + sex + (1|subjectid)'


# fit sc mod mixed model on module scores. 
mm1 = mm2 = list()
for (i in 1:length(ct.md)) {
  
  stopifnot( nrow(ct.md[[i]]) == nrow(mod_scores[[i]]) )
  
  # formula 1
  mm1[[i]] = FitLmerContrast(module_data_frame = mod_scores[[i]], 
                              celltype_column = 'celltype', 
                              metadata = ct.md[[i]], 
                              lmer_formula = f1, 
                              plotdatqc = FALSE, 
                              fixed_effects = NULL,
                              figpath = plot_savepath1)
  
  # formula 2 
  mm2[[i]] = FitLmerContrast(module_data_frame = mod_scores[[i]], 
                             celltype_column = 'celltype', 
                             metadata = ct.md[[i]], 
                             lmer_formula = f2, 
                             plotdatqc = FALSE, 
                             fixed_effects = NULL,
                             figpath = plot_savepath2)  

}

mm1 = do.call(rbind, mm1)
mm2 = do.call(rbind, mm2)

saveRDS(mm1,file = paste0(datapath, 'mm1.rds'))
saveRDS(mm2,file = paste0(datapath, 'mm2.rds'))


sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] SeuratObject_4.0.0 Seurat_4.0.1       scglmmr_0.1.0      magrittr_2.0.1     here_1.0.1        
# [6] forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0       
# [11] tidyr_1.1.2        tibble_3.0.6       ggplot2_3.3.3      tidyverse_1.3.0   
# 
# loaded via a namespace (and not attached):
# [1] estimability_1.3            scattermore_0.7             coda_0.19-4                
# [4] bit64_4.0.5                 irlba_2.3.3                 multcomp_1.4-16            
# [7] DelayedArray_0.16.3         rpart_4.1-15                data.table_1.14.0          
# [10] RCurl_1.98-1.3              generics_0.1.0              BiocGenerics_0.36.1        
# [13] cowplot_1.1.1               TH.data_1.0-10              RSQLite_2.2.7              
# [16] shadowtext_0.0.9            RANN_2.6.1                  future_1.21.0              
# [19] bit_4.0.4                   enrichplot_1.10.2           spatstat.data_2.1-0        
# [22] xml2_1.3.2                  lubridate_1.7.9.2           httpuv_1.5.5               
# [25] SummarizedExperiment_1.20.0 assertthat_0.2.1            viridis_0.5.1              
# [28] hms_1.0.0                   promises_1.2.0.1            caTools_1.18.1             
# [31] dbplyr_2.1.0                readxl_1.3.1                igraph_1.2.6               
# [34] DBI_1.1.1                   htmlwidgets_1.5.3           spatstat.geom_2.0-1        
# [37] stats4_4.0.5                ellipsis_0.3.1              ggpubr_0.4.0               
# [40] backports_1.2.1             annotate_1.68.0             deldir_0.2-10              
# [43] MatrixGenerics_1.2.1        vctrs_0.3.6                 Biobase_2.50.0             
# [46] ROCR_1.0-11                 abind_1.4-5                 cachem_1.0.4               
# [49] withr_2.4.3                 ggforce_0.3.3               packrat_0.7.0              
# [52] emmeans_1.5.4               sctransform_0.3.2           goftest_1.2-2              
# [55] cluster_2.1.2               DOSE_3.16.0                 lazyeval_0.2.2             
# [58] crayon_1.4.1                edgeR_3.32.1                pkgconfig_2.0.3            
# [61] labeling_0.4.2              tweenr_1.0.2                GenomeInfoDb_1.26.7        
# [64] nlme_3.1-152                rlang_0.4.10                globals_0.14.0             
# [67] lifecycle_1.0.0             miniUI_0.1.1.1              sandwich_3.0-0             
# [70] downloader_0.4              modelr_0.1.8                cellranger_1.1.0           
# [73] rprojroot_2.0.2             polyclip_1.10-0             GSVA_1.38.2                
# [76] matrixStats_0.58.0          lmtest_0.9-38               graph_1.68.0               
# [79] Matrix_1.3-2                carData_3.0-4               boot_1.3-27                
# [82] zoo_1.8-8                   reprex_1.0.0                ggridges_0.5.3             
# [85] pheatmap_1.0.12             png_0.1-7                   viridisLite_0.3.0          
# [88] bitops_1.0-6                KernSmooth_2.23-18          blob_1.2.1                 
# [91] qvalue_2.22.0               parallelly_1.23.0           rstatix_0.7.0              
# [94] S4Vectors_0.28.1            ggsignif_0.6.0              scales_1.1.1               
# [97] memoise_2.0.0               GSEABase_1.52.1             plyr_1.8.6                 
# [100] ica_1.0-2                   gplots_3.1.1                zlibbioc_1.36.0            
# [103] compiler_4.0.5              scatterpie_0.1.7            RColorBrewer_1.1-2         
# [106] lme4_1.1-26                 fitdistrplus_1.1-3          cli_2.5.0                  
# [109] XVector_0.30.0              listenv_0.8.0               patchwork_1.1.1            
# [112] pbapply_1.4-3               mgcv_1.8-34                 MASS_7.3-53.1              
# [115] tidyselect_1.1.0            stringi_1.5.3               GOSemSim_2.16.1            
# [118] locfit_1.5-9.4              ggrepel_0.9.1               GeneOverlap_1.26.0         
# [121] grid_4.0.5                  fastmatch_1.1-0             tools_4.0.5                
# [124] future.apply_1.7.0          parallel_4.0.5              rio_0.5.16                 
# [127] rstudioapi_0.13             foreign_0.8-81              gridExtra_2.3              
# [130] farver_2.0.3                Rtsne_0.15                  ggraph_2.0.5               
# [133] digest_0.6.27               rvcheck_0.1.8               BiocManager_1.30.10        
# [136] shiny_1.6.0                 Rcpp_1.0.6                  GenomicRanges_1.42.0       
# [139] car_3.0-10                  broom_0.7.5                 egg_0.4.5                  
# [142] later_1.1.0.1               RcppAnnoy_0.0.18            org.Hs.eg.db_3.12.0        
# [145] httr_1.4.2                  AnnotationDbi_1.52.0        colorspace_2.0-0           
# [148] tensor_1.5                  rvest_0.3.6                 XML_3.99-0.6               
# [151] fs_1.5.0                    reticulate_1.18             IRanges_2.24.1             
# [154] splines_4.0.5               uwot_0.1.10                 statmod_1.4.35             
# [157] spatstat.utils_2.1-0        graphlayouts_0.7.2          plotly_4.9.3               
# [160] xtable_1.8-4                jsonlite_1.7.2              nloptr_1.2.2.2             
# [163] tidygraph_1.2.0             ggfun_0.0.4                 R6_2.5.0                   
# [166] pillar_1.4.7                htmltools_0.5.2             mime_0.10                  
# [169] glue_1.4.2                  fastmap_1.1.0               minqa_1.2.4                
# [172] clusterProfiler_3.18.1      BiocParallel_1.24.1         codetools_0.2-18           
# [175] fgsea_1.16.0                mvtnorm_1.1-1               spatstat.sparse_2.0-0      
# [178] lattice_0.20-41             slanter_0.2-0               curl_4.3                   
# [181] leiden_0.3.7                gtools_3.8.2                zip_2.1.1                  
# [184] GO.db_3.12.1                openxlsx_4.2.3              survival_3.2-10            
# [187] limma_3.46.0                munsell_0.5.0               DO.db_2.9                  
# [190] GenomeInfoDbData_1.2.4      haven_2.3.1                 reshape2_1.4.4             
# [193] gtable_0.3.0                spatstat.core_2.0-0        



