# R version 4.0.5
# H1N1 differential expression testing and gene set enrichment analysis
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))

# parallel options for FseaList
BiocParallel::register(BiocParallel::SnowParam(4))
pparam = BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# set data path 
datapath = here("mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/")
dir.create(datapath)

# load pathways to be tested. 
sig_test = readRDS(file = here('signature_curation/combined_sig_sub.rds'))
core_sigs = readRDS(file = here('signature_curation/sig_test_sub.rds'))


# load each time statistical contrast model result extract contrast and rank genes by t statistic 
fit1 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/fit1.rds'))
fit7 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/fit7.rds'))
r1 = ExtractResult(model.fit.list = fit1, what = 'lmer.z.ranks', coefficient.number = 1, coef.name = 'L1')
r7 = ExtractResult(model.fit.list = fit7, what = 'lmer.z.ranks', coefficient.number = 1, coef.name = 'L1')

# run gene set enrichment Day 1 models 
# run unbiased modules and core signatures from past flu studies
g1c = FgseaList(rank.list.celltype = r1, pathways = core_sigs, BPPARAM = pparam)
g1f = FgseaList(rank.list.celltype = r1, pathways = sig_test, BPPARAM = pparam)

# day 7 
g7f = FgseaList(rank.list.celltype = r7, pathways = sig_test, BPPARAM = pparam)

# save 
saveRDS(object = g1c, file = paste0(datapath,'g1c.rds'))
saveRDS(object = g1f, file = paste0(datapath,'g1f.rds'))
saveRDS(object = g7f, file = paste0(datapath,'g7f.rds'))


sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
# [1] viridis_0.5.1            viridisLite_0.3.0        scglmmr_0.1.0            variancePartition_1.25.6
# [5] BiocParallel_1.24.1      limma_3.46.0             magrittr_2.0.1           here_1.0.1
# [9] forcats_0.5.1            stringr_1.4.0            dplyr_1.0.4              purrr_0.3.4
# [13] readr_1.4.0              tidyr_1.1.2              tibble_3.0.6             ggplot2_3.3.3
# [17] tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
# [1] tidyselect_1.1.0            lme4_1.1-26                 RSQLite_2.2.7               AnnotationDbi_1.52.0
# [5] grid_4.0.5                  scatterpie_0.1.7            munsell_0.5.0               codetools_0.2-18
# [9] statmod_1.4.35              withr_2.4.3                 colorspace_2.0-0            GOSemSim_2.16.1
# [13] Biobase_2.50.0              rstudioapi_0.13             stats4_4.0.5                ggsignif_0.6.0
# [17] DOSE_3.16.0                 labeling_0.4.2              MatrixGenerics_1.2.1        Rdpack_2.1.1
# [21] emmeans_1.5.4               GenomeInfoDbData_1.2.4      polyclip_1.10-0             pheatmap_1.0.12
# [25] bit64_4.0.5                 farver_2.0.3                rprojroot_2.0.2             downloader_0.4
# [29] coda_0.19-4                 vctrs_0.4.1                 generics_0.1.2              TH.data_1.0-10
# [33] R6_2.5.0                    doParallel_1.0.16           GenomeInfoDb_1.26.7         graphlayouts_0.7.2
# [37] locfit_1.5-9.4              bitops_1.0-6                cachem_1.0.4                fgsea_1.16.0
# [41] DelayedArray_0.16.3         assertthat_0.2.1            scales_1.1.1                multcomp_1.4-16
# [45] ggraph_2.0.5                enrichplot_1.10.2           gtable_0.3.0                egg_0.4.5
# [49] tidygraph_1.2.0             sandwich_3.0-0              rlang_1.0.2                 slanter_0.2-0
# [53] splines_4.0.5               rstatix_0.7.0               broom_0.7.5                 abind_1.4-5
# [57] BiocManager_1.30.10         reshape2_1.4.4              modelr_0.1.8                backports_1.2.1
# [61] qvalue_2.22.0               clusterProfiler_3.18.1      tools_4.0.5                 ellipsis_0.3.2
# [65] gplots_3.1.1                RColorBrewer_1.1-2          BiocGenerics_0.36.1         Rcpp_1.0.6
# [69] plyr_1.8.6                  progress_1.2.2              zlibbioc_1.36.0             RCurl_1.98-1.3
# [73] prettyunits_1.1.1           ggpubr_0.4.0                cowplot_1.1.1               S4Vectors_0.28.1
# [77] zoo_1.8-8                   SummarizedExperiment_1.20.0 haven_2.3.1                 ggrepel_0.9.1
# [81] fs_1.5.0                    data.table_1.14.0           DO.db_2.9                   openxlsx_4.2.3
# [85] reprex_1.0.0                mvtnorm_1.1-1               packrat_0.7.0               matrixStats_0.58.0
# [89] hms_1.0.0                   GSVA_1.38.2                 xtable_1.8-4                pbkrtest_0.5-0.1
# [93] RhpcBLASctl_0.21-247.1      XML_3.99-0.6                rio_0.5.16                  readxl_1.3.1
# [97] IRanges_2.24.1              gridExtra_2.3               compiler_4.0.5              KernSmooth_2.23-18
# [101] crayon_1.4.1                shadowtext_0.0.9            minqa_1.2.4                 ggfun_0.0.4
# [105] snow_0.4-3                  lubridate_1.7.9.2           DBI_1.1.1                   tweenr_1.0.2
# [109] dbplyr_2.1.0                MASS_7.3-53.1               boot_1.3-27                 Matrix_1.3-2
# [113] car_3.0-10                  cli_3.3.0                   rbibutils_2.0               parallel_4.0.5
# [117] igraph_1.2.6                GenomicRanges_1.42.0        pkgconfig_2.0.3             rvcheck_0.1.8
# [121] foreign_0.8-81              xml2_1.3.2                  foreach_1.5.1               annotate_1.68.0
# [125] XVector_0.30.0              GeneOverlap_1.26.0          estimability_1.3            rvest_0.3.6
# [129] digest_0.6.27               graph_1.68.0                cellranger_1.1.0            fastmatch_1.1-0
# [133] edgeR_3.32.1                GSEABase_1.52.1             curl_4.3                    gtools_3.8.2
# [137] nloptr_1.2.2.2              lifecycle_1.0.0             nlme_3.1-152                jsonlite_1.7.2
# [141] aod_1.3.1                   carData_3.0-4               pillar_1.4.7                lattice_0.20-41
# [145] fastmap_1.1.0               httr_1.4.2                  survival_3.2-10             GO.db_3.12.1
# [149] glue_1.6.2                  zip_2.1.1                   iterators_1.0.13            bit_4.0.4
# [153] ggforce_0.3.3               stringi_1.5.3               blob_1.2.1                  org.Hs.eg.db_3.12.0
# [157] caTools_1.18.1              memoise_2.0.0