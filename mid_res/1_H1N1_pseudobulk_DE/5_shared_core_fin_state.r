# R version 4.0.5
# H1N1 differential expression testing and gene set enrichment analysis
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
source('functions/MattPMutils.r')
# set fig path 
figpath = here("mid_res/1_H1N1_pseudobulk_DE/figuresV4/")

# load model fitting data and pseudobulk data saved in mixed model workflow (scipt 1)
d1d = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/d1d.rds'))

# log CPM the pseudobulk data 
lcpm = lapply(d1d, edgeR::cpm, log = TRUE)  

# read Core signature
ms1.names = readRDS(file = here("mid_res/1_H1N1_pseudobulk_DE/dataV4/ms1.names.rds"))

av2 = lapply(lcpm, function(x){ x[ms1.names, ] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    mutate(timepoint = str_sub(sample, -2, -1))
})

# for (i in 1:length(av2)) {
#   av2[[i]]$celltype = names(av2)[i]
# }
av2 = av2 %>% bind_rows(.id = 'celltype')
av3 = av2 %>% gather(gene, average, ADAR:ZBP1)
av4 = av3 %>% group_by(sample, celltype, timepoint) %>%
  summarize(core_ifn_score = mean(average)) 

av4 = av4 %>% mutate(subject = str_sub(sample, 1,3))
av4$celltype = str_replace_all(string = av4$celltype,pattern = '_', replacement = ' ')

# colors 
cu1 = unname(sapply(c("grey", "orange"), col.alpha, 0.5))
cu2 = c("grey", "orange")
p=
  ggplot(av4, aes(x = celltype, y = core_ifn_score ,  color = timepoint, fill = timepoint )) +
  theme_bw() +
  geom_boxplot() +
  #scale_x_discrete(position = "top") +
  ylab("shared day 1 induced IFN state") +
  xlab("")+
  scale_color_manual(values = cu2 ) +
  scale_fill_manual(values = cu1 ) +
  theme(axis.title.y  = element_text(size = 12, color = 'black' )) +
  theme(axis.text.x = element_text(size = 11, color = "black") )  +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, color = 'black'))+ 
  theme(legend.position = 'top')
p
ggsave(p, filename = paste0(figpath, "core_interferon_state.pdf"), width = 6, height = 5)



sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# 
# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ComplexHeatmap_2.6.2 scglmmr_0.1.0        here_1.0.1           forcats_0.5.1        stringr_1.4.0       
# [6] dplyr_1.0.4          purrr_0.3.4          readr_1.4.0          tidyr_1.1.2          tibble_3.0.6        
# [11] ggplot2_3.3.3        tidyverse_1.3.0     
# 
# loaded via a namespace (and not attached):
# [1] tidyselect_1.1.0            lme4_1.1-26                 RSQLite_2.2.7              
# [4] AnnotationDbi_1.52.0        BiocParallel_1.24.1         scatterpie_0.1.7           
# [7] munsell_0.5.0               codetools_0.2-18            statmod_1.4.35             
# [10] withr_2.4.3                 colorspace_2.0-0            GOSemSim_2.16.1            
# [13] Biobase_2.50.0              rstudioapi_0.13             stats4_4.0.5               
# [16] ggsignif_0.6.0              DOSE_3.16.0                 labeling_0.4.2             
# [19] Rdpack_2.1.1                MatrixGenerics_1.2.1        emmeans_1.5.4              
# [22] GenomeInfoDbData_1.2.4      polyclip_1.10-0             bit64_4.0.5                
# [25] farver_2.0.3                pheatmap_1.0.12             rprojroot_2.0.2            
# [28] downloader_0.4              coda_0.19-4                 vctrs_0.4.1                
# [31] generics_0.1.2              TH.data_1.0-10              doParallel_1.0.16          
# [34] R6_2.5.0                    GenomeInfoDb_1.26.7         clue_0.3-59                
# [37] graphlayouts_0.7.2          locfit_1.5-9.4              bitops_1.0-6               
# [40] cachem_1.0.4                fgsea_1.16.0                DelayedArray_0.16.3        
# [43] assertthat_0.2.1            scales_1.1.1                multcomp_1.4-16            
# [46] ggraph_2.0.5                enrichplot_1.10.2           gtable_0.3.0               
# [49] Cairo_1.5-12.2              egg_0.4.5                   tidygraph_1.2.0            
# [52] sandwich_3.0-0              rlang_1.0.2                 slanter_0.2-0              
# [55] GlobalOptions_0.1.2         splines_4.0.5               rstatix_0.7.0              
# [58] broom_0.7.5                 BiocManager_1.30.10         reshape2_1.4.4             
# [61] abind_1.4-5                 modelr_0.1.8                backports_1.2.1            
# [64] qvalue_2.22.0               clusterProfiler_3.18.1      tools_4.0.5                
# [67] ellipsis_0.3.2              gplots_3.1.1                RColorBrewer_1.1-2         
# [70] BiocGenerics_0.36.1         Rcpp_1.0.6                  plyr_1.8.6                 
# [73] progress_1.2.2              zlibbioc_1.36.0             RCurl_1.98-1.3             
# [76] prettyunits_1.1.1           ggpubr_0.4.0                GetoptLong_1.0.5           
# [79] viridis_0.5.1               cowplot_1.1.1               S4Vectors_0.28.1           
# [82] zoo_1.8-8                   SummarizedExperiment_1.20.0 haven_2.3.1                
# [85] ggrepel_0.9.1               cluster_2.1.2               fs_1.5.0                   
# [88] variancePartition_1.25.6    magrittr_2.0.1              data.table_1.14.0          
# [91] DO.db_2.9                   openxlsx_4.2.3              circlize_0.4.12            
# [94] reprex_1.0.0                mvtnorm_1.1-1               packrat_0.7.0              
# [97] matrixStats_0.58.0          hms_1.0.0                   GSVA_1.38.2                
# [100] xtable_1.8-4                pbkrtest_0.5-0.1            RhpcBLASctl_0.21-247.1     
# [103] XML_3.99-0.6                rio_0.5.16                  readxl_1.3.1               
# [106] IRanges_2.24.1              gridExtra_2.3               shape_1.4.6                
# [109] compiler_4.0.5              KernSmooth_2.23-18          crayon_1.4.1               
# [112] shadowtext_0.0.9            minqa_1.2.4                 ggfun_0.0.4                
# [115] lubridate_1.7.9.2           DBI_1.1.1                   tweenr_1.0.2               
# [118] dbplyr_2.1.0                MASS_7.3-53.1               boot_1.3-27                
# [121] Matrix_1.3-2                car_3.0-10                  cli_3.3.0                  
# [124] rbibutils_2.0               parallel_4.0.5              igraph_1.2.6               
# [127] GenomicRanges_1.42.0        pkgconfig_2.0.3             rvcheck_0.1.8              
# [130] foreign_0.8-81              foreach_1.5.1               xml2_1.3.2                 
# [133] annotate_1.68.0             XVector_0.30.0              GeneOverlap_1.26.0         
# [136] estimability_1.3            rvest_0.3.6                 digest_0.6.27              
# [139] graph_1.68.0                cellranger_1.1.0            fastmatch_1.1-0            
# [142] edgeR_3.32.1                GSEABase_1.52.1             curl_4.3                   
# [145] gtools_3.8.2                nloptr_1.2.2.2              rjson_0.2.20               
# [148] nlme_3.1-152                lifecycle_1.0.0             jsonlite_1.7.2             
# [151] aod_1.3.1                   carData_3.0-4               viridisLite_0.3.0          
# [154] limma_3.46.0                pillar_1.4.7                lattice_0.20-41            
# [157] fastmap_1.1.0               httr_1.4.2                  survival_3.2-10            
# [160] GO.db_3.12.1                glue_1.6.2                  zip_2.1.1                  
# [163] iterators_1.0.13            png_0.1-7                   bit_4.0.4                  
# [166] ggforce_0.3.3               stringi_1.5.3               blob_1.2.1                 
# [169] org.Hs.eg.db_3.12.0         caTools_1.18.1              memoise_2.0.0  







