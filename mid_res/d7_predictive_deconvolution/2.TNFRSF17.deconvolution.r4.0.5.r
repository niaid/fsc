# R 4.0.5
suppressMessages(library(tidyverse))
suppressMessages(library(scglmmr))
suppressMessages(library(here))
suppressMessages(library(magrittr))
source('functions/analysis_functions.R')

# save path 
figpath = here('mid_res/d7_predictive_deconvolution/figures/')

# single cell composition of m156 
core7 = readRDS(file = here('signature_curation/core_d7.rds'))
h1 = readRDS(file = here('data/h1h5_annotated_with_meta.rds'))

m156_sub = intersect(rownames(h1@data), core7$`LI.M156.0 plasma cell b cell Ig`)
gene = as.data.frame(as.matrix(t(as.matrix(h1@data[m156_sub, ]))))
#pdf = cbind(gene, h1@meta.data)
#pdf = pdf %>% gather(gene, normcount, m156_sub[1]:m156_sub[length(m156_sub)])

# composisiton of TNFRSF17
gene = as.data.frame(as.matrix(t(as.matrix(h1@raw.data[m156_sub, ]))))
pdf = cbind(gene, h1@meta.data)
pdf = pdf %>% gather(gene, count, m156_sub[1]:m156_sub[length(m156_sub)])
tnf = pdf %>% 
  filter(cohort == 'H1N1') %>% 
  filter(gene == "TNFRSF17" & timepoint %in% c("d0", "d7") ) %>% 
  group_by(celltype_joint, timepoint, gene) %>% 
  summarise(n = sum(count)) %>%
  ungroup() 


# visualize distribution
tnf$celltype_joint = str_replace_all(tnf$celltype_joint, pattern = '_',replacement = ' ')
p = ggplot(tnf, aes(x = reorder(celltype_joint, n), y = n, fill = timepoint)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  theme_bw() +
  coord_flip() + 
  theme(axis.text.x =  element_text(size = 9, color = 'black'), axis.title.x = element_text(color = 'black')) + 
  theme(axis.text.y =  element_text(size = 9, color = 'black'), axis.title.y = element_text(color = 'black')) + 
  ylab("UMI counts") + 
  scale_fill_manual(values = c("grey48", "red")) + 
  theme(axis.title.y  =  element_blank()) + 
  theme(legend.position = c(0.7, 0.2) ) + 
  ggtitle("TNFRSF17 Gene") 
p
ggsave(p, filename = paste0(figpath, "TNFRSF17_composition.pdf"),width = 3, height = 4.3)

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] magrittr_2.0.3    here_1.0.1        scglmmr_0.1.0     forcats_0.5.1     stringr_1.4.0     dplyr_1.0.4       purrr_0.3.4      
# [8] readr_1.4.0       tidyr_1.1.2       tibble_3.1.8      ggplot2_3.3.3     tidyverse_1.3.0   viridis_0.5.1     viridisLite_0.3.0
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.2                  tidyselect_1.2.0            lme4_1.1-26                 RSQLite_2.2.7              
# [5] AnnotationDbi_1.52.0        grid_4.0.5                  BiocParallel_1.24.1         scatterpie_0.1.7           
# [9] munsell_0.5.0               codetools_0.2-18            statmod_1.4.35              withr_2.4.3                
# [13] colorspace_2.0-0            GOSemSim_2.16.1             Biobase_2.50.0              rstudioapi_0.13            
# [17] ROCR_1.0-11                 stats4_4.0.5                ggsignif_0.6.0              DOSE_3.16.0                
# [21] labeling_0.4.2              MatrixGenerics_1.2.1        Rdpack_2.1.1                emmeans_1.5.4              
# [25] GenomeInfoDbData_1.2.4      polyclip_1.10-0             bit64_4.0.5                 farver_2.0.3               
# [29] pheatmap_1.0.12             rprojroot_2.0.2             downloader_0.4              coda_0.19-4                
# [33] vctrs_0.5.1                 generics_0.1.2              TH.data_1.0-10              R6_2.5.0                   
# [37] doParallel_1.0.16           GenomeInfoDb_1.26.7         graphlayouts_0.7.2          locfit_1.5-9.4             
# [41] bitops_1.0-6                cachem_1.0.4                fgsea_1.16.0                DelayedArray_0.16.3        
# [45] assertthat_0.2.1            scales_1.1.1                multcomp_1.4-16             ggraph_2.0.5               
# [49] enrichplot_1.10.2           gtable_0.3.0                egg_0.4.5                   tidygraph_1.2.0            
# [53] sandwich_3.0-0              rlang_1.0.6                 slanter_0.2-0               splines_4.0.5              
# [57] rstatix_0.7.0               broom_0.7.5                 BiocManager_1.30.10         reshape2_1.4.4             
# [61] abind_1.4-5                 modelr_0.1.8                backports_1.2.1             qvalue_2.22.0              
# [65] clusterProfiler_3.18.1      tools_4.0.5                 ellipsis_0.3.2              gplots_3.1.1               
# [69] RColorBrewer_1.1-2          BiocGenerics_0.36.1         Rcpp_1.0.9                  plyr_1.8.6                 
# [73] progress_1.2.2              zlibbioc_1.36.0             RCurl_1.98-1.3              prettyunits_1.1.1          
# [77] ggpubr_0.4.0                cowplot_1.1.1               S4Vectors_0.28.1            zoo_1.8-8                  
# [81] SummarizedExperiment_1.20.0 haven_2.4.3                 ggrepel_0.9.1               fs_1.5.0                   
# [85] variancePartition_1.25.6    data.table_1.14.0           DO.db_2.9                   openxlsx_4.2.3             
# [89] RANN_2.6.1                  reprex_1.0.0                mvtnorm_1.1-1               packrat_0.7.0              
# [93] matrixStats_0.58.0          hms_1.0.0                   GSVA_1.38.2                 xtable_1.8-4               
# [97] pbkrtest_0.5-0.1            RhpcBLASctl_0.21-247.1      XML_3.99-0.6                rio_0.5.16                 
# [101] readxl_1.3.1                IRanges_2.24.1              gridExtra_2.3               compiler_4.0.5             
# [105] KernSmooth_2.23-18          crayon_1.4.1                shadowtext_0.0.9            minqa_1.2.4                
# [109] ggfun_0.0.4                 lubridate_1.8.0             DBI_1.1.1                   tweenr_1.0.2               
# [113] dbplyr_2.1.0                MASS_7.3-53.1               boot_1.3-27                 Matrix_1.4-1               
# [117] car_3.0-10                  cli_3.4.1                   rbibutils_2.0               parallel_4.0.5             
# [121] igraph_1.2.6                GenomicRanges_1.42.0        pkgconfig_2.0.3             rvcheck_0.1.8              
# [125] foreign_0.8-81              xml2_1.3.2                  foreach_1.5.1               annotate_1.68.0            
# [129] XVector_0.30.0              GeneOverlap_1.26.0          estimability_1.3            rvest_0.3.6                
# [133] digest_0.6.27               graph_1.68.0                cellranger_1.1.0            fastmatch_1.1-0            
# [137] edgeR_3.32.1                GSEABase_1.52.1             curl_4.3                    gtools_3.8.2               
# [141] nloptr_1.2.2.2              lifecycle_1.0.3             nlme_3.1-152                jsonlite_1.7.2             
# [145] aod_1.3.1                   carData_3.0-4               limma_3.46.0                fansi_0.4.2                
# [149] pillar_1.8.1                lattice_0.20-41             fastmap_1.1.0               httr_1.4.2                 
# [153] survival_3.2-10             GO.db_3.12.1                glue_1.6.2                  zip_2.1.1                  
# [157] iterators_1.0.13            bit_4.0.4                   ggforce_0.3.3               stringi_1.5.3              
# [161] blob_1.2.1                  org.Hs.eg.db_3.12.0         caTools_1.18.1              memoise_2.0.0