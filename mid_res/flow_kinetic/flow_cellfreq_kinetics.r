# flow activated monocyte kinetic 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(emmeans))
source(here('functions/MattPMutils.r'))

# save path 
figpath = here('mid_res/flow_kinetic/figures/'); 
dir.create(figpath,recursive = TRUE)

# test only innate subsets based on hypothesis generated from CITE-seq data
# load flow data 
fp = data.table::fread('data/CHI_H1N1_data/flow/flow_annotation.txt', header = TRUE)
id.select = paste('ID', 64:78,sep = '') %>%  as.character()
fp = fp %>% filter(`Population ID` %in% id.select)

# load flow data day 1 fold changes 
fd = 
  data.table::fread(file = here('data/CHI_H1N1_data/flow/day1-day0.log10.txt'),header = TRUE) %>% 
  filter(ID %in% fp$`Population ID`) %>%  
  mutate(subset_name = plyr::mapvalues(ID, from = fp$`Population ID`,to = fp$`Subset name`)) %>%
  select(-ID) %>% 
  column_to_rownames('subset_name') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Subject')


# map titers to data 
titer = data.table::fread(file = here('data/CHI_H1N1_data/titer/titer_processed.txt'))
fd$adjMFC_class = plyr::mapvalues(x = fd$Subject, from = titer$Subject,to = titer$adjMFC_class )
fd = fd %>% select(adjMFC_class, everything())
fd$adjMFC_class = factor(fd$adjMFC_class, levels = c('0','1','2'))
fd = fd[!is.na(fd$adjMFC_class), ]

# baseline 
# these are raw percentages so use non parametric rank stats 
fd3 = 
  data.table::fread(file = here('data/CHI_H1N1_data/flow/day0.raw.txt'),header = TRUE) %>% 
  filter(ID %in% fp$`Population ID`) %>%  
  mutate(subset_name = plyr::mapvalues(ID, from = fp$`Population ID`,to = fp$`Subset name`)) %>%
  select(-ID) %>% 
  column_to_rownames('subset_name') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Subject')
# map adj mfc class 
fd3$adjMFC_class = plyr::mapvalues(x = fd3$Subject, from = titer$Subject,to = titer$adjMFC_class )
fd3 = fd3 %>% select(adjMFC_class, everything())
fd3$adjMFC_class = factor(fd3$adjMFC_class, levels = c('0','1','2'))
fd3 = fd3[!is.na(fd3$adjMFC_class), ]
fd3 = fd3[!fd3$adjMFC_class == '1', ]

wilcox.res3 = apply(
  X =  fd3[, 3:ncol(fd3)],
  MARGIN = 2,
  FUN = function(x) {
    wilcox.test(x ~ fd3$adjMFC_class) %>%  broom::tidy()
  }) %>% 
  bind_rows(.id = 'subset')
wilcox.res3 %>%  filter(p.value < 0.1)
saveRDS(wilcox.res3,file = paste0(datapath,'wilcox.res3.rds'))


# comparison 
flow_compare = list(c('2','0'))

# color specification
cu1 = sapply(c('dodgerblue', 'red'), col.alpha, 0.2) %>% unname()
cu2 = c('dodgerblue', 'red')

# theme 
mtheme = list(
  theme_bw(), 
  geom_boxplot(show.legend = FALSE, outlier.shape = 21), 
  ggpubr::stat_compare_means(comparisons = flow_compare,method = 'wilcox', paired = FALSE),
  scale_y_continuous( breaks= scales::pretty_breaks(), expand = c(0.15,0)), 
  theme(axis.text.y = element_text(size = 8), 
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = 'black')),
  scale_x_discrete(labels = c("low", "high")),
  xlab("Antibody \n Response"),
  scale_fill_manual(values = cu1),
  scale_color_manual(values = cu2) 
)

# plot 
p3 = ggplot(fd3, aes(x = adjMFC_class, y = `activated monocyte HLA-DR+`,
                     color = adjMFC_class,
                     fill = adjMFC_class)) +  
  mtheme
ggsave(p3,filename = paste0(figpath, 'mono_HLADR.pdf'), width = 2, height = 3)



# longitudinal analysis of cell population frequency
# connection between baseline (-7 and 0) and day 1 kinetics 
d01 = data.table::fread(file = here('data/CHI_H1N1_data/flow/longitudinal/pre7.raw.txt'), header = TRUE)
d02 = data.table::fread(file = here('data/CHI_H1N1_data/flow/longitudinal/day0.raw.txt'), header = TRUE)
d1 = data.table::fread(file = here('data/CHI_H1N1_data/flow/longitudinal/day1.raw.txt'), header = TRUE)
d7 = data.table::fread(file = here('data/CHI_H1N1_data/flow/longitudinal/day7.raw.txt'), header = TRUE)
d70 = data.table::fread(file = here('data/CHI_H1N1_data/flow/longitudinal/day70.raw.txt'), header = TRUE)

d.list = list('t01' = d01, 't02' = d02, 't1' = d1, 't7' = d7, 't70' = d70 )

# format, combine and label by time 
d = lapply(d.list, function(x) 
  x %>%  
    filter(ID %in% fp$`Population ID`) %>%  
    mutate(subset_name = plyr::mapvalues(ID, from = fp$`Population ID`,to = fp$`Subset name`)) %>%
    select(-ID) %>% 
    column_to_rownames('subset_name') %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column('Subject') ) %>% 
  bind_rows(.id = 'timepoint')
# map adj mfc class 
d$adjMFC_class = plyr::mapvalues(x = d$Subject, from = titer$Subject,to = titer$adjMFC_class )
d = d %>%  filter(adjMFC_class %in% c('0','2'))
d$adjMFC_class = factor(d$adjMFC_class, levels = c('0','2'))
d$timepoint = factor(d$timepoint, levels = c("t01", "t02", "t1",  "t7",  "t70"))

# plot time course 
mtheme2 = list(
  theme_bw(), 
  theme(axis.text.y = element_text(size = 8), 
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = 'black')),
  scale_fill_manual(values = cu1),
  scale_color_manual(values = cu2) 
)

# monocyte
d2 = d[!is.na(d$`activated monocyte HLA-DR+`), ]
p3=
  ggplot(d2, aes(x = timepoint, y = `activated monocyte HLA-DR+`, color = adjMFC_class, group = Subject)) + 
  geom_line(size = 0.5, alpha = 0.2, show.legend = FALSE) + 
  geom_point(size = 0.5, alpha = 0.2, shape = 21, show.legend = FALSE) + 
  geom_smooth(data = d2, size  = 2.5, 
              method = 'loess',
              aes(x = timepoint,  y = `activated monocyte HLA-DR+`, 
                  color = adjMFC_class,
                  fill = adjMFC_class,
                  group = adjMFC_class), alpha = 0.2,
              se = TRUE,  show.legend = FALSE) +
  scale_x_discrete(expand = c(0,0.1)) +
  mtheme2 
p3
ggsave(p3,filename = paste0(figpath, 'monoDRfreqtime.pdf'), width = 3, height = 3)

### mixed model 
d2 = d %>% 
  select(timepoint, subjectid = Subject,
         drmono = `activated monocyte HLA-DR+`,
         adjmfc = adjMFC_class) %>% 
  mutate(timepoint= as.character(timepoint)) %>% 
  mutate(timepoint = ifelse(timepoint %in% c('t01','t02'),yes =  't0',no =  timepoint)) %>% 
  mutate(timepoint = factor(timepoint , levels = c( "t0",  "t1" , "t7"  ,"t70")))

# fit model 
m = lme4::lmer(drmono ~ timepoint + (1|subjectid),data = d2)
m2 = lme4::lmer(drmono ~ timepoint*adjmfc + (1|subjectid),data = d2)
anova(m,m2)
# Data: d2
# Models:
#   m: drmono ~ timepoint + (1 | subjectid)
# m2: drmono ~ timepoint * adjmfc + (1 | subjectid)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# m     6 615.95 632.58 -301.98   603.95                       
# m2   10 615.02 642.73 -297.51   595.02 8.9311  4    0.06284 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
emm1 = emmeans(m2, revpairwise ~ timepoint|adjmfc)
p = emm1$contrasts[c(1,7), ] %>%  plot() + theme_bw() + 
  xlab('effect size') + 
  ylab('group') + 
  theme(text = element_text(size = 5))
p
ggsave(p,filename = paste0(figpath, 'monoDREMM.pdf'), width = 3, height = 1)


sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] emmeans_1.5.4   scglmmr_0.1.0   forcats_0.5.1   stringr_1.4.0   dplyr_1.0.4     purrr_0.3.4     readr_1.4.0     tidyr_1.1.2    
# [9] tibble_3.1.8    ggplot2_3.3.3   tidyverse_1.3.0 here_1.0.1     
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.2                  tidyselect_1.2.0            lme4_1.1-26                 htmlwidgets_1.5.3          
# [5] RSQLite_2.2.7               AnnotationDbi_1.52.0        grid_4.0.5                  BiocParallel_1.24.1        
# [9] scatterpie_0.1.7            munsell_0.5.0               codetools_0.2-18            statmod_1.4.35             
# [13] withr_2.4.3                 colorspace_2.0-0            GOSemSim_2.16.1             Biobase_2.50.0             
# [17] knitr_1.39                  rstudioapi_0.13             stats4_4.0.5                ggsignif_0.6.0             
# [21] DOSE_3.16.0                 labeling_0.4.2              MatrixGenerics_1.2.1        Rdpack_2.1.1               
# [25] GenomeInfoDbData_1.2.4      polyclip_1.10-0             bit64_4.0.5                 farver_2.0.3               
# [29] pheatmap_1.0.12             rprojroot_2.0.2             downloader_0.4              coda_0.19-4                
# [33] vctrs_0.5.1                 generics_0.1.2              TH.data_1.0-10              xfun_0.30                  
# [37] R6_2.5.0                    doParallel_1.0.16           GenomeInfoDb_1.26.7         graphlayouts_0.7.2         
# [41] locfit_1.5-9.4              bitops_1.0-6                cachem_1.0.4                fgsea_1.16.0               
# [45] DelayedArray_0.16.3         assertthat_0.2.1            scales_1.1.1                nnet_7.3-15                
# [49] multcomp_1.4-16             ggraph_2.0.5                enrichplot_1.10.2           gtable_0.3.0               
# [53] egg_0.4.5                   tidygraph_1.2.0             sandwich_3.0-0              rlang_1.0.6                
# [57] slanter_0.2-0               splines_4.0.5               rstatix_0.7.0               checkmate_2.0.0            
# [61] broom_0.7.5                 BiocManager_1.30.10         reshape2_1.4.4              abind_1.4-5                
# [65] modelr_0.1.8                backports_1.2.1             Hmisc_4.5-0                 qvalue_2.22.0              
# [69] clusterProfiler_3.18.1      tools_4.0.5                 ellipsis_0.3.2              gplots_3.1.1               
# [73] RColorBrewer_1.1-2          BiocGenerics_0.36.1         Rcpp_1.0.9                  plyr_1.8.6                 
# [77] base64enc_0.1-3             progress_1.2.2              zlibbioc_1.36.0             RCurl_1.98-1.3             
# [81] prettyunits_1.1.1           rpart_4.1-15                ggpubr_0.4.0                viridis_0.5.1              
# [85] cowplot_1.1.1               S4Vectors_0.28.1            zoo_1.8-8                   cluster_2.1.2              
# [89] SummarizedExperiment_1.20.0 haven_2.4.3                 ggrepel_0.9.1               fs_1.5.0                   
# [93] variancePartition_1.25.6    magrittr_2.0.3              data.table_1.14.0           DO.db_2.9                  
# [97] openxlsx_4.2.3              reprex_1.0.0                mvtnorm_1.1-1               packrat_0.7.0              
# [101] matrixStats_0.58.0          hms_1.0.0                   GSVA_1.38.2                 xtable_1.8-4               
# [105] pbkrtest_0.5-0.1            RhpcBLASctl_0.21-247.1      XML_3.99-0.6                jpeg_0.1-8.1               
# [109] rio_0.5.16                  readxl_1.3.1                IRanges_2.24.1              gridExtra_2.3              
# [113] compiler_4.0.5              KernSmooth_2.23-18          crayon_1.4.1                shadowtext_0.0.9           
# [117] htmltools_0.5.2             minqa_1.2.4                 mgcv_1.8-34                 ggfun_0.0.4                
# [121] Formula_1.2-4               lubridate_1.8.0             DBI_1.1.1                   corrplot_0.84              
# [125] tweenr_1.0.2                dbplyr_2.1.0                MASS_7.3-53.1               boot_1.3-27                
# [129] Matrix_1.4-1                car_3.0-10                  cli_3.4.1                   rbibutils_2.0              
# [133] parallel_4.0.5              igraph_1.2.6                GenomicRanges_1.42.0        pkgconfig_2.0.3            
# [137] rvcheck_0.1.8               foreign_0.8-81              xml2_1.3.2                  foreach_1.5.1              
# [141] annotate_1.68.0             XVector_0.30.0              GeneOverlap_1.26.0          estimability_1.3           
# [145] rvest_0.3.6                 digest_0.6.27               graph_1.68.0                cellranger_1.1.0           
# [149] fastmatch_1.1-0             htmlTable_2.1.0             edgeR_3.32.1                GSEABase_1.52.1            
# [153] curl_4.3                    gtools_3.8.2                nloptr_1.2.2.2              lifecycle_1.0.3            
# [157] nlme_3.1-152                jsonlite_1.7.2              aod_1.3.1                   carData_3.0-4              
# [161] viridisLite_0.3.0           limma_3.46.0                fansi_0.4.2                 pillar_1.8.1               
# [165] lattice_0.20-41             fastmap_1.1.0               httr_1.4.2                  survival_3.2-10            
# [169] GO.db_3.12.1                glue_1.6.2                  UpSetR_1.4.0                zip_2.1.1                  
# [173] png_0.1-7                   iterators_1.0.13            bit_4.0.4                   ggforce_0.3.3              
# [177] stringi_1.5.3               blob_1.2.1                  org.Hs.eg.db_3.12.0         latticeExtra_0.6-29        
# [181] caTools_1.18.1              memoise_2.0.0   