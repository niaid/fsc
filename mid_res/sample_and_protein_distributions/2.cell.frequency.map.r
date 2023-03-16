# this script is run in R 4.0.5
# sample barplot
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(here))

figpath = here("mid_res/sample_and_protein_distributions/figures/")


s = readRDS(file = "data/h1h5_annotated_with_meta.rds")

freq_plot = s@meta.data %>%
  select(celltype_joint, sample) %>%
  group_by(celltype_joint, sample) %>%
  summarise(n = n()) %>%
  group_by(sample) %>%
  mutate(log_cell = log10(n)) %>%
  select(sample, celltype_joint, log_cell) %>%
  spread(celltype_joint, log_cell) %>%
  mutate(sample = if_else(str_sub(sample, 1, 2) == 'H5',
                          true = str_sub(sample, -6, -1),
                          false = sample)) %>%
  column_to_rownames("sample") %>%
  t()

annotation = read_delim(file = "data/full_metadata/full_sample_metadata.txt", delim = "\t")

md = s@meta.data %>% 
  select(sample) %>%
  group_by(sample) %>%
  summarise(n_cells = log10(n())) %>%
  mutate(subjectID = str_sub(sample, -6, -4)) %>%
  mutate(timepoint = str_sub(sample, -2, -1)) %>%
  mutate(group = plyr::mapvalues(
    x = subjectID,
    from = annotation$subjectid,
    to = annotation$adjMFC)) %>%
  select(sample, group, timepoint) %>%
  mutate(group = if_else(str_sub(sample, 1, 2) == 'H5', true = "adjuvant", false = group)) %>%
  mutate(sample = if_else(
    str_sub(sample, 1, 2) == 'H5',
    true = str_sub(sample, -6, -1),
    false = sample)) %>%
  column_to_rownames("sample")
  
# quant palette 
mat_colors = list(
  group = c("grey", "red", "deepskyblue3"),
  timepoint = c("grey", "orange", "black")
)
names(mat_colors$timepoint) = unique(md$timepoint)
names(mat_colors$group) = unique(md$group)

# plot 
rownames(freq_plot) = str_replace_all(rownames(freq_plot),pattern = '_', replacement = ' ')
pheatmap::pheatmap(freq_plot,
                   annotation_col = md,
                   display_numbers = FALSE, 
                   color = grDevices::colorRampPalette(colors = c("gray100", "black" ))(100),
                   cluster_cols = F, border_color = NA,
                   width = 10, height = 5,treeheight_row = 20,
                   annotation_colors = mat_colors,
                   filename = paste0(figpath, 'sample_celltype_map.pdf')
)
 

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] here_1.0.1         forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0       
# [7] tidyr_1.1.2        tibble_3.1.8       ggplot2_3.3.3      tidyverse_1.3.0    sp_1.4-5           SeuratObject_4.1.0
# [13] Seurat_4.0.1      
# 
# loaded via a namespace (and not attached):
# [1] Rtsne_0.15            colorspace_2.0-0      deldir_1.0-6          ellipsis_0.3.2        ggridges_0.5.3       
# [6] rprojroot_2.0.2       fs_1.5.0              rstudioapi_0.13       spatstat.data_2.1-0   farver_2.0.3         
# [11] leiden_0.3.7          listenv_0.8.0         ggrepel_0.9.1         lubridate_1.8.0       fansi_0.4.2          
# [16] xml2_1.3.2            codetools_0.2-18      splines_4.0.5         polyclip_1.10-0       jsonlite_1.7.2       
# [21] packrat_0.7.0         broom_0.7.5           ica_1.0-2             cluster_2.1.2         dbplyr_2.1.0         
# [26] png_0.1-7             rgeos_0.5-9           pheatmap_1.0.12       uwot_0.1.10           shiny_1.6.0          
# [31] sctransform_0.3.2     spatstat.sparse_2.0-0 compiler_4.0.5        httr_1.4.2            backports_1.2.1      
# [36] assertthat_0.2.1      Matrix_1.4-1          fastmap_1.1.0         lazyeval_0.2.2        cli_3.4.1            
# [41] later_1.1.0.1         htmltools_0.5.2       tools_4.0.5           igraph_1.2.6          gtable_0.3.0         
# [46] glue_1.6.2            RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.9            scattermore_0.7      
# [51] cellranger_1.1.0      vctrs_0.5.1           nlme_3.1-152          progressr_0.10.0      lmtest_0.9-38        
# [56] globals_0.14.0        rvest_0.3.6           mime_0.10             miniUI_0.1.1.1        lifecycle_1.0.3      
# [61] irlba_2.3.3           goftest_1.2-2         future_1.21.0         MASS_7.3-53.1         zoo_1.8-8            
# [66] scales_1.1.1          spatstat.core_2.0-0   hms_1.0.0             promises_1.2.0.1      spatstat.utils_2.3-0 
# [71] parallel_4.0.5        RColorBrewer_1.1-2    reticulate_1.18       pbapply_1.4-3         gridExtra_2.3        
# [76] rpart_4.1-15          stringi_1.5.3         rlang_1.0.6           pkgconfig_2.0.3       matrixStats_0.58.0   
# [81] lattice_0.20-41       ROCR_1.0-11           tensor_1.5            patchwork_1.1.1       htmlwidgets_1.5.3    
# [86] cowplot_1.1.1         tidyselect_1.2.0      parallelly_1.23.0     RcppAnnoy_0.0.18      plyr_1.8.6           
# [91] magrittr_2.0.3        R6_2.5.0              generics_0.1.2        DBI_1.1.1             withr_2.4.3          
# [96] pillar_1.8.1          haven_2.4.3           mgcv_1.8-34           fitdistrplus_1.1-3    survival_3.2-10      
# [101] abind_1.4-5           future.apply_1.7.0    crayon_1.4.1          modelr_0.1.8          KernSmooth_2.23-18   
# [106] utf8_1.2.2            spatstat.geom_2.4-0   plotly_4.9.3          readxl_1.3.1          grid_4.0.5           
# [111] data.table_1.14.0     reprex_1.0.0          digest_0.6.27         xtable_1.8-4          httpuv_1.5.5         
# [116] munsell_0.5.0         viridisLite_0.3.0    

