suppressMessages(library(ComplexHeatmap))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(here))
# functions
source("functions/analysis_functions.R")
source("functions/protein_annotation_functions.r")

# save path 
figpath = here("mid_res/aggregated_protein_libraries/figures/")
dir.create(figpath)

# cluster information combined heatmap
h1 = readRDS(file = here("data/h1h5_annotated_with_meta.rds")) %>%
  SetAllIdent(id = "celltype_joint") %>% 
  SubsetData(ident.remove = "DOUBLET")

# specify subsets of proteins for fig 1 
prot_order = c(
  # B 
  "CD19", "CD20", "IgM", "IgD", "CD40", "CD185", 
  # pdc 
  "CD123", "CD303", 
  # my lin / dc 
  "HLA-DR", "CD11c",  
  # DC and  mono 
  "CD71", "CD14", "CD16",  "CD11b", "CD1d", "CD1c", 
  # hsc
  "CD34",
  # T 
  "CD3", "CD45RO", "CD45RA" , "CD62L", "CD27", 
  "CD28", "CD279", "CD4", "CD8", "CD161",
  # nk
  "CD244", "KLRG1", "CD127", "CD56", "CD57", "CD38",
  # state 
  "CD103" ,  "CD196", "CD195", "CD25", "CD86", 
  "CD69",     "CD31"
)


# aggregate protein data 
# single cell data 
prot.dat = as.data.frame(t(h1@assay$CITE@data))
# replact ADT string 
colnames(prot.dat) = str_sub(colnames(prot.dat), start = 1, end = -6)

# aggregate (mean)
prot_data = cbind(prot.dat, h1@meta.data) %>%
  group_by(sample, subject_id = sampleid , timepoint,  
           time_cohort, batch, age, gender, celltype =  celltype_joint,
           antibody_response =  adjmfc.group)  %>%
  summarize_at(.vars = colnames(prot.dat), .funs = median) %>%
  ungroup() %>%
  select(sample, celltype, prot_order) %>%
  arrange(celltype, sample) %>%
  unite(col = "sample_celltype", sample:celltype, sep = "_") %>%
  column_to_rownames("sample_celltype") %>% 
  t()


# cell frequency 
md = h1@meta.data
df = md %>% 
  group_by(sample, subject_id = sampleid, timepoint, 
           time_cohort, batch, 
           age, gender, celltype = celltype_joint, antibody_response =  adjmfc.group) %>% 
  summarize(count = n(), log_lib_size = log10(sum(nUMI))) %>% 
  group_by(sample) %>% 
  mutate(cell_freq=count/sum(count)*100) %>% 
  arrange(celltype, sample) %>% 
  mutate(log_cell_count = log10(count))
cellfreq = df$cell_freq

# celltype 
cellt = df$celltype
ha = HeatmapAnnotation(celltype = cellt, 
                       col = list(celltype = c(
                         "BC_Mem" = "lightslateblue",
                         "BC_Naive" = "#2B3D26",       
                         "CD103_Tcell" = "#E25822",       
                         "CD14_Mono"= "red",       
                         "CD16_Mono"  = "firebrick4",       
                         "CD38_Bcell" = "#882D17",       
                         "CD4_CD161_Mem_Tcell" = "navy",       
                         "CD4_CD25_Tcell"= "#B3446C",       
                         "CD4_CD56_Tcell" = "maroon1",       
                         "CD4_CD57_Tcell" = "#604E97",       
                         "CD4_Efct_Mem_Tcell" ="#F99379",       
                         "CD4Naive_Tcell" = "#0067A5",       
                         "CD8_CD161_Tcell" = "olivedrab", 
                         "CD8_Mem_Tcell" = "#008856",       
                         "CD8_Naive_Tcell" = "#848482",       
                         "CD8_NKT" = "#C2B280",       
                         "HSC" = "#BE0032",       
                         "IgA_CD14_Mono" = "#A1CAF1",       
                         "MAIT_Like" = "#F38400",       
                         "mDC" = "#875692",       
                         "NK" = "#F3C300",     
                         "pDC" = "#222222"))
)

# Create annotations 
libraries_map = columnAnnotation(
  log_library_size = column_anno_points(
    df$log_lib_size, 
    size = unit(0.3, 'mm'),
    pch = 21, axis = TRUE, border = TRUE,
    gp = gpar(color = "black")
  ),
  height = unit(1.8, units = "cm")
)


# matrix color values 
col_fun = circlize::colorRamp2(breaks = c(-1,0,2,4,8,12,16,20),
                               colors = viridis::viridis(n = 8, option = "B"))

# organize proteins by lineages
rownames(prot_data) = factor(rownames(prot_data),levels = prot_order)

# cluster by column; save heatmap 
pdf(paste0(figpath,"heatmap3.pdf"), width = 7, height = 6.5)
draw(
  ComplexHeatmap::Heatmap(
    matrix = prot_data, 
    name = "", 
    col = col_fun, 
    row_names_gp = gpar(color = "black", fontsize = 10),
    top_annotation = ha,
    bottom_annotation = libraries_map,
    show_column_names = FALSE, 
    cluster_rows = FALSE, 
    cluster_columns = TRUE,
    use_raster = TRUE), show_annotation_legend = FALSE)
dev.off()

sessionInfo()
# R version 3.5.3 Patched (2019-03-11 r77192)
# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ComplexHeatmap_1.20.0 viridis_0.5.1         viridisLite_0.3.0     here_0.1              forcats_0.4.0         stringr_1.4.0         dplyr_0.8.5           purrr_0.3.3          
# [9] readr_1.3.1           tidyr_1.0.2           tibble_2.1.1          tidyverse_1.2.1       Seurat_2.3.4          Matrix_1.2-15         cowplot_0.9.4         ggplot2_3.1.1        
# 
# loaded via a namespace (and not attached):
# [1] circlize_0.4.10         readxl_1.3.1            snow_0.4-3              backports_1.1.4         Hmisc_4.2-0             plyr_1.8.4              igraph_1.2.4.1         
# [8] lazyeval_0.2.2          splines_3.5.3           crosstalk_1.0.0         digest_0.6.25           foreach_1.4.4           htmltools_0.3.6         lars_1.2               
# [15] fansi_0.4.0             gdata_2.18.0            magrittr_2.0.1          checkmate_1.9.3         cluster_2.0.7-1         mixtools_1.1.0          ROCR_1.0-7             
# [22] modelr_0.1.4            R.utils_2.8.0           colorspace_1.4-1        rvest_0.3.4             haven_2.1.0             xfun_0.7                crayon_1.3.4           
# [29] jsonlite_1.6            survival_2.43-3         zoo_1.8-6               iterators_1.0.10        ape_5.3                 glue_1.3.1              pals_1.5               
# [36] gtable_0.3.0            webshot_0.5.1           GetoptLong_1.0.2        kernlab_0.9-27          shape_1.4.4             prabclus_2.3-1          DEoptimR_1.0-8         
# [43] maps_3.3.0              scales_1.0.0            pheatmap_1.0.12         mvtnorm_1.0-10          bibtex_0.4.2            miniUI_0.1.1.1          Rcpp_1.0.1             
# [50] metap_1.1               dtw_1.20-1              xtable_1.8-4            htmlTable_1.13.1        reticulate_1.12         foreign_0.8-71          bit_1.1-14             
# [57] mapproj_1.2.6           proxy_0.4-23            mclust_5.4.5            SDMTools_1.1-221.1      Formula_1.2-3           stats4_3.5.3            tsne_0.1-3             
# [64] htmlwidgets_1.3         httr_1.4.0              gplots_3.0.1.1          RColorBrewer_1.1-2      fpc_2.2-1               ellipsis_0.3.0          acepack_1.4.1          
# [71] modeltools_0.2-22       ica_1.0-2               pkgconfig_2.0.2         R.methodsS3_1.7.1       flexmix_2.3-15          nnet_7.3-12             utf8_1.1.4             
# [78] manipulateWidget_0.10.0 tidyselect_0.2.5        labeling_0.3            rlang_0.4.5             reshape2_1.4.3          later_0.8.0             munsell_0.5.0          
# [85] cellranger_1.1.0        tools_3.5.3             cli_1.1.0               generics_0.0.2          broom_0.5.2             ggridges_0.5.1          npsurv_0.4-0           
# [92] knitr_1.23              bit64_0.9-7             fitdistrplus_1.0-14     robustbase_0.93-5       rgl_0.100.30            caTools_1.17.1.2        RANN_2.6.1             
# [99] packrat_0.5.0           pbapply_1.4-0           nlme_3.1-137            mime_0.6                R.oo_1.22.0             xml2_1.2.0              hdf5r_1.2.0            
# [106] compiler_3.5.3          rstudioapi_0.10         png_0.1-7               lsei_1.2-0              stringi_1.4.3           lattice_0.20-38         ggsci_2.9              
# [113] vctrs_0.2.4             pillar_1.4.1            lifecycle_0.1.0         GlobalOptions_0.1.2     Rdpack_0.11-0           lmtest_0.9-37           data.table_1.12.2      
# [120] bitops_1.0-6            irlba_2.3.3             gbRd_0.4-11             httpuv_1.5.1            R6_2.4.0                latticeExtra_0.6-28     promises_1.0.1         
# [127] KernSmooth_2.23-15      gridExtra_2.3           sessioninfo_1.1.1       codetools_0.2-16        dichromat_2.0-0         MASS_7.3-51.1           gtools_3.8.1           
# [134] assertthat_0.2.1        rjson_0.2.20            rprojroot_1.3-2         withr_2.1.2             diptest_0.75-7          parallel_3.5.3          doSNOW_1.0.16          
# [141] hms_0.4.2               rpart_4.1-13            class_7.3-15            segmented_0.5-4.0       Rtsne_0.15              shiny_1.3.2             lubridate_1.7.4        
# [148] base64enc_0.1-3        
