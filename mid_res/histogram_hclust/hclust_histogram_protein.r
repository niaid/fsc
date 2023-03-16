# high resolution histogram heatmaps
# script uses R 3.5.1
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(ggridges))
suppressMessages(library(ggsci))
suppressMessages(library(viridis))
suppressMessages(library(here))
source(file = "functions/analysis_functions.R")

# save paths 
figpath = here("mid_res/histogram_hclust/figures/")
dir.create(figpath) 


# Define proteins for hclust and visualization 

# T cell markers
tc_markers = c("CD3_PROT", "CD4_PROT", "CD8_PROT", "CD45RA_PROT", "CD45RO_PROT",
               "CD161_PROT", "CD127_PROT", "CD57_PROT", "CD27_PROT", "CD62L_PROT",
               "KLRG1_PROT", "CD103_PROT",  "CD25_PROT", "CD31_PROT")

# B cell markers 
bc_markers = c("CD20_PROT", "CD38_PROT", "IgD_PROT", "CD133_PROT", "IgM_PROT", "CD40_PROT")

# monocyte / dc markers 
mono_markers = c("CD33_PROT", "CD14_PROT", "CD16_PROT", "CD141_PROT", "CD11b_PROT")

# NK markers 
nk_markers = c("CD56_PROT") 

# CD markers 
dc_markers = c("CD1c_PROT")

# rare cell markers 
rare_markers = c("CD303_PROT", "CD123_PROT", "CD34_PROT")

# cell activation markers 
activation = c("CD71_PROT", "CD183_PROT", "CD184_PROT", "CD185_PROT", "CD39_PROT",
               "CD279_PROT", "CD278 _PROT","CD194_PROT", "CD195_PROT", "CD196_PROT",
               "CD117_PROT", "CD244_PROT")

prot_use = c(tc_markers, 
             bc_markers, 
             mono_markers, 
             nk_markers,  
             dc_markers, 
             rare_markers, 
             activation) 
prot_use_plot = str_replace(prot_use, pattern = "_PROT", replacement = "")

c("#000000", "#E69F00", "#56B4E9", "#009E73", 
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# match a color palette to these markers 
my_pal = rev(c(rep("red3", length(tc_markers)), 
               rep("royalblue1", length(bc_markers)), 
               rep("#009E73", length(mono_markers)),
               rep("#0072B2", length(nk_markers)), 
               rep("#D55E00", length(dc_markers)), 
               rep("#CC79A7", length(rare_markers)),
               rep("black", length(activation))
)) 


# read in H1 Seurat object 
h1 = ReadCohort(joint_object_dir = "data/h1h5_annotated_with_meta.rds", cohort = "H1N1")
h1 = h1 %>% SetAllIdent(id = "celltype_joint") %>% SubsetData(ident.remove = "DOUBLET", subset.raw = TRUE)

# get vector of all clusters 
celltypes = h1@meta.data$celltype_joint %>% unique 
h1 = SetAllIdent(h1,id = "celltype_joint") %>%
  SubsetData(max.cells.per.ident = 1000, random.seed = 1, subset.raw = TRUE)


# convert to tidy ; aggregate as the mean of proteins 
adt = h1@assay$CITE@data %>% t %>% as.data.frame() %>% rownames_to_column("cell")
md = h1@meta.data %>% select(celltype = celltype_joint)
adt = cbind(adt, md)
mean_mtx = adt %>% 
  select(celltype, everything()) %>% 
  group_by(celltype) %>% 
  summarize_at(.vars = prot_use, .funs = base::mean) %>% 
  column_to_rownames("celltype") %>% 
  t %>% 
  as.data.frame 

# index for tidying 
index1 = rownames(mean_mtx)[1]
index2 = rownames(mean_mtx)[length(rownames(mean_mtx))]

# order by lineage 
celltype_order = h1@meta.data$celltype_joint %>% unique() %>% sort()
celltype_order = celltype_order[c(12,11,7,8,9,10,13,19,14:16,3,17,21,1,2,6,4,5,18,20,22)]

# alt (not used)
# use hclust within pheatmap to get ordered of clustered protein and celltypes
#x = pheatmap::pheatmap(mean_mtx, silent = TRUE)
#celltype_order = colnames(mean_mtx[ ,x$tree_col$order]) %>% rev

# convert tidy and reorder based on hclust 
adt.l = adt %>% 
  select(prot_use, celltype) %>% 
  gather(key = prot, value = dsb_count, index1:index2) %>%
  mutate(prot = str_sub(prot, 1, -6)) %>% 
  mutate(prot = factor(prot, levels = rev(prot_use_plot))) %>% 
  mutate(celltype  = factor(celltype, levels = celltype_order))

# plot 
col_split = length(celltype_order) %>% as.numeric()
adt.l = 
  adt.l %>% filter(dsb_count > -5) %>% 
  filter(!celltype=="DOUBLET" )
  
p = ggplot(adt.l, aes(x = dsb_count, y = prot, color = prot, fill = prot)) + 
  geom_density_ridges2(show.legend = F, inherit.aes = T, size = 0.1) + 
  theme_minimal() +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) + 
  geom_vline(xintercept = 0, color = "red", size=0.3) +
  xlab("dsb normalized protein expression") + 
  theme(axis.title.x = element_text(size = 15)) + 
  facet_wrap(~celltype, ncol = col_split, scales = "free_x") + 
  theme(panel.spacing.x = unit(0.1, "lines"))+
  theme(strip.background = element_blank()) +
  theme(strip.text = element_text(colour = 'black', size = 10, angle = 90, hjust = 0)) + 
  theme(axis.text.x = element_text(size = 5,  family = "Helvetica", color = "black")) +
  theme(axis.text.y = element_text(size = 10,  family = "Helvetica", color = "black")) 
ggsave(p, filename = paste0(figpath,"H1_cluster_histogram_heatmap.pdf"), width = 14.5,  height =10)


# R version 3.5.3 Patched (2019-03-11 r77192)
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] viridis_0.5.1     viridisLite_0.3.0 ggsci_2.9         ggridges_0.5.1    magrittr_1.5      forcats_0.4.0     stringr_1.4.0    
# [8] dplyr_0.8.5       purrr_0.3.3       readr_1.3.1       tidyr_1.0.2       tibble_2.1.1      tidyverse_1.2.1   Seurat_2.3.4     
# [15] Matrix_1.2-15     cowplot_0.9.4     ggplot2_3.1.1    
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15          colorspace_1.4-1    ellipsis_0.3.0      class_7.3-15        modeltools_0.2-22   mclust_5.4.5       
# [7] htmlTable_1.13.1    base64enc_0.1-3     rstudioapi_0.10     proxy_0.4-23        npsurv_0.4-0        flexmix_2.3-15     
# [13] bit64_0.9-7         lubridate_1.7.4     mvtnorm_1.0-10      xml2_1.2.0          codetools_0.2-16    splines_3.5.3      
# [19] R.methodsS3_1.7.1   lsei_1.2-0          robustbase_0.93-5   knitr_1.23          Formula_1.2-3       jsonlite_1.6       
# [25] broom_0.5.2         ica_1.0-2           cluster_2.0.7-1     kernlab_0.9-27      png_0.1-7           R.oo_1.22.0        
# [31] pheatmap_1.0.12     compiler_3.5.3      httr_1.4.0          backports_1.1.4     assertthat_0.2.1    lazyeval_0.2.2     
# [37] cli_1.1.0           lars_1.2            acepack_1.4.1       htmltools_0.3.6     tools_3.5.3         igraph_1.2.4.1     
# [43] gtable_0.3.0        glue_1.3.1          RANN_2.6.1          reshape2_1.4.3      Rcpp_1.0.1          cellranger_1.1.0   
# [49] vctrs_0.2.4         gdata_2.18.0        ape_5.3             nlme_3.1-137        iterators_1.0.10    fpc_2.2-1          
# [55] gbRd_0.4-11         lmtest_0.9-37       xfun_0.7            rvest_0.3.4         lifecycle_0.1.0     irlba_2.3.3        
# [61] gtools_3.8.1        DEoptimR_1.0-8      MASS_7.3-51.1       zoo_1.8-6           scales_1.0.0        hms_0.4.2          
# [67] doSNOW_1.0.16       parallel_3.5.3      RColorBrewer_1.1-2  reticulate_1.12     pbapply_1.4-0       gridExtra_2.3      
# [73] rpart_4.1-13        segmented_0.5-4.0   latticeExtra_0.6-28 stringi_1.4.3       foreach_1.4.4       checkmate_1.9.3    
# [79] caTools_1.17.1.2    bibtex_0.4.2        Rdpack_0.11-0       SDMTools_1.1-221.1  rlang_0.4.5         pkgconfig_2.0.2    
# [85] dtw_1.20-1          prabclus_2.3-1      bitops_1.0-6        lattice_0.20-38     ROCR_1.0-7          labeling_0.3       
# [91] htmlwidgets_1.3     bit_1.1-14          tidyselect_0.2.5    plyr_1.8.4          R6_2.4.0            generics_0.0.2     
# [97] snow_0.4-3          gplots_3.0.1.1      Hmisc_4.2-0         haven_2.1.0         pillar_1.4.1        foreign_0.8-71     
# [103] withr_2.1.2         fitdistrplus_1.0-14 mixtools_1.1.0      survival_2.43-3     nnet_7.3-12         tsne_0.1-3         
# [109] modelr_0.1.4        crayon_1.3.4        hdf5r_1.2.0         KernSmooth_2.23-15  readxl_1.3.1        grid_3.5.3         
# [115] data.table_1.12.2   metap_1.1           digest_0.6.19       diptest_0.75-7      R.utils_2.8.0       stats4_3.5.3       
# [121] munsell_0.5.0    