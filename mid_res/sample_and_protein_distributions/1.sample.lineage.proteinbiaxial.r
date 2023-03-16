## Must be run in R 3.5.1 
# umap of joint clustering results 
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(here))
source("functions/analysis_functions.R")

# Set path 
figpath = here("mid_res/sample_and_protein_distributions/figures/"); 
dir.create(figpath, recursive = TRUE)

# full sample bar plot 
md = readRDS(file = here("data/h1h5_annotated_with_meta.rds"))@meta.data
celltypes = md$celltype_joint %>% unique() %>% sort()
t4 = celltypes[c(7:12)]
t8 = celltypes[c(13:16)]
myeloid = celltypes[c(4:5, 18, 19, 21, 23)]
bc = celltypes[c(1,2,6)]
nk = celltypes[c(22)]
unconventionalT = celltypes[c(3,20)]

md = md %>% mutate(lineage = 
            if_else(celltype_joint %in% t4, "CD4 T Cell",
            if_else(celltype_joint %in% t8, "CD8 T cell",
            if_else(celltype_joint %in% myeloid, "Myeloid", 
            if_else(celltype_joint %in% bc, "B cell", 
            if_else(celltype_joint %in% nk, "NK cell",
            if_else(celltype_joint %in% unconventionalT, "unconventional T",
                    false = "other")))))))
md2 = md %>% filter(!celltype_joint %in% "DOUBLET")

# calc and vis fraction of total 
d = md2 %>% group_by(lineage, sample) %>% tally

p = 
  ggplot(d, aes(x = sample, y = n, fill = lineage)) + 
  geom_bar(position = 'fill', stat = 'identity',  show.legend = TRUE) +
  ggsci::scale_fill_jama() + 
  ylab("percent of total") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + 
  theme(axis.text.x = element_text(size =6)) + 
  theme(axis.title.x = element_blank())
ggsave(p, filename = paste0(figpath, "LINEAGE_fullsamplebarplot.pdf"), width = 9, height = 4)

#### #Manual gate plots 
figpath = paste0(figpath,"mgplots/") ; dir.create(figpath)

# Day 1 cohort CD14 monocyte data
cite = as.data.frame(t(readRDS(file = here("h1h5_annotated.rds"))@assay$CITE@data))
colnames(cite) = str_sub(colnames(cite), 1, -6)
mdf = cbind(md,cite)
h1md = mdf %>% filter(cohort == "H1N1")

# match colors in umap 
celltypes = readRDS(here("data/celltypes_vector_ordered.rds"))
cu = pals::kelly(n = 22) %>% as.vector()
cu = cu[-1]
cu = c("midnightblue", cu, "lightslateblue")
cu[15] = "maroon1"
cu[11] = "darkseagreen1"
cu = cu[-1]
cu = rev(cu)
h1md$celltype_joint = factor(h1md$celltype_joint, levels = celltypes)

################### 
# B cells 
p = 
  ggplot(h1md %>% filter(lineage %in% "B cell"), aes(x = CD27, CD38, color = celltype_joint)) +
  theme_bw(base_size = 12) + 
  scale_color_manual(values = c("lightslateblue","#2B3D26", "#882D17")) + 
  geom_density_2d() + 
  labs(color="celltype") + 
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) + 
  theme(legend.position = c(0.85, 0.4)) + 
  theme(legend.key.size = unit(0.35, units='cm'))
ggsave(p, filename = paste0(figpath,'Bcell.pdf'), width = 3.5, height = 3)

p = 
  ggplot(h1md %>% filter(lineage %in% "B cell"), aes(x = CD27, CD38)) + 
  theme_bw(base_size = 12) + 
  geom_bin2d(bins = 200) + 
  scale_fill_viridis(option = "B") + 
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  theme(legend.position = c(0.85, 0.4)) + 
  theme(legend.key.size = unit(0.45, units='cm'))
ggsave(p, filename = paste0(figpath,'Bcell_density.pdf'), width = 3.5, height = 3)



################### 
### Myeloid Cells 
p = 
  ggplot(h1md %>% filter(lineage %in% "Myeloid") %>% filter(!celltype_joint %in% 'IgA_CD14_Mono'),
         aes(x = CD16, CD303, color = celltype_joint)) + 
  theme_bw(base_size = 12) + 
  scale_color_manual(values = c("#654522", "#8DB600", "#BE0032","#875692","#222222"))+
  geom_density_2d() + 
  labs(color="celltype") + 
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  theme(legend.position = c(0.75, 0.6)) 
p
ggsave(p, filename = paste0(figpath,'myeloid1.pdf'), width = 3.5, height = 3)

p = 
  ggplot(h1md %>% filter(lineage %in% "Myeloid") %>% filter(!celltype_joint %in% 'IgA_CD14_Mono'),
         aes(x = CD16, CD303)) + 
  theme_bw(base_size = 12) + 
  geom_bin2d(bins = 200) + 
  scale_fill_viridis(option = "B") + 
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  theme(legend.position = c(0.75, 0.6)) + 
  theme(legend.key.size = unit(0.45, units = 'cm'))
p
ggsave(p, filename = paste0(figpath,'myeloid1_density.pdf'), width = 3.5, height = 3)

# prots 2 
p = 
  ggplot(h1md %>% filter(lineage %in% "Myeloid") %>% filter(!celltype_joint %in% 'IgA_CD14_Mono'),
         aes(x = CD16, CD14, color = celltype_joint)) + 
  theme_bw(base_size = 12) + 
  scale_color_manual(values = c("#654522", "#8DB600", "#BE0032","#875692","#222222"))+
  geom_density_2d() + 
  labs(color="celltype") + 
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  theme(legend.position = c(0.8, 0.75))  + 
  theme(legend.key.size = unit(0.45, units = 'cm'))
p
ggsave(p, filename = paste0(figpath,'myeloid2.pdf'), width = 3.5, height = 3)

# 
# ### T clels 
p = ggplot(h1md %>% filter(lineage %in% "CD8 T cell"), aes(x = CD161, CD45RO, color = celltype_joint)) +
  theme_bw(base_size = 12) +
  ggsci::scale_color_jco() +
  geom_density_2d() + 
  labs(color="celltype") + 
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) 
p
ggsave(p, filename = paste0(figpath,'cd8Tcell.pdf'), width = 5, height = 3)



# R version 3.5.3 Patched (2019-03-11 r77192)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] viridis_0.5.1     viridisLite_0.3.0 here_0.1          forcats_0.4.0     stringr_1.4.0     dplyr_0.8.5       purrr_0.3.3       readr_1.3.1       tidyr_1.0.2       tibble_2.1.1      tidyverse_1.2.1  
# [12] Seurat_2.3.4      Matrix_1.2-15     cowplot_0.9.4     ggplot2_3.1.1    
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1            snow_0.4-3              backports_1.1.4         Hmisc_4.2-0             plyr_1.8.4              igraph_1.2.4.1          lazyeval_0.2.2          splines_3.5.3          
# [9] crosstalk_1.0.0         digest_0.6.25           foreach_1.4.4           htmltools_0.3.6         lars_1.2                gdata_2.18.0            magrittr_2.0.1          checkmate_1.9.3        
# [17] cluster_2.0.7-1         mixtools_1.1.0          ROCR_1.0-7              modelr_0.1.4            R.utils_2.8.0           colorspace_1.4-1        rvest_0.3.4             haven_2.1.0            
# [25] xfun_0.7                crayon_1.3.4            jsonlite_1.6            survival_2.43-3         zoo_1.8-6               iterators_1.0.10        ape_5.3                 glue_1.3.1             
# [33] pals_1.5                gtable_0.3.0            webshot_0.5.1           kernlab_0.9-27          prabclus_2.3-1          DEoptimR_1.0-8          maps_3.3.0              scales_1.0.0           
# [41] mvtnorm_1.0-10          bibtex_0.4.2            miniUI_0.1.1.1          Rcpp_1.0.1              metap_1.1               dtw_1.20-1              xtable_1.8-4            htmlTable_1.13.1       
# [49] reticulate_1.12         foreign_0.8-71          bit_1.1-14              mapproj_1.2.6           proxy_0.4-23            mclust_5.4.5            SDMTools_1.1-221.1      Formula_1.2-3          
# [57] stats4_3.5.3            tsne_0.1-3              htmlwidgets_1.3         httr_1.4.0              gplots_3.0.1.1          RColorBrewer_1.1-2      fpc_2.2-1               acepack_1.4.1          
# [65] modeltools_0.2-22       ica_1.0-2               pkgconfig_2.0.2         R.methodsS3_1.7.1       flexmix_2.3-15          nnet_7.3-12             manipulateWidget_0.10.0 tidyselect_0.2.5       
# [73] labeling_0.3            rlang_0.4.5             reshape2_1.4.3          later_0.8.0             munsell_0.5.0           cellranger_1.1.0        tools_3.5.3             cli_1.1.0              
# [81] generics_0.0.2          broom_0.5.2             ggridges_0.5.1          npsurv_0.4-0            knitr_1.23              bit64_0.9-7             fitdistrplus_1.0-14     robustbase_0.93-5      
# [89] rgl_0.100.30            caTools_1.17.1.2        RANN_2.6.1              packrat_0.5.0           pbapply_1.4-0           nlme_3.1-137            mime_0.6                R.oo_1.22.0            
# [97] xml2_1.2.0              hdf5r_1.2.0             compiler_3.5.3          rstudioapi_0.10         png_0.1-7               lsei_1.2-0              stringi_1.4.3           lattice_0.20-38        
# [105] ggsci_2.9               vctrs_0.2.4             pillar_1.4.1            lifecycle_0.1.0         Rdpack_0.11-0           lmtest_0.9-37           data.table_1.12.2       bitops_1.0-6           
# [113] irlba_2.3.3             gbRd_0.4-11             httpuv_1.5.1            R6_2.4.0                latticeExtra_0.6-28     promises_1.0.1          KernSmooth_2.23-15      gridExtra_2.3          
# [121] codetools_0.2-16        dichromat_2.0-0         MASS_7.3-51.1           gtools_3.8.1            assertthat_0.2.1        rprojroot_1.3-2         withr_2.1.2             diptest_0.75-7         
# [129] parallel_3.5.3          doSNOW_1.0.16           hms_0.4.2               grid_3.5.3              rpart_4.1-13            class_7.3-15            segmented_0.5-4.0       Rtsne_0.15             
# [137] shiny_1.3.2             lubridate_1.7.4         base64enc_0.1-3  

