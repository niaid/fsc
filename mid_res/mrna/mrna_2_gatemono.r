# R 4.0.5
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(Matrix))
suppressMessages(library(ggsci))
library(magrittr)
library(dsb)

# set paths
datapath = file.path(here('mid_res/mrna/generated_data/'))
figpath = file.path(here('mid_res/mrna/figures/'))

# load data 
rna = readRDS(file = paste0(datapath,'rna.rds'))
md = readRDS(file = paste0(datapath,'md.rds'))
adt = readRDS(file = paste0(datapath,'adt.rds'))

#slim 
rna <- as(object = rna, Class = "dgCMatrix")
adt <- as(object = adt, Class = "dgCMatrix")

# seurat workflow 
s = CreateSeuratObject(counts = rna, min.cells = 20,  meta.data = md)

# normalize ADT with dsb function ModelNegativeADTnorm
iso = c("Isotype1_ADT", "Isotype2_ADT", "Isotype3_ADT", "Isotype4_ADT")
adt_norm = ModelNegativeADTnorm(cell_protein_matrix = adt, 
                                denoise.counts = TRUE, 
                                use.isotype.control = TRUE, 
                                isotype.control.name.vec = iso, 
                                quantile.clipping = TRUE, 
                                return.stats = TRUE)

# define phenotyping antibodies
rownames(adt_norm$dsb_normalized_matrix) = 
  str_replace_all(
    string = rownames(adt_norm$dsb_normalized_matrix),
    pattern = '_', replacement = '-'
    )

s[['CITE']] = CreateAssayObject(counts = adt)
s = SetAssayData(object = s,
                 slot = 'data', 
                 new.data = adt_norm$dsb_normalized_matrix,
                 assay = 'CITE')

# save processed object 
saveRDS(s,file = paste0(datapath, 's.rds'))

# gate out monocytes 
library(scales)
d = cbind(s@meta.data, data.frame(t(as.matrix(s@assays$CITE@data))))
p = ggplot(d, aes(x = CD14.ADT, y = CD3.ADT)) + geom_point(size  =0.1, alpha = 0.4) +
  scale_y_continuous( breaks=pretty_breaks() )  + 
  scale_x_continuous( breaks=pretty_breaks() ) + 
  geom_abline(slope = 1,intercept = -0.9, color = 'red') + 
  geom_vline(xintercept = 1.5, color = 'red') + 
  geom_hline(yintercept = 1.5, color = 'red')

# manually create triangular gate 
d$tx = d$CD14.ADT > 1.5
d$ty = d$CD3.ADT < 1.5
d$tlm = d$CD14*1 + -0.9
# add gate info 
d$pp3 = ifelse(d$tx==TRUE & d$ty==TRUE & d$CD3.ADT < d$tlm, yes = '1',no = '0')

# plot with gated cells highlighted 
p = ggplot(d, aes(x = CD14.ADT, y = CD3.ADT, color = pp3)) + 
  geom_point(size  =0.1, alpha = 0.4) +
  theme_bw() + 
  scale_y_continuous( breaks=pretty_breaks())  + 
  scale_x_continuous( breaks=pretty_breaks()) + 
  geom_abline(slope = 1,intercept = -0.9, color = 'red') + 
  geom_vline(xintercept = 1.5, color = 'red') + 
  geom_hline(yintercept = 1.5, color = 'red') + 
  ylab('CD3 ADT dsb:ModelNegative') + xlab('CD14 ADT dsb::ModelNegative')
p
ggsave(p,filename = paste0(figpath,'monogate.png'), width = 9, height = 8)

# define monocytes
dmono = d[d$pp3==1, ] %>% rownames() 
saveRDS(dmono,file = paste0(datapath, 'dmono.rds'))

# subset and save monocyte Seurat object. 
s.mono  = subset(s,cells = dmono)
saveRDS(s.mono,file = paste0(datapath, 's.mono.rds'))


# gate out mdc 
p = ggplot(d, aes(x = CD11c.ADT, y = CD1c.BDCA1.ADT)) + 
  geom_point(size  =0.1, alpha = 0.4) +
  geom_vline(xintercept = 2.2, color = 'red') + 
  geom_hline(yintercept = 1.5, color = 'red')

# define mDC 
d$mdc = ifelse(d$CD11c.ADT>2.2 & d$CD1c.BDCA1.ADT>1.5, yes = '1',no = '0')

#plot gated cells  
p = ggplot(d, aes(x = CD11c.ADT, y = CD1c.BDCA1.ADT, color = mdc)) + 
  theme_bw() + 
  geom_point(size  =0.1, alpha = 0.4) +
  geom_vline(xintercept = 2.2, color = 'red') + 
  geom_hline(yintercept = 1.5, color = 'red') 
ggsave(p,filename = paste0(figpath,'mdc_gate.png'), width = 9, height = 8)


# subst mDC 
mdc.cells = d[d$mdc=="1", ] %>% rownames()

# subset mDC 
s.mdc = subset(s,cells = mdc.cells)
saveRDS(s.mdc,file = paste0(datapath, 's.mdc.rds'))

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] magrittr_2.0.1     scales_1.1.1       plotly_4.9.3       shiny_1.6.0       
# [5] dsb_1.0.2          ggsci_2.9          Matrix_1.3-2       here_1.0.1        
# [9] forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4       
# [13] readr_1.4.0        tidyr_1.1.2        tibble_3.0.6       ggplot2_3.3.3     
# [17] tidyverse_1.3.0    SeuratObject_4.0.0 Seurat_4.0.1      
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-0      deldir_1.0-6         
# [4] ellipsis_0.3.2        ggridges_0.5.3        rsconnect_0.8.25     
# [7] mclust_5.4.7          rprojroot_2.0.2       fs_1.5.0             
# [10] rstudioapi_0.13       spatstat.data_2.1-0   farver_2.0.3         
# [13] leiden_0.3.7          listenv_0.8.0         ggrepel_0.9.1        
# [16] lubridate_1.7.9.2     xml2_1.3.2            codetools_0.2-18     
# [19] splines_4.0.5         cachem_1.0.4          polyclip_1.10-0      
# [22] jsonlite_1.7.2        packrat_0.7.0         broom_0.7.5          
# [25] ica_1.0-2             cluster_2.1.2         dbplyr_2.1.0         
# [28] png_0.1-7             pheatmap_1.0.12       uwot_0.1.10          
# [31] sctransform_0.3.2     spatstat.sparse_2.0-0 compiler_4.0.5       
# [34] httr_1.4.2            backports_1.2.1       assertthat_0.2.1     
# [37] fastmap_1.1.0         lazyeval_0.2.2        limma_3.46.0         
# [40] cli_3.3.0             later_1.1.0.1         htmltools_0.5.2      
# [43] tools_4.0.5           igraph_1.2.6          gtable_0.3.0         
# [46] glue_1.6.2            RANN_2.6.1            reshape2_1.4.4       
# [49] Rcpp_1.0.6            scattermore_0.7       jquerylib_0.1.3      
# [52] cellranger_1.1.0      vctrs_0.4.1           nlme_3.1-152         
# [55] crosstalk_1.1.1       lmtest_0.9-38         globals_0.14.0       
# [58] rvest_0.3.6           mime_0.10             miniUI_0.1.1.1       
# [61] lifecycle_1.0.0       irlba_2.3.3           goftest_1.2-2        
# [64] FSA_0.9.0             future_1.21.0         MASS_7.3-53.1        
# [67] zoo_1.8-8             spatstat.core_2.0-0   hms_1.0.0            
# [70] promises_1.2.0.1      spatstat.utils_2.3-0  parallel_4.0.5       
# [73] RColorBrewer_1.1-2    yaml_2.2.1            reticulate_1.18      
# [76] pbapply_1.4-3         gridExtra_2.3         sass_0.4.0           
# [79] rpart_4.1-15          stringi_1.5.3         rlang_1.0.2          
# [82] pkgconfig_2.0.3       matrixStats_0.58.0    lattice_0.20-41      
# [85] ROCR_1.0-11           tensor_1.5            labeling_0.4.2       
# [88] patchwork_1.1.1       htmlwidgets_1.5.3     cowplot_1.1.1        
# [91] tidyselect_1.1.0      parallelly_1.23.0     RcppAnnoy_0.0.18     
# [94] plyr_1.8.6            R6_2.5.0              generics_0.1.2       
# [97] DBI_1.1.1             withr_2.4.3           pillar_1.4.7         
# [100] haven_2.3.1           mgcv_1.8-34           fitdistrplus_1.1-3   
# [103] survival_3.2-10       abind_1.4-5           future.apply_1.7.0   
# [106] modelr_0.1.8          crayon_1.4.1          KernSmooth_2.23-18   
# [109] spatstat.geom_2.4-0   viridis_0.5.1         grid_4.0.5           
# [112] readxl_1.3.1          data.table_1.14.0     reprex_1.0.0         
# [115] digest_0.6.27         xtable_1.8-4          httpuv_1.5.5         
# [118] munsell_0.5.0         viridisLite_0.3.0     bslib_0.3.1  
 
