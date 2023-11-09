# R 4.0.5
# testing out mrna data 
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(Matrix))

# save paths 
datapath = file.path(here('mid_res/mrna/generated_data/')) 
dir.create(datapath)

# read sparse matrix 
mtx = Matrix::readMM(file = here('data/GSE171964/matrix.mtx'))

# reformat because not formatted in geo for Seurat::Read10X()
cells = read.delim(file = here('data/GSE171964/barcodes.tsv'))
cells = cells %>% separate(x, into = c('space', 'barcode'),sep =' ')
cells = cells$barcode
saveRDS(cells, file = paste0(datapath,'cells.rds'))
cells = readRDS(file = here('mid_res/mrna/generated_data/cells.rds'))


# reformat features 
features = read.delim(file = here('data/GSE171964/features.tsv'))
features = features %>% separate(x, into = c('space', 'feature'),sep =' ')
features = features$feature
saveRDS(features, file = paste0(datapath,'features.rds'))
features = readRDS(file = here('mid_res/mrna/generated_data/features.rds'))

# set names of matrix 
colnames(mtx) = cells
rownames(mtx) = features

# load metadata - updated author correction data 
p = read.delim(file = here('data/GSE171964/GSE171964_geo_pheno_v2.csv'), sep = ',')

# define the day 0 day 1 cells and format for sc meta
psub = p %>% filter(day %in% c('0', '1', '21','22'))
cell_sub = psub$barcode
md = psub %>% column_to_rownames('barcode')

# subset matrix 
mtx = mtx[ ,cell_sub]

# get ADT data 
proteins = features[grep(pattern = '_ADT',x = features)]
adt = mtx[proteins, ]
rna = mtx[rownames(mtx)[!rownames(mtx) %in% proteins], ]

# save
saveRDS(rna, file = paste0(datapath,'rna.rds'))
saveRDS(adt, file = paste0(datapath,'adt.rds'))
saveRDS(md, file = paste0(datapath,'md.rds'))

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] Matrix_1.3-2       here_1.0.1         forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4       
# [6] purrr_0.3.4        readr_1.4.0        tidyr_1.1.2        tibble_3.0.6       ggplot2_3.3.3     
# [11] tidyverse_1.3.0    SeuratObject_4.0.0 Seurat_4.0.1      
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-0      deldir_0.2-10         ellipsis_0.3.1       
# [5] ggridges_0.5.3        rprojroot_2.0.2       fs_1.5.0              rstudioapi_0.13      
# [9] spatstat.data_2.1-0   leiden_0.3.7          listenv_0.8.0         ggrepel_0.9.1        
# [13] lubridate_1.7.9.2     xml2_1.3.2            codetools_0.2-18      splines_4.0.5        
# [17] polyclip_1.10-0       jsonlite_1.7.2        broom_0.7.5           ica_1.0-2            
# [21] cluster_2.1.2         dbplyr_2.1.0          png_0.1-7             uwot_0.1.10          
# [25] shiny_1.6.0           sctransform_0.3.2     spatstat.sparse_2.0-0 compiler_4.0.5       
# [29] httr_1.4.2            backports_1.2.1       assertthat_0.2.1      fastmap_1.1.0        
# [33] lazyeval_0.2.2        cli_2.5.0             later_1.1.0.1         htmltools_0.5.1.1    
# [37] tools_4.0.5           igraph_1.2.6          gtable_0.3.0          glue_1.4.2           
# [41] RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.6            scattermore_0.7      
# [45] cellranger_1.1.0      vctrs_0.3.6           nlme_3.1-152          lmtest_0.9-38        
# [49] globals_0.14.0        rvest_0.3.6           mime_0.10             miniUI_0.1.1.1       
# [53] lifecycle_1.0.0       irlba_2.3.3           goftest_1.2-2         future_1.21.0        
# [57] MASS_7.3-53.1         zoo_1.8-8             scales_1.1.1          spatstat.core_2.0-0  
# [61] hms_1.0.0             promises_1.2.0.1      spatstat.utils_2.1-0  parallel_4.0.5       
# [65] RColorBrewer_1.1-2    reticulate_1.18       pbapply_1.4-3         gridExtra_2.3        
# [69] rpart_4.1-15          stringi_1.5.3         rlang_0.4.10          pkgconfig_2.0.3      
# [73] matrixStats_0.58.0    lattice_0.20-41       ROCR_1.0-11           tensor_1.5           
# [77] patchwork_1.1.1       htmlwidgets_1.5.3     cowplot_1.1.1         tidyselect_1.1.0     
# [81] parallelly_1.23.0     RcppAnnoy_0.0.18      plyr_1.8.6            magrittr_2.0.1       
# [85] R6_2.5.0              generics_0.1.0        DBI_1.1.1             withr_2.4.1          
# [89] pillar_1.4.7          haven_2.3.1           mgcv_1.8-34           fitdistrplus_1.1-3   
# [93] survival_3.2-10       abind_1.4-5           future.apply_1.7.0    modelr_0.1.8         
# [97] crayon_1.4.1          KernSmooth_2.23-18    spatstat.geom_2.0-1   plotly_4.9.3         
# [101] grid_4.0.5            readxl_1.3.1          data.table_1.14.0     reprex_1.0.0         
# [105] digest_0.6.27         xtable_1.8-4          httpuv_1.5.5          munsell_0.5.0        
# [109] viridisLite_0.3.0    

