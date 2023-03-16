# Analysis of plasmablast and Activated_B.
# script uses R 3.5.1 
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(here))

# source functions
theme_set(theme_bw())
source("functions/analysis_functions.R")

# file path
figpath = here("mid_res/pblast_abc_integration/figures/")
datapath = here("mid_res/pblast_abc_integration/generated_data/")
dir.create(figpath); dir.create(datapath)

# load h1 annotated data 
h1 = ReadCohort(joint_object_dir = "data/h1h5_annotated_with_meta.rds", cohort = "H1N1") %>% 
  SubsetData(accept.high = 28, subset.name = "CD3_PROT", subset.raw = T)

#manually gate activated memory B cells and plasmablasts
p = GenePlot4(h1, gene1 = "CD19_PROT", gene2 = "CD14_PROT", pt.size = 0.1)
ggsave(p, filename = paste0(figpath,"/cd19cells.png"),width = 4, height = 3)
p = GenePlot4(h1, gene1 = "CD19_PROT", gene2 = "CD3_PROT", pt.size = 0.1) + 
  geom_vline(xintercept = 8) +
  geom_hline(yintercept = 5)
ggsave(p, filename = paste0(figpath,"/cd19cells_2.pdf"),width = 4, height = 3)

############## Pt 1 Manual gate asc abc 
GateBC = function(SeuratObject, return.seurat = F) {
  adt = as.matrix(SeuratObject@assay$CITE@data)
  cells <- names(which(adt["CD19_PROT", ] > 8 &
                         adt["CD3_PROT", ] < 5 &
                         adt["CD14_PROT", ] < 5  ))
  if(return.seurat == TRUE) {
    sub = SubsetData(SeuratObject, cells.use = cells, subset.raw = TRUE)
    return(sub)
  } else { return(cells) }
}

bcells_gate = GateBC(SeuratObject = h1, return.seurat = F)
cd19 = SubsetData(h1, cells.use = bcells_gate, subset.raw = T)

# abc asc 
p = GenePlot4(cd19, gene1 = "CD71_PROT", gene2 = "IgD_PROT",pt.size = 0.1) + 
 geom_vline(xintercept = 5) + geom_hline(yintercept = 10)
ggsave(p, filename =  paste0(figpath,"/cd71igd_cells.pdf"),width = 4, height = 3)

# naive memory 
p = GenePlot4(cd19, gene1 = "CD27_PROT", gene2 = "IgD_PROT",pt.size = 0.1) + 
   geom_vline(xintercept = 4) + geom_hline(yintercept = 10)
ggsave(p, filename =  paste0(figpath,"/cd27_igd_cells.pdf"),width = 4, height = 3)

# activated B and asc gate 
Gate_Activated_BASC =  function(SeuratObject) {
  adt = as.matrix(SeuratObject@assay$CITE@data)
  cells <- names(which(adt["CD71_PROT", ] > 5 & adt["IgD_PROT", ] < 10))
  return(cells) 
}
Activated_B_asc_cells = Gate_Activated_BASC(SeuratObject = cd19)
Activated_Basc =  SubsetData(cd19, cells.use = Activated_B_asc_cells, subset.raw = T)

# plot gates
p = GenePlot4(Activated_Basc, gene1 = "CD20_PROT", gene2 = "CD38_PROT", pt.size = 0.6) +
   geom_vline(xintercept = 8) + 
  geom_hline(yintercept = 10)
ggsave(p, filename =  paste0(figpath,"/Activated_B_pblast.pdf"),width = 4, height = 3)

bmd = cbind(Activated_Basc@meta.data, as.data.frame(t(Activated_Basc@assay$CITE@data)))
bmd = bmd %>%  mutate(
  b_type =
    if_else(
      CD19_PROT > 8 &
        CD71_PROT > 5 &
        IgD_PROT  < 10 &
        CD20_PROT < 8 & CD38_PROT > 10,
      true = "Plasmablast",
      if_else(
        CD19_PROT > 8 &
          CD71_PROT > 5 &
          IgD_PROT  < 10 &
          CD20_PROT > 8 &
          CD38_PROT  < 10,
        true = "Activated_Bcell",
        false = "NA"
      )
    )
) %>%
  filter(b_type %in% c("Plasmablast", "Activated_Bcell")) %>%
  select(b_type, barcode_check) %>%
  column_to_rownames("barcode_check")
bsub = SubsetData(Activated_Basc, cells.use = rownames(bmd), subset.raw = TRUE) %>% AddMetaData(metadata = bmd)

########## Pt 2 add module scores for ellebedy gene sets. 
bcgenes = read.table(file = here("signature_curation/ellebedy_genes.txt"), sep = "\t", header = T)
Activated_B.genes = bcgenes %>% filter(celltype == "ABC-")
asc.genes = bcgenes %>% filter(celltype == "ASC-")
module.list = list(as.character(Activated_B.genes$Gene), as.character(asc.genes$Gene))
names(module.list) = c("Activated_B_module", "ASC_module")
saveRDS(module.list, file =  paste0(datapath,"/ellebedy_bcell.rds"))

# pt 2 module s
# add module score for ellebedy genes 
bsub = AddModuleScore(bsub, 
                      genes.list = module.list, 
                      enrich.name = names(module.list), 
                      random.seed = 1)

names(bsub@meta.data)[c(33, 34)] = c("Activated_Bcell_Gene_Score", "Plasmablast_Gene_score")
bsubdf = bsub@meta.data %>% select(Activated_Bcell_Gene_Score, Plasmablast_Gene_score, b_type)

# plot module scores. 
bsubdf = bsub@meta.data 
p = ggpubr::ggviolin(data = bsubdf, x = "b_type", 
                     y = c("Plasmablast_Gene_score", "Activated_Bcell_Gene_Score"), 
                     combine = TRUE, 
                     fill = "b_type", 
                     palette = "d3") 
p = p %>% ggpubr::ggadd(add = "jitter", jitter = 0.35, alpha = 0.4, size = 1, shape = 16)  
p = p +
  theme(legend.position = "none") +
  ylab("module score") + xlab("") + 
  theme(strip.background = element_blank()) + 
  ggpubr::stat_compare_means(method = "wilcox")
ggsave(p, filename =paste0(figpath,"/pb_asc_modules.pdf"), height = 3, width = 4.5)

sessionInfo()
# R version 3.5.3 Patched (2019-03-11 r77192)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] viridis_0.5.1     viridisLite_0.3.0 here_0.1          forcats_0.4.0     stringr_1.4.0     dplyr_0.8.5      
# [7] purrr_0.3.3       readr_1.3.1       tidyr_1.0.2       tibble_2.1.1      tidyverse_1.2.1   Seurat_2.3.4     
# [13] Matrix_1.2-15     cowplot_0.9.4     ggplot2_3.1.1    
# 
# loaded via a namespace (and not attached):
# [1] Rtsne_0.15          colorspace_1.4-1    class_7.3-15        modeltools_0.2-22   ggridges_0.5.1      rprojroot_1.3-2    
# [7] mclust_5.4.5        htmlTable_1.13.1    base64enc_0.1-3     rstudioapi_0.10     proxy_0.4-23        ggpubr_0.2         
# [13] npsurv_0.4-0        flexmix_2.3-15      bit64_0.9-7         lubridate_1.7.4     mvtnorm_1.0-10      xml2_1.2.0         
# [19] codetools_0.2-16    splines_3.5.3       R.methodsS3_1.7.1   lsei_1.2-0          robustbase_0.93-5   knitr_1.23         
# [25] Formula_1.2-3       jsonlite_1.6        packrat_0.5.0       broom_0.5.2         ica_1.0-2           cluster_2.0.7-1    
# [31] kernlab_0.9-27      png_0.1-7           R.oo_1.22.0         compiler_3.5.3      httr_1.4.0          backports_1.1.4    
# [37] assertthat_0.2.1    lazyeval_0.2.2      cli_1.1.0           lars_1.2            acepack_1.4.1       htmltools_0.3.6    
# [43] tools_3.5.3         igraph_1.2.4.1      gtable_0.3.0        glue_1.3.1          RANN_2.6.1          reshape2_1.4.3     
# [49] Rcpp_1.0.1          cellranger_1.1.0    vctrs_0.2.4         gdata_2.18.0        ape_5.3             nlme_3.1-137       
# [55] iterators_1.0.10    fpc_2.2-1           gbRd_0.4-11         lmtest_0.9-37       xfun_0.7            rvest_0.3.4        
# [61] lifecycle_0.1.0     irlba_2.3.3         gtools_3.8.1        DEoptimR_1.0-8      MASS_7.3-51.1       zoo_1.8-6          
# [67] scales_1.0.0        hms_0.4.2           doSNOW_1.0.16       parallel_3.5.3      RColorBrewer_1.1-2  reticulate_1.12    
# [73] pbapply_1.4-0       gridExtra_2.3       rpart_4.1-13        segmented_0.5-4.0   latticeExtra_0.6-28 stringi_1.4.3      
# [79] foreach_1.4.4       checkmate_1.9.3     caTools_1.17.1.2    bibtex_0.4.2        Rdpack_0.11-0       SDMTools_1.1-221.1 
# [85] rlang_0.4.5         pkgconfig_2.0.2     dtw_1.20-1          prabclus_2.3-1      bitops_1.0-6        lattice_0.20-38    
# [91] ROCR_1.0-7          labeling_0.3        htmlwidgets_1.3     bit_1.1-14          tidyselect_0.2.5    plyr_1.8.4         
# [97] magrittr_2.0.1      R6_2.4.0            generics_0.0.2      snow_0.4-3          gplots_3.0.1.1      Hmisc_4.2-0        
# [103] haven_2.1.0         pillar_1.4.1        foreign_0.8-71      withr_2.1.2         fitdistrplus_1.0-14 mixtools_1.1.0     
# [109] survival_2.43-3     nnet_7.3-12         tsne_0.1-3          modelr_0.1.4        crayon_1.3.4        hdf5r_1.2.0        
# [115] KernSmooth_2.23-15  readxl_1.3.1        grid_3.5.3          data.table_1.12.2   metap_1.1           digest_0.6.25      
# [121] diptest_0.75-7      R.utils_2.8.0       stats4_3.5.3        munsell_0.5.0  