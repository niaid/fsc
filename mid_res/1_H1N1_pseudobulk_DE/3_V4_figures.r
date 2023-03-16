# R version 4.0.5
# H1N1 differential expression testing and gene set enrichment analysis
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))

# set fig path 
figpath = here("mid_res/1_H1N1_pseudobulk_DE/figuresV4/")
datapath = here("mid_res/1_H1N1_pseudobulk_DE/dataV4/")

# heirarchical signal visualization 

# day 1 gsea result 
g1c = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/g1c.rds'))
p = PlotFgsea(gsea_result_list = g1c, p.threshold = 0.05)
#new
mann = data.table::fread(file = here('signature_curation/sig_test_sub_annotation.txt'))
categ = unique(mann$annotation)
pd = p$data 
pd$annotation = plyr::mapvalues(pd$pathway, from = mann$pathway, to = mann$annotation)
pd$annotation = factor(pd$annotation,levels = categ)

# add annotation to plot data env
p$data = pd
p$data$celltype = str_replace_all(p$data$celltype , pattern = '_', replacement  = ' ')
p$data = p$data %>% filter(!pathway == 'btm M4.1 cell cycle (I)')

p$data
# save plot
high.col = ggsci::pal_jama()(2)[2]
low.col = ggsci::pal_jama()(2)[1]
p = p + 
  facet_grid(vars(annotation),  scales = 'free', space = 'free', switch = 'y', ) + 
  theme(strip.text.y = element_text(angle = 0)) +
  theme(strip.placement = "outside") +
  scale_size_area(max_size = 4) +
  theme(legend.position = 'right', legend.justification = 'bottom') +
  theme(strip.background = element_blank()) + 
  scale_fill_gradient2(low = low.col, mid = 'white', high = high.col) + 
  theme(legend.position = 'bottom') 
p
ggsave(p,filename = paste0(figpath, 'g1c.gsea.pdf'),width = 6.3, height = 6.2)




# day 1 heatmap fc 
# highlight genes that have a relative signal enriched in certain cell types 
# define core day 1 induced signature 
fit1 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/fit1.rds'))
d1res = ExtractResult(model.fit.list = fit1, coefficient.number = 1, coef.name = 'L1')

# get leading edge genes of all enrichments by cell type
li = LeadingEdgeIndexed(gsea.result.list = g1c,padj.threshold = 0.05)
li.full = unique(unlist(li, use.names = FALSE))

# logFold Change estimates from mixed model 
gm = GetGeneMatrix(result.list = d1res, 
                   gene_subset = li.full, 
                   pvalfilter = 0.05, 
                   stat_for_matrix = 'logFC', 
                   logfcfilter = 0.1)

# get rid of underscore in celltype name 
colnames(gm) = str_replace_all(colnames(gm), pattern = '_', replacement  = ' ')

# function to define number of columns that have a non zero value
nnzero <- apply(gm, 1, function(x) sum(! x == 0 ))

# Define shared induced state
ms1 = gm[which(nnzero >= 5),]
ms1.names = rownames(ms1)

# save genes defining the shared state
saveRDS(ms1.names, file = paste0(datapath, 'ms1.names.rds'))

# heatmap of core shared day 1 interferon state 
pdf(file = paste0(figpath, 'ms1.mat.pdf'),width = 5, height = 5)
pheatmap::pheatmap(ms1, treeheight_row = 10, treeheight_col = 10,
                   color = viridis::viridis(n = 10, option = "inferno"),clustering_method = 'ward.D',
                   breaks = seq(from = 0, to = 1.5, length.out = 10),
                   fontsize_row = 6,  fontsize_col = 12, border_color = NA)
dev.off()

# celltype specific signals 
asub = c(
  # bc mem and naive 
  'TRIM56','MYD88', 'TBK1', 'CD69',
  # CD14 mono 
  'CAMK2D','IFI27' ,'IL15RA' ,'IL1RN' ,'INSIG1','LCK' ,
  'LYN','SERPING1','SLC16A6' ,'SLC25A1','TGFB1',
  'CCL2' ,'FCGR1B',
  # mdc cd16 cd14
  'WARS' ,'PARP14', 'IFITM2' ,'FBXO6' ,'PSMA4' ,'ACTR2', 'ICAM1', 'IL15', 
  #CD16 mono 
  'LAP3' ,'CYBA' ,'FYB' ,'FCGR1A',
  # CD4 naive 
  'IRF4' ,'IL2RB' ,'IRF8' ,'MXD1' , 'RELA' ,
  # mDC 
  'TYK2' ,'TNFSF10',
  # CD4 EM 
  'IRS2','PTPN2',
  # CD8 naive t 
  'EIF2AK2', 'EIF4G2' ,'PIK3CG',
  # nk 
  'IFI44' ,
  # CD8 CD161 
  'MX2','PIK3R5' , 'CALR'
)


library(ComplexHeatmap)
gmd= gm[! rownames(gm) %in% c(ms1.names), ]
gmd[gmd<0] <- 0
gmd = slanter::slanted_reorder(gmd)

# select labes 
rlab = which(rownames(gmd) %in% asub)
ha = rowAnnotation(foo = anno_mark(
  at = rlab,
  labels = asub,
  labels_gp = gpar(color = "black", fontsize = 8)
))

# rlab2 = which(rownames(gmd) %in% asub2)
# ha2 = rowAnnotation(foo = anno_mark(
#   at = rlab2,
#   labels = asub2,
#   labels_gp = gpar(color = "red", fontsize = 8)
# ))

# color map 
col_fun = circlize::colorRamp2(
  breaks = seq(0,1.2,  length.out = 10),
  colors = viridis::viridis(n = 10, option = 'inferno')
  )

# draw heatmap and rasterize 
pdf(file = paste0(figpath, 'gmd.mat.pdf'),width = 5, height = 7.5)
ComplexHeatmap::Heatmap(matrix = gmd, 
                        right_annotation = ha, 
                        
                        show_row_names = FALSE,
                        # do not cluster use slanter 
                        cluster_columns = FALSE,  cluster_rows = FALSE,
                        col = col_fun, 
                        use_raster = TRUE)
dev.off()


# heirarchical signal deconvolution of reactome interferon signature 
# g1c object loaded in first line above 
# g1c = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/g1c.rds'))

# get all leading edge genes across subsets 
rifn = g1c %>% 
  bind_rows() %>% 
  filter( pathway == 'reactome interferon signaling' & padj < 0.05) %$% leadingEdge %>% 
  unlist() %>% 
  unique()

# extract matrix of log fold change estimate across donors 
mtx = GetGeneMatrix(result.list = d1res, 
                    stat_for_matrix = "logFC", 
                    gene_subset = rifn,
                    pvalfilter = 0.05,
                    logfcfilter = 0.25)
# remove underscores
colnames(mtx) = str_replace_all(colnames(mtx),pattern = '_', replacement = " ")

# draw heatmap 
pheatmap::pheatmap(
  mtx ,
  border_color = NA,
  color = viridis::viridis(n = 18, option = "A"),
  breaks = seq(from = 0, to = 1.5, length.out = 19),
  treeheight_col = 10, treeheight_row = 0,
  fontsize = 6, fontsize_row = 5, fontsize_col = 6,
  width = 2.3, height = 5,
  clustering_method = 'complete',
  filename = paste0(figpath, "d1_reactomeLeadingedge_cluster_genes_heatmap_full.pdf")
)


# now show the main perturbed cell type (mono 14) fold changes in pseudobulk logcpm
# across each individual donor 
# load model fitting data and pseudobulk data saved in mixed model workflow (scipt 1)
d1d = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/d1d.rds'))
d1 = readRDS(file = here("mid_res/1_H1N1_pseudobulk_DE/dataV4/d1.rds"))

# log CPM the pseudobulk data 
lcpm = lapply(d1d, edgeR::cpm, log = TRUE)  

# scglmmr function to extract leading edge genes 
lexp1 = LeadEdgeTidySampleExprs(av.exprs.list = lcpm, gsea.list = g1c, padj.filter = 0.05, NES.filter = 0)

# annotate time and batch on the heatmap
heatmap_anno = d1[, c('batch', 'timepoint')]
anno_color = list(
  timepoint = c('d1' = "orange", "d0" = "white"),
  batch = c('1' = "black", '2' = "white")
)

# define color vector
cu = c("#053061", "#1E61A5", "#3C8ABE", "#7CB7D6", "#BAD9E9", "#E5EEF3", 
       "#F9EAE1", "#F9C7AD", "#EB9273", "#CF5246", "#AB1529", "#67001F")

# scglmmr function for leading edge gene matrix across donors
mat2 = scglmmr::LeadEdgeSampleHeatmap(tidy.exprs.list = lexp1, 
                                 modulename = "reactome interferon signaling",
                                 celltype_plot = 'CD14_Mono',
                                 metadata = meta, 
                                 metadata_annotate = c('adjMFC', 'batch'),
                                 sample_column = 'sample',
                                 returnmat = TRUE)
# draw heatmap 
pheatmap::pheatmap(mat2, 
                   border_color = NA,
                   treeheight_row = 0, treeheight_col = 10,
                   annotation = heatmap_anno,
                   annotation_colors = anno_color,
                   color = cu,
                   width = 5,  height = 7.6,
                   scale = "row",
                   filename = paste0(figpath, "mono_d1_ReactomeIFN.pdf")
                   )

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
