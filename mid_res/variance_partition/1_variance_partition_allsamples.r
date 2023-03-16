# R4 
# initialize 
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
library(BiocParallel)
library(variancePartition)
library(scglmmr)

register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# figpath
figpath = here('mid_res/variance_partition/figures/'); dir.create(figpath, recursive = TRUE)
datapath = here('mid_res/variance_partition/generated_data/'); dir.create(datapath, recursive = TRUE)

# load data 
s = readRDS(file = here('data/h1h5_annotated_with_meta.rds'))
table(s@meta.data$time_cohort, s@meta.data$sampleid)
table(s@meta.data$timepoint, s@meta.data$sampleid)

# pb data 
meta = s@meta.data
umi = s@raw.data
pdf()
tab = scglmmr::SubjectCelltypeTable(metadata = meta, celltype_column = "celltype_joint", sample_column = "sample")

# remove cells prior to pseudobulk analysis 
meta = meta[!meta$celltype_joint %in% c(tab$celltypes_remove, 'DOUBLET'), ]
umi = umi[ ,rownames(meta)]

# make pseudobulk data 
pb = scglmmr::PseudobulkList(rawcounts = umi,
                             metadata = meta, 
                             sample_col = "sample",
                             celltype_col = "celltype_joint", 
                             avg_or_sum = "sum")

# add cell type to column names of pb list 
for (i in 1:length(pb)) {
  colnames(pb[[i]]) = paste(colnames(pb[[i]]), names(pb[i]), sep = "~")
}
saveRDS(pb,file = paste0(figpath, 'pb_vp.rds'))
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))

# merge pseudobulk data 
pd = do.call(cbind, pb)

# make celltype ~ sample metadata 
csd = 
  data.frame(sid = colnames(pd)) %>% 
  mutate(sample_celltype = sid) %>% 
  separate(sid, into = c('sample', 'celltype'), sep = '~') 
  #column_to_rownames('sample_celltype') 

# make subject metadata
samplemd = 
  meta %>% 
  group_by(sample) %>% 
  select(age, sampleid, gender, batch, time_cohort, timepoint, adjmfc.group) %>% 
  summarise_each(funs = unique) %>% 
  ungroup() %>% 
  as.data.frame()
saveRDS(samplemd, file = paste0(here('data/samplemd.rds')))

# sample_celltype metadata for design matrices 
cf = full_join(csd, samplemd, by = 'sample') %>%
  column_to_rownames('sample_celltype')


#############
# process bulk data 
#filter lowly expressed (in this case basically unexpressed genes)
gene_keep = edgeR::filterByExpr(pd, 
                                min.count = 1,
                                min.total.count = 1,
                                min.prop = 0.5,
                                group = as.factor(cf$celltype))
# FALSE  TRUE 
# 5314 14320 
pd = pd[gene_keep, ]

# normalize bulk data 
pd = edgeR::DGEList(counts = pd, samples = cf)
pd = edgeR::calcNormFactors(object = pd)


##############
# Get voom observational weights 
# these precision weights for every gene for every sample model uncertainty
design <- model.matrix(~celltype, cf)
v <- voom(pd, design)


############
# variance partition model 

# specify mixed effects interacion model
f = ~ age + (1|gender) + (1|sampleid) + (1|celltype) + (1|timepoint) + (1|adjmfc.group) + (1|celltype:timepoint) 

# run model on each gene extract varinace explained 
vp <- fitExtractVarPartModel(exprObj = v, formula = f, data = cf, REML = FALSE, BPPARAM = pparam)
saveRDS(vp, file = paste0(datapath, 'vp.rds'))
vp  = readRDS(file = here('mid_res/variance_partition/generated_data/vp.rds'))

# plot 
p = plotVarPart(vp)
dat$variable %>% str
dat = p$data
levels(dat$variable) = list(`response group` = 'adjmfc.group', `cell type`  = 'celltype',
                            `celltype:timepoint` = 'celltype:timepoint', sex = 'gender', 
                            subjectID = 'sampleid', timepoint = 'timepoint', age = 'age', 
                            residuals = 'Residuals')

dat = dat %>%  filter(!variable == 'Residuals')
#dat$variable[dat$variable == 'gender'] = 'sex'
p = ggplot(dat, aes(x = reorder(variable, value), y = value, fill = variable , color = variable)) + 
  theme_bw() + 
  theme(axis.text = element_text(color = 'black')) + 
  geom_boxplot(outlier.color = 'red', outlier.alpha = 0.2, outlier.shape = 21, show.legend = FALSE) + 
  ggsci::scale_fill_npg(alpha = 0.5) +
  ggsci::scale_color_npg() +
  theme(axis.text = element_text(color = 'black')) + 
  ylab('% variance explained') + xlab('') +
  coord_flip()
p
ggsave(p, filename = paste0(figpath,'fullvp.pdf'), width = 4.6, height = 1.5)
  

# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scglmmr_0.1.0            forcats_0.5.1            stringr_1.4.0            dplyr_1.0.4              purrr_0.3.4              readr_1.4.0             
# [7] tidyr_1.1.2              tibble_3.0.6             tidyverse_1.3.0          variancePartition_1.20.0 Biobase_2.50.0           BiocGenerics_0.36.1     
# [13] scales_1.1.1             BiocParallel_1.24.1      limma_3.46.0             ggplot2_3.3.3            SeuratObject_4.0.0       Seurat_4.0.1            
# [19] here_1.0.1              
# 
# loaded via a namespace (and not attached):
#   [1] estimability_1.3            scattermore_0.7             coda_0.19-4                 bit64_4.0.5                 multcomp_1.4-16            
# [6] irlba_2.3.3                 DelayedArray_0.16.3         data.table_1.14.0           rpart_4.1-15                RCurl_1.98-1.3             
# [11] doParallel_1.0.16           generics_0.1.0              snow_0.4-3                  TH.data_1.0-10              callr_3.7.0                
# [16] cowplot_1.1.1               usethis_2.0.1               RSQLite_2.2.7               shadowtext_0.0.9            RANN_2.6.1                 
# [21] future_1.21.0               bit_4.0.4                   enrichplot_1.10.2           spatstat.data_2.1-0         xml2_1.3.2                 
# [26] lubridate_1.7.9.2           httpuv_1.5.5                ggsci_2.9                   SummarizedExperiment_1.20.0 assertthat_0.2.1           
# [31] viridis_0.5.1               hms_1.0.0                   promises_1.2.0.1            fansi_0.4.2                 progress_1.2.2             
# [36] caTools_1.18.1              dbplyr_2.1.0                readxl_1.3.1                igraph_1.2.6                DBI_1.1.1                  
# [41] htmlwidgets_1.5.3           spatstat.geom_2.0-1         stats4_4.0.5                ellipsis_0.3.1              ggpubr_0.4.0               
# [46] backports_1.2.1             annotate_1.68.0             deldir_0.2-10               MatrixGenerics_1.2.1        vctrs_0.3.6                
# [51] remotes_2.4.0               ROCR_1.0-11                 abind_1.4-5                 cachem_1.0.4                withr_2.4.1                
# [56] ggforce_0.3.3               emmeans_1.5.4               sctransform_0.3.2           prettyunits_1.1.1           goftest_1.2-2              
# [61] cluster_2.1.2               DOSE_3.16.0                 lazyeval_0.2.2              crayon_1.4.1                labeling_0.4.2             
# [66] edgeR_3.32.1                pkgconfig_2.0.3             tweenr_1.0.2                GenomeInfoDb_1.26.7         nlme_3.1-152               
# [71] pkgload_1.2.1               devtools_2.4.2              rlang_0.4.10                globals_0.14.0              lifecycle_1.0.0            
# [76] miniUI_0.1.1.1              sandwich_3.0-0              downloader_0.4              modelr_0.1.8                cellranger_1.1.0           
# [81] rprojroot_2.0.2             polyclip_1.10-0             GSVA_1.38.2                 matrixStats_0.58.0          lmtest_0.9-38              
# [86] graph_1.68.0                Matrix_1.3-2                carData_3.0-4               boot_1.3-27                 zoo_1.8-8                  
# [91] reprex_1.0.0                pheatmap_1.0.12             ggridges_0.5.3              processx_3.5.2              png_0.1-7                  
# [96] viridisLite_0.3.0           bitops_1.0-6                KernSmooth_2.23-18          blob_1.2.1                  qvalue_2.22.0              
# [101] parallelly_1.23.0           rstatix_0.7.0               ggsignif_0.6.0              S4Vectors_0.28.1            memoise_2.0.0              
# [106] GSEABase_1.52.1             magrittr_2.0.1              plyr_1.8.6                  ica_1.0-2                   gplots_3.1.1               
# [111] zlibbioc_1.36.0             compiler_4.0.5              scatterpie_0.1.7            RColorBrewer_1.1-2          lme4_1.1-26                
# [116] fitdistrplus_1.1-3          cli_2.5.0                   XVector_0.30.0              listenv_0.8.0               patchwork_1.1.1            
# [121] pbapply_1.4-3               ps_1.5.0                    MASS_7.3-53.1               mgcv_1.8-34                 tidyselect_1.1.0           
# [126] stringi_1.5.3               GOSemSim_2.16.1             locfit_1.5-9.4              ggrepel_0.9.1               grid_4.0.5                 
# [131] fastmatch_1.1-0             tools_4.0.5                 rio_0.5.16                  future.apply_1.7.0          rstudioapi_0.13            
# [136] foreign_0.8-81              foreach_1.5.1               gridExtra_2.3               farver_2.0.3                Rtsne_0.15                 
# [141] ggraph_2.0.5                digest_0.6.27               rvcheck_0.1.8               BiocManager_1.30.10         shiny_1.6.0                
# [146] Rcpp_1.0.6                  car_3.0-10                  GenomicRanges_1.42.0        broom_0.7.5                 egg_0.4.5                  
# [151] later_1.1.0.1               RcppAnnoy_0.0.18            org.Hs.eg.db_3.12.0         httr_1.4.2                  AnnotationDbi_1.52.0       
# [156] colorspace_2.0-0            rvest_0.3.6                 XML_3.99-0.6                fs_1.5.0                    tensor_1.5                 
# [161] reticulate_1.18             IRanges_2.24.1              splines_4.0.5               uwot_0.1.10                 statmod_1.4.35             
# [166] spatstat.utils_2.1-0        graphlayouts_0.7.2          plotly_4.9.3                sessioninfo_1.1.1           xtable_1.8-4               
# [171] jsonlite_1.7.2              nloptr_1.2.2.2              tidygraph_1.2.0             ggfun_0.0.4                 testthat_3.0.2             
# [176] R6_2.5.0                    pillar_1.4.7                htmltools_0.5.1.1           mime_0.10                   glue_1.4.2                 
# [181] fastmap_1.1.0               minqa_1.2.4                 clusterProfiler_3.18.1      codetools_0.2-18            fgsea_1.16.0               
# [186] utf8_1.1.4                  pkgbuild_1.2.0              mvtnorm_1.1-1               lattice_0.20-41             spatstat.sparse_2.0-0      
# [191] pbkrtest_0.5-0.1            curl_4.3                    leiden_0.3.7                colorRamps_2.3              gtools_3.8.2               
# [196] zip_2.1.1                   openxlsx_4.2.3              GO.db_3.12.1                survival_3.2-10             desc_1.3.0                 
# [201] munsell_0.5.0               DO.db_2.9                   GenomeInfoDbData_1.2.4      iterators_1.0.13            haven_2.3.1                
# [206] reshape2_1.4.4              gtable_0.3.0                spatstat.core_2.0-0  
# 
# 
# 
