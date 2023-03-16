# b cell figures 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
suppressMessages(library(magrittr))
suppressMessages(library(emmeans))
suppressMessages(library(Seurat))
source(here('functions/MattPMutils.r'))
figpath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/figuresV4/bsig/")
dir.create(figpath)
datapath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/bsig_data/")
dir.create(datapath)


# B cell signals from CITE-seq cohort 
gc = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/gc.rds'))
d = gc %>% bind_rows(.id = 'celltype')  %>% 
  filter(celltype == 'BC_Naive') %>% 
  filter(padj < 0.1) 

###### gsea plot subset
mtheme1 = list(
  theme_bw(base_size = 10.5), 
  theme(text = element_text(color = 'black')),
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold",family = "Helvetica"), 
        axis.text.y =  element_text(size = 12, color = 'black'))
)
p = ggplot(d, aes(x = NES, y = reorder(pathway, NES),  
                fill = celltype, size = -log10(padj)), group = celltype ) + 
  mtheme1 +
  theme(axis.text.y  = element_text(size = 9))  + 
  ylab("") +
  xlab('Normalized Enrichment Score') + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  geom_point(shape = 21 , fill = 'deepskyblue3') 
ggsave(p,filename = paste0(figpath, 'BCNaive.as03.enrichment.pdf'), width = 6, height = 2)


# Load day 1 object for both cohorts bcells   
s = readRDS(file = "data/h1h5_annotated_with_meta.rds")
md = s@meta.data %>% 
  filter(celltype_joint == 'BC_Naive') %>% 
  filter(time_cohort == 'd1')
umi = s@raw.data[ ,md$barcode_check]
adt = s@assay$CITE@data[ ,md$barcode_check]

# log normalize rna 
s = CreateSeuratObject(counts = umi, meta.data = md)
s = NormalizeData(s,normalization.method = 'LogNormalize')

# plot B cell protein distributions 
d = cbind(s@meta.data, as.data.frame(t(adt)))
prot_vis= c("CD19_PROT",  "CD20_PROT", "IgD_PROT",  "CD27_PROT","IgM_PROT", 
            "CD21_PROT", "CD40_PROT", "CD38_PROT", "CD24_PROT", "CD14_PROT", 
            "CD3_PROT")
dpl = d %>% 
  filter(celltype_joint == "BC_Naive") %>% 
  select(all_of(prot_vis), sample, cohort) %>% 
  gather(protein, dsb_norm_value, prot_vis[1]:prot_vis[length(prot_vis)])
dpl$protein = factor(dpl$protein, levels = rev(prot_vis))
dpl$protein = str_sub(dpl$protein, 1, -6)
dpl$cohort[dpl$cohort == 'H5N1'] = 'AS03'
dpl$cohort[dpl$cohort == 'H1N1'] = 'No AS03'
p = ggplot(dpl, aes(x = dsb_norm_value, y = reorder(protein, dsb_norm_value), color = cohort, fill = cohort )) + 
  ggridges::geom_density_ridges2(show.legend = FALSE, size = 0.3 ) +
  theme_bw() +
  facet_wrap(~cohort) + 
  geom_vline(xintercept = 0, color = 'black', linetype  = 'dashed') + 
  scale_fill_manual(values = c(col.alpha("grey", 0.8), col.alpha("deepskyblue3", 0.8))) + 
  scale_color_manual(values =c(col.alpha("grey", 0.8), col.alpha("deepskyblue3", 0.8))) + 
  ggtitle("Naive B cell cluster") + 
  theme(strip.background = element_blank()) + 
  theme(axis.text.y = element_text(color = "black")) + 
  ylab("") + xlab("dsb normalized protein")  
p
ggsave(p, filename = paste0(figpath, "BCN_cohort_proteindistributions.pdf"), width = 3, height = 3.8)


# B cell state signature analysis 
# extract signature geens  
gsea1 = readRDS(here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/gc.rds'))
mods = c("CD40_ACT", "REACTOME_ACTIVATION_OF_BH3_ONLY_PROTEINS", 
         "KEGG_P53_SIGNALING_PATHWAY", "LI.M160 leukocyte differentiation")
cd40 = readRDS('signature_curation/combined_sig_sub.rds')['CD40_ACT']

# derive apoptosis signature
gsea1$BC_Naive %>% 
  filter(pathway %in% mods) %$% 
  leadingEdge
apoptosis.signature =
  list('apoptosis.signature' = 
         gsea1$BC_Naive %>%
         filter(pathway %in% mods[2:4]) %$% leadingEdge %>%
         unlist(use.names = FALSE) %>%  
         unique())
sig.test = c(cd40, apoptosis.signature)
saveRDS(sig.test,file = paste0(datapath,'sig.test.rds'))


##################
# fit single cell model
# score modules 
ms = WeightedCellModuleScore(gene_matrix = s@assays$RNA@data, 
                             module_list = sig.test, 
                             cellwise_scaling = FALSE,
                             return_weighted = FALSE )
# combine score and metadata 
d = cbind(s@meta.data, ms)
index1 = names(sig.test)[1]; 
index2 = names(sig.test)[length(sig.test)]

# Calculate d1 FC of average module expression 
ddf = d %>% 
  group_by(sample, sampleid, cohort, timepoint,  celltype_joint) %>% 
  summarise_at(.vars = names(sig.test), .funs = mean) %>% 
  ungroup() %>% 
  gather(module, average, index1:index2) %>% 
  mutate(celltype_module = paste(celltype_joint, module, sep = "~")) %>% 
  arrange(celltype_joint, sampleid) %>% 
  mutate(fold_change = lead(average) - average) 


scale.simple = function(x){ (x - mean(x))/ sd(x)}
signal_cor = 
  ddf %>% 
  filter(timepoint == "d0") %>% 
  filter(module %in% c( 'CD40_ACT', 'apoptosis.signature')) %>% 
  select(sample, cohort,  module, fold_change) %>% 
  spread(module, fold_change) 
signal_cor$apoptosis.signature = scale.simple(signal_cor$apoptosis.signature)
signal_cor$CD40_ACT = scale.simple(signal_cor$CD40_ACT)

p = 
  ggplot(signal_cor %>% mutate(timepoint = str_sub(sample, -2, -1)), 
         aes(x = apoptosis.signature, y = CD40_ACT)) + 
  theme_bw() +  
  geom_smooth(method = "lm", color = col.alpha('black', 0.8))  + 
  xlab('B cell apoptosis signature fold change') + 
  ylab('CD40 Activation signature fold change') + 
  geom_point(aes(fill = cohort), size = 3, shape = 21, show.legend = FALSE) + 
  scale_fill_manual(values = c(col.alpha("grey", 0.8), col.alpha("deepskyblue3",0.8))) + 
  ggpubr::stat_cor(method = "pearson", label.x.npc = 0.01, label.y.npc = 0.01) + 
  ggtitle("Naive B cells")
ggsave(p, filename = paste0(figpath, "CD40score_vs_apoptosissig.pdf"), width = 3.2, height = 3.2)  
saveRDS(signal_cor, file = paste0(datapath, 'signalcor.rds'))


# Fit mixed model to apoptosis signature. 
d$cohort_timepoint = factor(d$cohort_timepoint, levels = c("H1N1_d0", "H1N1_d1", "H5N1_d0", "H5N1_d1"))
d$sex = factor(d$gender)
c00 = c(1,0,0,0); 
c01 = c(0,1,0,0); 
c10 = c(0,0,1,0); 
c11 = c(0,0,0,1) 
contrast_2 = list("time1vs0_group2vs1" = ((c11 - c10) - (c01 - c00)), "time0_group2vs1" = (c10 - c00))
f1 = 'apoptosis.signature ~ 0 + cohort_timepoint + age + sex + (1|sampleid)'
m1 = lme4::lmer(formula = f1, data = d)
emm1 = emmeans(object = m1, specs = ~ cohort_timepoint, data = d, lmer.df = "asymptotic")
contrast_fit = emmeans::contrast(emm1, method = contrast_2)
msummary1 = summary(contrast_fit,infer = c(TRUE, TRUE))
msummary1$module = 'apoptosis.signature'
# contrast           estimate      SE  df asymp.LCL asymp.UCL z.ratio p.value module             
# time1vs0_group2vs1  -0.1168 0.00655 Inf  -0.12961   -0.1039 -17.835 <.0001  apoptosis.signature
# time0_group2vs1      0.0172 0.01293 Inf  -0.00814    0.0425   1.331 0.1833  apoptosis.signature
# 
# Results are averaged over the levels of: sex 
# Degrees-of-freedom method: asymptotic 
# Confidence level used: 0.95
saveRDS(msummary1, file = paste0(datapath,"apoptosis_signature_singlecellmodel_result.rds"))


# visualize 
# plotsingle cell distributionn and emmeans contrasts 
cu = c("grey48", "grey", "grey48",  "deepskyblue3")
cu.alpha = sapply(cu, col.alpha, alpha = 0.8) %>% unname()

# set theme 
plot.aes = list(theme_bw(), 
              theme(axis.title.x = element_text(size = 15),
                    axis.title.y = element_text(size = 15)), 
              scale_color_manual('grey'))

em_aes = list(theme_bw(), 
              coord_flip(), 
              theme(axis.title.y = element_text(size = 7),
                    axis.text.y = element_text(size = 7)), 
              scale_color_manual('grey'))

# combined signature change emm in p1 and change y value in p0
p0 = ggplot(d, aes(x = cohort_timepoint, y = apoptosis.signature, fill = cohort_timepoint )) + 
  geom_violin(show.legend = F,trim = TRUE) + 
  plot.aes + 
  ylab('apoptosis signature') + 
  xlab('vaccine group ~ time') + 
  scale_fill_manual(values = cu.alpha) +
  ggtitle('Naive B cells') +
  theme(axis.title.x = element_text(size = 12))
ggsave(p0, filename = paste0(figpath, 'apoptosis.sig.cells.pdf'), width = 4, height = 3.5)
p1 = plot(emm1) +
  em_aes + 
  theme(axis.text.x = element_blank())
ggsave(p1, filename = paste0(figpath, 'apoptosis.sig.emmeans.pdf'), width = 1.2, height =3 )
p2 = plot(msummary1) + 
  theme_bw() +
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  ggtitle(unique(msummary1$module))
ggsave(p2, filename = paste0(figpath, 'contrast.emmeans.pdf'), width = 4, height = 1.2)


sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] emmeans_1.5.4      SeuratObject_4.0.0 Seurat_4.0.1       magrittr_2.0.1     scglmmr_0.1.0      here_1.0.1         forcats_0.5.1     
# [8] stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0        tidyr_1.1.2        tibble_3.0.6       ggplot2_3.3.3     
# [15] tidyverse_1.3.0   
# 
# loaded via a namespace (and not attached):
# [1] estimability_1.3            scattermore_0.7             coda_0.19-4                 knitr_1.39                  bit64_4.0.5                
# [6] FSA_0.9.0                   irlba_2.3.3                 multcomp_1.4-16             DelayedArray_0.16.3         rpart_4.1-15               
# [11] data.table_1.14.0           RCurl_1.98-1.3              doParallel_1.0.16           generics_0.1.2              BiocGenerics_0.36.1        
# [16] RhpcBLASctl_0.21-247.1      cowplot_1.1.1               TH.data_1.0-10              RSQLite_2.2.7               shadowtext_0.0.9           
# [21] RANN_2.6.1                  future_1.21.0               bit_4.0.4                   enrichplot_1.10.2           spatstat.data_2.1-0        
# [26] httpuv_1.5.5                xml2_1.3.2                  lubridate_1.7.9.2           SummarizedExperiment_1.20.0 assertthat_0.2.1           
# [31] viridis_0.5.1               xfun_0.30                   hms_1.0.0                   promises_1.2.0.1            fansi_0.4.2                
# [36] progress_1.2.2              caTools_1.18.1              dbplyr_2.1.0                readxl_1.3.1                htmlwidgets_1.5.3          
# [41] igraph_1.2.6                DBI_1.1.1                   spatstat.geom_2.4-0         stats4_4.0.5                ellipsis_0.3.2             
# [46] ggpubr_0.4.0                backports_1.2.1             annotate_1.68.0             aod_1.3.1                   deldir_1.0-6               
# [51] MatrixGenerics_1.2.1        vctrs_0.4.1                 Biobase_2.50.0              ROCR_1.0-11                 abind_1.4-5                
# [56] cachem_1.0.4                withr_2.4.3                 ggforce_0.3.3               packrat_0.7.0               checkmate_2.0.0            
# [61] sctransform_0.3.2           prettyunits_1.1.1           goftest_1.2-2               cluster_2.1.2               DOSE_3.16.0                
# [66] lazyeval_0.2.2              crayon_1.4.1                labeling_0.4.2              edgeR_3.32.1                pkgconfig_2.0.3            
# [71] tweenr_1.0.2                GenomeInfoDb_1.26.7         nlme_3.1-152                nnet_7.3-15                 rlang_1.0.2                
# [76] globals_0.14.0              lifecycle_1.0.0             miniUI_0.1.1.1              sandwich_3.0-0              downloader_0.4             
# [81] modelr_0.1.8                cellranger_1.1.0            rprojroot_2.0.2             polyclip_1.10-0             GSVA_1.38.2                
# [86] matrixStats_0.58.0          lmtest_0.9-38               graph_1.68.0                Matrix_1.3-2                carData_3.0-4              
# [91] boot_1.3-27                 zoo_1.8-8                   base64enc_0.1-3             reprex_1.0.0                ggridges_0.5.3             
# [96] pheatmap_1.0.12             png_0.1-7                   viridisLite_0.3.0           bitops_1.0-6                KernSmooth_2.23-18         
# [101] blob_1.2.1                  qvalue_2.22.0               parallelly_1.23.0           jpeg_0.1-8.1                rstatix_0.7.0              
# [106] S4Vectors_0.28.1            ggsignif_0.6.0              scales_1.1.1                memoise_2.0.0               GSEABase_1.52.1            
# [111] plyr_1.8.6                  ica_1.0-2                   gplots_3.1.1                zlibbioc_1.36.0             compiler_4.0.5             
# [116] scatterpie_0.1.7            RColorBrewer_1.1-2          lme4_1.1-26                 fitdistrplus_1.1-3          cli_3.3.0                  
# [121] XVector_0.30.0              listenv_0.8.0               pbapply_1.4-3               patchwork_1.1.1             htmlTable_2.1.0            
# [126] Formula_1.2-4               mgcv_1.8-34                 MASS_7.3-53.1               tidyselect_1.1.0            stringi_1.5.3              
# [131] GOSemSim_2.16.1             locfit_1.5-9.4              latticeExtra_0.6-29         ggrepel_0.9.1               GeneOverlap_1.26.0         
# [136] grid_4.0.5                  fastmatch_1.1-0             tools_4.0.5                 future.apply_1.7.0          parallel_4.0.5             
# [141] rio_0.5.16                  rstudioapi_0.13             foreach_1.5.1               foreign_0.8-81              gridExtra_2.3              
# [146] farver_2.0.3                Rtsne_0.15                  ggraph_2.0.5                digest_0.6.27               rvcheck_0.1.8              
# [151] BiocManager_1.30.10         shiny_1.6.0                 Rcpp_1.0.6                  GenomicRanges_1.42.0        car_3.0-10                 
# [156] broom_0.7.5                 egg_0.4.5                   later_1.1.0.1               RcppAnnoy_0.0.18            org.Hs.eg.db_3.12.0        
# [161] httr_1.4.2                  AnnotationDbi_1.52.0        Rdpack_2.1.1                colorspace_2.0-0            tensor_1.5                 
# [166] rvest_0.3.6                 XML_3.99-0.6                fs_1.5.0                    reticulate_1.18             IRanges_2.24.1             
# [171] splines_4.0.5               uwot_0.1.10                 statmod_1.4.35              spatstat.utils_2.3-0        graphlayouts_0.7.2         
# [176] plotly_4.9.3                xtable_1.8-4                jsonlite_1.7.2              nloptr_1.2.2.2              tidygraph_1.2.0            
# [181] ggfun_0.0.4                 R6_2.5.0                    Hmisc_4.5-0                 mime_0.10                   htmltools_0.5.2            
# [186] pillar_1.4.7                glue_1.6.2                  fastmap_1.1.0               minqa_1.2.4                 clusterProfiler_3.18.1     
# [191] BiocParallel_1.24.1         codetools_0.2-18            fgsea_1.16.0                utf8_1.1.4                  mvtnorm_1.1-1              
# [196] spatstat.sparse_2.0-0       lattice_0.20-41             pbkrtest_0.5-0.1            slanter_0.2-0               curl_4.3                   
# [201] leiden_0.3.7                gtools_3.8.2                zip_2.1.1                   GO.db_3.12.1                openxlsx_4.2.3             
# [206] survival_3.2-10             limma_3.46.0                munsell_0.5.0               DO.db_2.9                   GenomeInfoDbData_1.2.4     
# [211] iterators_1.0.13            variancePartition_1.25.6    haven_2.3.1                 reshape2_1.4.4              gtable_0.3.0               
# [216] spatstat.core_2.0-0         rbibutils_2.0  

