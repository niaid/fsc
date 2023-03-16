Multiscale integration of human and single-cell immune variations reveals unadjuvanted
vaccine high responders are naturally adjuvanted
================

-   <a href="#table-of-contents" id="toc-table-of-contents">Table of
    Contents</a>
    -   <a href="#instructions-for-analysis-workflow-"
        id="toc-instructions-for-analysis-workflow-">Instructions for analysis
        workflow. <a name="instructions"></a></a>
    -   <a href="#fig1--fig-s1-sample-cell-frequency-and-protein-distributions-"
        id="toc-fig1--fig-s1-sample-cell-frequency-and-protein-distributions-">Fig1
        &amp; Fig S1. sample, cell frequency, and protein distributions
        <a name="Fig.1.1"></a></a>
    -   <a
        href="#transcriptome-analysis-of-manually-gated-plasmablasts-and-activated-b-cells-"
        id="toc-transcriptome-analysis-of-manually-gated-plasmablasts-and-activated-b-cells-">transcriptome
        analysis of manually gated plasmablasts and activated B cells.
        <a name="fig1.2"></a></a>
    -   <a
        href="#fig-1-multivariate-analysis-of-human-and-cell-type-variations-"
        id="toc-fig-1-multivariate-analysis-of-human-and-cell-type-variations-">Fig
        1. multivariate analysis of human and cell type variations.
        <a name="fig1.3"></a></a>
    -   <a
        href="#fig-2--fig-s2-mixed-effects-timed-vaccination-response-model--unadjuvanted-cohort-"
        id="toc-fig-2--fig-s2-mixed-effects-timed-vaccination-response-model--unadjuvanted-cohort-">Fig
        2 &amp; Fig S2. mixed effects timed vaccination response model –
        unadjuvanted cohort. <a name="fig2.1"></a></a>
    -   <a
        href="#fig-s2-visualization-of-day-7-post-vaccination-phenotypes-and-predictive-signature-deconvolution-"
        id="toc-fig-s2-visualization-of-day-7-post-vaccination-phenotypes-and-predictive-signature-deconvolution-">Fig
        S2 visualization of day 7 post vaccination phenotypes and predictive
        signature deconvolution <a name="fig2.2"></a></a>
    -   <a
        href="#fig-2-bottom-up-single-cell-reconstruction-of-single-cell-monocyte-pseudotime-"
        id="toc-fig-2-bottom-up-single-cell-reconstruction-of-single-cell-monocyte-pseudotime-">Fig
        2. bottom up single cell reconstruction of single cell monocyte
        pseudotime <a name="fig2.3"></a></a>
    -   <a
        href="#fig-3--figs3-mixed-effects-timed-vaccination-response-model--as03-cite-seq-cohort-"
        id="toc-fig-3--figs3-mixed-effects-timed-vaccination-response-model--as03-cite-seq-cohort-">Fig
        3. &amp; FigS3 mixed effects timed vaccination response model – AS03
        CITE-seq cohort <a name="fig3.1"></a></a>
    -   <a href="#fig-3--figs3-b-cell-as03-phenotype-analysis-"
        id="toc-fig-3--figs3-b-cell-as03-phenotype-analysis-">Fig 3. &amp; FigS3
        B cell AS03 phenotype analysis <a name="fig3.2"></a></a>
    -   <a
        href="#fig-3--figs3-mixed-effects-timed-vaccination-response-model--as03-validatiton-cohort-"
        id="toc-fig-3--figs3-mixed-effects-timed-vaccination-response-model--as03-validatiton-cohort-">Fig
        3. &amp; FigS3 mixed effects timed vaccination response model – AS03
        Validatiton cohort <a name="fig3.3"></a></a>
    -   <a
        href="#fig-3--figs3-as03-specific-cell-phenotypes-validation-comparison-figures-"
        id="toc-fig-3--figs3-as03-specific-cell-phenotypes-validation-comparison-figures-">Fig
        3. &amp; FigS3 AS03 specific cell phenotypes validation comparison
        figures <a name="fig3.4"></a></a>
    -   <a
        href="#fig-4-define-high-responder-baseline-cell-phenotypes-from-multivariate-model-with-enrichment-"
        id="toc-fig-4-define-high-responder-baseline-cell-phenotypes-from-multivariate-model-with-enrichment-">Fig
        4. Define high responder baseline cell phenotypes from multivariate
        model with enrichment <a name="fig4.1"></a></a>
    -   <a
        href="#fig-4-correlate-expression-of-baseline-high-responder-phenotypes-with-plasmablast-response-"
        id="toc-fig-4-correlate-expression-of-baseline-high-responder-phenotypes-with-plasmablast-response-">Fig
        4. Correlate expression of baseline high responder phenotypes with
        plasmablast response <a name="fig4.2"></a></a>
    -   <a
        href="#fig-4--figs4-construct-and-visualize-high-responder-multicellular-network-phenotypes"
        id="toc-fig-4--figs4-construct-and-visualize-high-responder-multicellular-network-phenotypes">Fig
        4. &amp; FigS4 Construct and visualize high responder multicellular
        network phenotypes<a name="fig4.3"></a></a>
    -   <a href="#fig-4-early-kinetics-of-baseline-states-"
        id="toc-fig-4-early-kinetics-of-baseline-states-">Fig 4. Early kinetics
        of baseline states <a name="fig4.4"></a></a>
    -   <a
        href="#fig-4-analysis-of-mrna-vaccine-data-to-define-induction-of-high-responder-phenotypes-"
        id="toc-fig-4-analysis-of-mrna-vaccine-data-to-define-induction-of-high-responder-phenotypes-">Fig
        4. Analysis of mRNA vaccine data to define induction of high responder
        phenotypes <a name="fig4.5"></a></a>
    -   <a
        href="#fig5-define-and-test-as03-specific-cell-phenotypes-in-high-responders-at-baseline-"
        id="toc-fig5-define-and-test-as03-specific-cell-phenotypes-in-high-responders-at-baseline-">Fig.5.
        Define and test AS03 specific cell phenotypes in high responders at
        baseline <a name="fig5.1"></a></a>
    -   <a
        href="#fig5-analysis-of-cell-frequency-of-activated-monocyte-phenotypes-in-flow-cytometry-data-"
        id="toc-fig5-analysis-of-cell-frequency-of-activated-monocyte-phenotypes-in-flow-cytometry-data-">Fig.5.
        Analysis of cell frequency of activated monocyte phenotypes in flow
        cytometry data <a name="fig5.2"></a></a>
    -   <a href="#fig5-analysis-of-cytof-stimulation-phenotypes-"
        id="toc-fig5-analysis-of-cytof-stimulation-phenotypes-">Fig.5. Analysis
        of CyTOF stimulation phenotypes <a name="fig5.3"></a></a>
    -   <a href="#write-output-" id="toc-write-output-">Write output
        <a name="output"></a></a>
    -   <a href="#low-level-bioinformatic-processing-to-generate-starting-data-"
        id="toc-low-level-bioinformatic-processing-to-generate-starting-data-">Low-level
        bioinformatic processing to generate starting data
        <a name="preprocessing"></a></a>

**Code to reproduce all manuscript results and figures**  
<br> Matt Mulè <br>

The data associated with this manuscript for use with the analysis code
here is found here doi: 10.5281/zenodo.7365959.  
**The zenodo repository also includes CITE-seq data in additional
formats for reuse with other packages.**  
`flu_vacc_CITEseq_Seurat4.rds` - a Seurat version 4 object (separate
assays for RNA and ADT).  
`flu_vacc_CITEseq_combined_assay_SingleCellExperiment.rds` - a
[SingleCellExperiment](https://www.bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)
object (RNA + ADT combined in single matrix).  
`flu_vacc_CITEseq_combinedassay.h5ad` - an Anndata object for analysis
in python / [scanpy](https://scanpy.readthedocs.io/en/latest/). (RNA +
ADT combined in single matrix).

## Table of Contents

1.  [Instructions for analysis workflow](#instructions)  
2.  [Fig1 & Fig S1. sample, cell frequency, and protein
    distributions](#fig1.1)  
3.  [Fig S1. transcriptome analysis of manually gated plasmablasts and
    activated B cells](#fig1.2)  
4.  [Fig 1. multivariate analysis of human and cell type
    variations](#fig1.3)  
5.  [Fig 2 & Fig S2. mixed effects timed vaccination response model –
    unadjuvanted cohort](#fig2.1)  
6.  [Fig S2 visualization of day 7 post vaccination phenotypes and
    predictive signature deconvolution](#fig2.2)
7.  [Fig 2. bottom up single cell reconstruction of single cell monocyte
    pseudotime](#fig2.3)  
8.  [Fig 3. & FigS3 mixed effects timed vaccination response model –
    AS03 CITE-seq cohort](#fig3.1)  
9.  [Fig 3. & FigS3 B cell AS03 phenotype analysis](#fig3.2)  
10. [Fig 3. & FigS3 mixed effects timed vaccination response model –
    AS03 Validatiton cohort](#fig3.3)  
11. [Fig 3. & FigS3 AS03 specific cell phenotypes figure
    generation](#fig3.4)  
12. [Fig 4. Define high responder baseline cell phenotypes from
    multivariate model with enrichment](#fig4.1)  
13. [Fig 4. Correlate expression of baseline high responder phenotypes
    with plasmablast response](#fig4.2)  
14. [Fig 4. & FigS4 Construct and visualize high responder multicellular
    network phenotypes](#fig4.3)  
15. [Fig 4. Early kinetics of baseline states](#fig4.4)  
16. [Fig 4. Analysis of mRNA vaccine data to define induction of high
    responder phenotypes](#fig4.5)
17. [Fig.5. Define and test AS03 specific cell phenotypes in high
    responders at baseline](#fig5.1)
18. [Fig.5. Analysis of cell frequency of activated monocyte phenotypes
    in flow cytometry data](#fig5.2)
19. [Fig.5. Analysis of CyTOF stimulation phenotypes](#fig5.3)
20. [Write output](#output)
21. [Low-level bioinformatic processing to generate starting
    data](#preprocessing)

### Instructions for analysis workflow. <a name="instructions"></a>

Many scripts in this workflow require the
[scglmmr](https://github.com/MattPM/scglmmr) package associated with
this manuscript. The functions for this package can also be found in
`functions/scglmmr_functions/`.

This repository documents all of the code needed to reproduce the
analysis results of the manuscript. The data for this analysis in the
Zenodo repository below should be added to a data/ directory in the
project root folder. The analysis can now run as is without
modification–Each R script is self-contained, reading data from the
/data folder and writing to figures or results files within each
analysis subdirectory relative to the root directory using the R package
`here`. Unless otherwise noted, scripts were run with R 4.0.5. Package
versions are listed in the table.

**Starting raw data for analysis: 10.5281/zenodo.7365959**

non CITEseq data used in the analysis is saved in the following
directories:  
data  
–vand  
–CHI_H1N1_data  
–full_metadata  
–CHI_H5N1_data  
–GSE171964  
–stim

Annotated CITE-seq data in a Seurat object is saved in the project
directory as h1h5_annotated_with_meta.rds  
\# a version with added sample metadata on each cell is created as below
and saved in the first script: data –h1h5_annotated_with_meta.rds

### Fig1 & Fig S1. sample, cell frequency, and protein distributions <a name="Fig.1.1"></a>

*This section is run with R 3.5.1*  
Calculate and visualize distribution of individuals across protein based
clusters.  
Biaxial plots of unsupervised cluster dsb normalized protein
expression.  
mid_res/sample_and_protein_distributions/1.sample.lineage.proteinbiaxial.r

``` r
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
```

Calculate and visualize distribution of immune cell frequency across
individuals.

mid_res/sample_and_protein_distributions/2.cell.frequency.map.r

``` r
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
```

Clustered protein histogram distribution across cell types

mid_res/histogram_hclust/hclust_histogram_protein.r

``` r
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
```

Aggregated protein expression heatmap across 780 libraries by sample x
timepoint

mid_res/histogram_hclust/hclust_histogram_protein.r

``` r
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
```

### transcriptome analysis of manually gated plasmablasts and activated B cells. <a name="fig1.2"></a>

Analyze transcriptome of manually gated cells using cell type specific
gene signatures  
mid_res/pblast_abc_integration/1_pbasc_analysis_gates_modscoresv2.r

``` r
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
```

### Fig 1. multivariate analysis of human and cell type variations. <a name="fig1.3"></a>

Analyze variance explained by each factor of a multivariate model fit
across all cell types. This model includes call type as a variable.
mid_res/variance_partition/1_variance_partition_allsamples.r

``` r
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
```

Analyze variance explained by each factor of multivariate models fit
within each cell protein based type.
mid_res/variance_partition/2_variance_partition_withincelltype.r

``` r
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
tab = scglmmr::SubjectCelltypeTable(metadata = meta, celltype_column = "celltype_joint", sample_column = "sample")

# remove cells prior to pseudobulk analysis 
meta = meta[!meta$celltype_joint %in% c(tab$celltypes_remove, 'DOUBLET'), ]
umi = umi[ ,rownames(meta)]

# creat sample metadata 
samplemd = 
  meta %>% 
  group_by(sample) %>% 
  select(age, gender, sampleid, batch, time_cohort, timepoint, adjmfc.group) %>% 
  summarise_each(funs = unique) %>% 
  ungroup() %>% 
  as.data.frame()


# read pb data from vpar1 script 1 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))

# run variance partition workflow for each cell subset separately
for (i in 1:length(pb)) {
  
  csd = 
    data.frame(sid = colnames(pb[[i]])) %>% 
    mutate(sample_celltype = sid) %>% 
    separate(sid, into = c('sample', 'celltype'), sep = '~') 
  
  # sample_celltype metadata for design matrices 
  cf = full_join(csd, samplemd, by = 'sample') %>%
    column_to_rownames('sample_celltype')
  
  #############
  # process bulk data indexed over each celltype 
  #filter lowly expressed (in this case basically unexpressed genes)
  pd = pb[[i]]
  gene_keep = edgeR::filterByExpr(pd, 
                                  min.count = 5,
                                  min.total.count = 2,
                                  min.prop = 0.5,
                                  group = as.factor(cf$sample))
  print(table(gene_keep))
  pd = pd[gene_keep, ]
  
  # normalize bulk data 
  pd = edgeR::DGEList(counts = pd, samples = cf)
  pd = edgeR::calcNormFactors(object = pd)
  
  ##############
  # Get voom observational weights 
  # these precision weights for every gene for every sample model uncertainty
  design <- model.matrix(~timepoint, cf)
  v <- voom(pd, design, plot = TRUE)
  
  # specify mixed effects interacion model
  f = ~ age + (1|gender)  + (1|sampleid) + (1|timepoint) + (1|adjmfc.group) + (1|timepoint:adjmfc.group) 
  
  # run model on each gene extract varinace explained 
  vp <- fitExtractVarPartModel(exprObj = v, formula = f, data = cf, REML = FALSE, BPPARAM = pparam)
  saveRDS(vp, file = paste0(datapath, names(pb)[i], 'vp.rds'))
  
  # plot 
  p = plotVarPart(vp) 
  p = ggplot(p$data, aes(x = reorder(variable, value), y = value, fill = variable )) + 
    theme_bw() + 
    theme(axis.text = element_text(color = 'black')) + 
    geom_boxplot(outlier.color = 'red', outlier.alpha = 0.2, outlier.shape = 21, show.legend = FALSE) + 
    ggsci::scale_fill_jama() + 
    ylab('% variance explained') + xlab('') +
    coord_flip()
  ggsave(p,filename = paste0(figpath,names(pb)[i],'vp.pdf'), width = 4.5, height = 1.4)
}




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
```

Figure generation from analysis above
mid_res/variance_partition/3_figures_variance_partition_withincelltype.R

``` r
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(variancePartition))
source('functions/MattPMutils.r')
suppressMessages(library(magrittr))
# figpath
figpath = here('mid_res/variance_partition/figures_vars/'); dir.create(figpath, recursive = TRUE)
datapath = here('mid_res/variance_partition/generated_data2/'); dir.create(datapath, recursive = TRUE) 

# parallel opts
# register(SnowParam(4))
# pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)
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
ccu = structure(col[[1]]) 
names(ccu) = str_replace_all(string = names(ccu), pattern = '_', replacement = ' ')
ccu2 = sapply(ccu, col.alpha, 0.5)


# load bulk data 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))
dl = list.files(path = here('mid_res/variance_partition/generated_data/'),
                pattern = '.rds', recursive = TRUE,full.names = TRUE)
dl = dl[-c(15,17)] # remove total bulk 

# get cell type names (file names)
cts = list.files(path = here('mid_res/variance_partition/generated_data/'),
                 pattern = '.rds', recursive = TRUE,full.names = FALSE)
cts = cts[-c(15,17)]
cts = str_replace_all(string = cts,pattern = 'vp.rds', replacement = '')


# read and format variance partition results 
vl = lapply(dl, readRDS)
names(vl) = cts
dl = list()
for (i in 1:length(vl)) {
  p = plotVarPart(vl[[i]])
  p$data$celltype = names(vl[i])
  d = p$data 
  dl[[i]] = d
}

# combine results across cell types 
dl[[i]] %>% head 
test = reduce(dl, .f = rbind)
test2 = test %>% select(-c(gene))
test2$celltype = factor(test2$celltype, levels = cts)

# rename metadata vars 
levels(test2$variable) =
  list(
    `response group` = 'adjmfc.group',
    sex = 'gender',
    subjectID = 'sampleid',
    timepoint = 'timepoint',
    `timepoint:response` = "timepoint:adjmfc.group",
    age = 'age',
    residuals = 'Residuals'
  )


# visualize full results 
test2$celltype = str_replace_all(test2$celltype, pattern = '_',replacement = ' ')
d = test2 %>% filter(variable %in% c('age', 'sex', 'subjectID', 'timepoint'))
d$variable = factor(d$variable,levels = c('subjectID', 'timepoint', 'age', 'sex'))

# add outlier designation 
d = d %>% 
  group_by(celltype, variable) %>%
  mutate(outlier = value > quantile(value, 0.75) + IQR(value) * 1.5) %>%
  ungroup


p = ggplot(d, aes(x = reorder(celltype,value), y = value , color = celltype, fill= celltype)) +
  facet_wrap(~variable, nrow = 1) +
  theme_bw() + 
  theme(axis.text = element_text(color = 'black', size = 10)) +
  geom_boxplot(varwidth = TRUE, 
               outlier.color = 'red', 
               outlier.alpha = 0.2, 
               outlier.shape = NA, 
               show.legend = FALSE, size = 0.3) + 
  geom_point(data = function(x) dplyr::filter_(x, ~ outlier),
             position = 'jitter', 
             shape = 21, size = 1.1, stroke = 0, alpha = 1/3,
             show.legend = FALSE) + 
  scale_color_manual(values = ccu) +
  scale_fill_manual(values = ccu2) +
  ylab('% variance explained') + xlab('') +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 18, color = 'black'),
        axis.text.y = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 16, color = 'black')) + 
  coord_flip() 
p
ggsave(p, filename = paste0(figpath, 'vpartfull_.pdf'), width = 10, height = 5.1) 
ggsave(p, filename = paste0(figpath, 'vpartfull_.png'), width = 10, height = 5.1) 



######## 
# monnocyte 
m = vl$CD14_Mono %>%
  as.data.frame() %>%
  rownames_to_column('gene')

colnames(m) = c(
  'gene',
  'response.group',
  'sex',
  'subjectID',
  'timepoint',
  'timepoint:response',
  'age',
  'residuals'
)

# rank genes by vars 
ds = m[order(desc(m$sampleid)), ]

# get sampele meta data to add to gene data 
samplemd = readRDS(file = here('data/samplemd.rds'))
mdat = pb$CD14_Mono
mdat = edgeR::cpm(mdat, log = TRUE)

# plot genes 
p = plotPercentBars(vl$CD14_Mono[c('DDX3Y', 'TMEM176B',  'STAT1','PPARGC1',    'TP53RK'), ] ) 
levels(p$data$variable) = c('response.group', 'sex', 'SubjectID', 'timepoint', 'timepoint:response', 'age', 'residuals')
levels(p$data$variable) = c('response.group', 'sex', 'subjectID', 'timepoint', 'timepoint:response', 'age', 'residuals')
p = p + 
  ggsci::scale_fill_jama(alpha = 0.9) +
  theme_bw() + 
  theme(axis.text = element_text(size = 15, color = 'black')) + 
  theme(axis.title = element_text(size = 20, color = 'black'))+ 
  theme(legend.position = 'top') +
  guides(fill=guide_legend(nrow=4, byrow=TRUE)) + 
  theme(legend.text = element_text(size = 18), legend.key.size = unit(0.8,units = 'cm')) + 
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(hjust=1)) 
ggsave(p, filename = paste0(figpath,'genesubpct.pdf'), width = 5.6, height = 5)


# genes 
mgene.highlight = c(
  'PPARGC1B',
  'TMEM176B',
  'LILRA3',
  'TP53RK',
  'PRPF19',
  'HLA-DRB5',
  'GBP2',
  'PSME2',
  'VAMP5',
  'STAT1',
  'CD69',
  'MAP3K8',
  'DDX3Y'
)
# make matrix
d2 = as.data.frame(as.matrix(t(mdat[ mgene.highlight, ])))
d2 = d2 %>%
  rownames_to_column('sid') %>% 
  separate(sid, into = c('sample', 'celltype'), sep = '~') 
d2 = full_join(d2, samplemd, by = 'sample')
d2$sampleid = str_sub(d2$sampleid, -3,-1)

# subject 
p = 
  ggplot(d2, aes(x = reorder(sampleid, TMEM176B), y = TMEM176B, fill = timepoint)) +
  theme_bw() +
  xlab('Subject ID') +
  geom_point(shape = 21, size = 2.5, fill = col.alpha('black', 0.7) , color = 'white') + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 5))  + 
  theme(axis.title = element_text(color = 'black', size = 14)) +
  theme(legend.position = c(0.74, 0.24)) +
  theme(legend.key.size = unit(0.2,units = 'cm')) +
  scale_fill_manual(values = c('black', 'black', 'black'))
p
ggsave(p, filename = paste0(figpath, 'TMEM176B.pdf'), width = 2.5, height = 2.5)


# sex 
p = 
  ggplot(d2, aes(x = gender, y = DDX3Y)) +
  theme_bw() + 
  xlab('Sex') +
  theme(axis.title = element_text(color = 'black', size = 14)) +
  geom_boxplot(show.legend = FALSE, fill = col.alpha(acol = 'black', alpha = 0.3)) 
 ggsave(p, filename = paste0(figpath, 'sex_gene.pdf'), width = 2.5, height = 2.5)

# time 
p = 
  ggplot(d2, aes(x = timepoint, y = STAT1)) +
  theme_bw() + 
  xlab('Time') +
  theme(axis.title = element_text(color = 'black', size = 14)) +
  geom_boxplot(show.legend = FALSE, fill = col.alpha(acol = 'black', alpha = 0.3)) 
ggsave(p, filename = paste0(figpath, 'timegene.pdf'), width = 2.5, height = 2.5)

# Age 
p = 
  ggplot(d2, aes(x = age, y = TP53RK)) +
  theme_bw() + 
  theme(axis.title = element_text(color = 'black', size = 14)) +
  geom_point(shape = 21, size = 2.5, fill = col.alpha('black', 0.7) , color = 'white') + 
  scale_fill_manual(values = c('black')) + 
  geom_smooth(method = 'lm', color= 'black') + 
  ggpubr::stat_cor(label.x.npc = 0.1, label.y.npc = 0.1) 
p
ggsave(p, filename = paste0(figpath, 'agegene2.pdf'), width = 2.5, height = 2.5)

# Age 2 
p = 
  ggplot(d2, aes(x = age, y = PPARGC1B)) +
  theme_bw() + 
  theme(axis.title = element_text(color = 'black', size = 14)) +
  geom_point(shape = 21, size = 2.5, fill = col.alpha('black', 0.7) , color = 'white') + 
  scale_fill_manual(values = c('black')) + 
  geom_smooth(method = 'lm', color= 'black') + 
  ggpubr::stat_cor(label.x.npc = 0.1, label.y.npc = 0.1) 
p
ggsave(p, filename = paste0(figpath, 'agegene.pdf'), width = 2.5, height = 2.5)
```

Gene set enrichment of genes ranked by variance explained by age in CD8
T cell subsets.  
mid_res/variance_partition/4_age_variancefraction_enrichment.r

``` r
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(variancePartition))
source('functions/MattPMutils.r')
library(magrittr)
# figpath
figpath = here('mid_res/variance_partition/figures_vars/'); dir.create(figpath, recursive = TRUE)
datapath = here('mid_res/variance_partition/generated_data2/'); dir.create(datapath, recursive = TRUE) 

# parallel opts
# register(SnowParam(4))
# pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)
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
ccu = structure(col[[1]]) 
names(ccu) = str_replace_all(string = names(ccu), pattern = '_', replacement = ' ')
ccu2 = sapply(ccu, col.alpha, 0.5)

# set theme 
mtheme = list(
  theme_bw(), 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
)

# load bulk data 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))
dl = list.files(path = here('mid_res/variance_partition/generated_data/'),
                pattern = '.rds', recursive = TRUE,full.names = TRUE)
dl = dl[-c(15,17)] # remove total bulk 

# get cell type names (file names)
cts = list.files(path = here('mid_res/variance_partition/generated_data/'),
                 pattern = '.rds', recursive = TRUE,full.names = FALSE)
cts = cts[-c(15,17)]
cts = str_replace_all(string = cts,pattern = 'vp.rds', replacement = '')


# read and format variance partition results 
vl = lapply(dl, readRDS)
names(vl) = cts
dl = list()
for (i in 1:length(vl)) {
  p = plotVarPart(vl[[i]])
  p$data$celltype = names(vl[i])
  d = p$data 
  dl[[i]] = d
}


# rank genes by variance fraction assciated with age 
hlmk = readRDS(file = here('signature_curation/hallmark.rds'))
# parallel opts
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)


agelist = slist = list()
#vars = c('age', 'timepoint', 'gender', 'sampleid')
for (u in 1:length(vl)) {
  
  # get variance fractions
  m = vl[[u]] %>% 
    as.data.frame() %>%  
    rownames_to_column('gene') 
  
  # rank genes by sex 
  da = m[order(desc(m$age)), ] 
  rank.age = structure(da$age, names= da$gene)
  
  # rank genes by subject 
  ds = m[order(desc(m$sampleid)), ] 
  rank.subject = structure(ds$sampleid, names= ds$gene)
  
  # format lists 
  agelist[[u]] = rank.age
  slist[[u]] = rank.subject
}

# run gsea for age and subject 
age.gs = scglmmr::FgseaList(rank.list.celltype = age.gsea,pathways = hlmk, scoreType = "pos", BPPARAM= pparam)

age.gsea = sgsea = list()
for (u in 1:length(agelist)) {
  age.gsea[[u]] = fgsea::fgsea(hlmk, agelist[[u]], scoreType = "pos",  BPPARAM = pparam)
  sgsea[[u]] = fgsea::fgsea(hlmk, slist[[u]], scoreType = "pos",  BPPARAM = pparam)
}
for (u in 1:length(agelist)) {
  age.gsea[[u]]$celltype = names(vl)[u]
}
p = scglmmr::PlotFgsea(gsea_result_list = age.gsea, padj_filter = 0.05 )


age.gsea.sub = lapply(age.gsea, function(x)
  x %>%  filter(
    pathway %in% c(
      'HALLMARK_IL2_STAT5_SIGNALING',
      'HALLMARK_IL6_JAK_STAT3_SIGNALING',
      'HALLMARK_INFLAMMATORY_RESPONSE',
      'HALLMARK_ALLOGRAFT_REJECTION'
    )
  ))

li = scglmmr::LeadingEdgeIndexed(gsea.result.list = age.gsea.sub,padj.threshold = 0.05)
li = Filter(li, f = length)
li$CD8_CD161_Tcell


dage = bind_rows(age.gea.sub, .id = 'celltype')
dage = dage %>%  filter(padj < 0.05)



dage$pathway = str_replace_all(dage$pathway,pattern = 'HALLMARK_',replacement = '')
dage$pathway = str_replace_all(dage$pathway,pattern = '_',replacement = ' ')
cu = c("olivedrab", "#848482")
cu2 = sapply(cu, col.alpha, 0.5)
p = 
  ggplot(dage %>% 
           filter(celltype %in% c('CD8_Naive_Tcell', 'CD8_CD161_Tcell')), 
         aes(x = NES, y = pathway, label = celltype, 
             group = celltype, fill = celltype, color = celltype)) + 
  mtheme + 
  geom_linerange(aes(x = NES, color = celltype, xmin = 0, xmax = NES),
                 position = position_dodge(width = 0.35)) +
  geom_point(aes(x = NES, color = celltype),position = position_dodge(width = 0.35)) +
  geom_point(shape = 21, size = 2.5, position = position_dodge(width = 0.35)) +
  ggsci::scale_fill_npg() + 
  theme(axis.text = element_text(color = 'black')) + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  ylab('') + 
  xlab('Normalized Enrichment Score') + 
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle('Age associated variance enrichment') 

p
ggsave(p,filename = paste0(figpath, 'tcell_age.pdf'), width = 6, height =2.2)




# CD8 Naive 
# logcpm matrix
cd8.genes = unique(unlist(li$CD8_Naive_Tcell, use.names = FALSE))
mdat = pb$CD8_Naive_Tcell
mdat = edgeR::cpm(mdat, log = TRUE)
d2 = as.data.frame(as.matrix(t(mdat[cd8.genes,])))
d2 = d2 %>% 
  rownames_to_column('sid') %>% 
  separate(sid, into = c('sample', 'celltype'), sep = '~')
d2 = full_join(d2, samplemd, by = 'sample')
d2$sampleid = str_sub(d2$sampleid, -3, -1)

# scale age
scale.simple = function(x) {
  (x - mean(x)) / sd(x)
}
d2$age = scale.simple(d2$age)

# fit models
dmat = d2 %>% select(age, all_of(cd8.genes))
age.scaled = d2$age

age.coef = apply(dmat[, 2:ncol(dmat)], MARGIN = 2, function(x) {
  y = lm(x ~ 0 + age.scaled)
  return(y)
})

age.res = lapply(age.coef, broom::tidy)
age.res = bind_rows(age.res,.id = 'gene')

age.pos = age.res %>% 
  filter(estimate > 0) %$% 
  gene
plot(agelist[[11]][age.pos])

age.var = data.frame(age.var = agelist[[11]][age.pos]) %>% 
  rownames_to_column('gene')

p = 
  ggplot(age.var, aes(y = reorder(gene, age.var) , x = age.var*100, label = gene )) +
  mtheme + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab('') + 
  xlab('% variance explained by age') +
  xlim(c(-5, 40)) +
  geom_point(color = "#4DBBD5FF") +
  ggrepel::geom_text_repel(data = age.var %>% 
                             filter(gene %in% c('KLF6', 'RGS16','TNF', 'IL12A', 'HLA-DQA1', 
                                                'ABCA1', 'CD70','GZMB', 'IL10A', "FASLG", 
                                                'CCL5', 'NOD2', 'CD74', 'IFNG', 'CCL4', 'FAS')),
                           size = 3,
                           force        = 1,
                           nudge_x      = -20,
                           direction    = "y",
                           hjust        = 1,
                           segment.size = 0.001) + 
  ggtitle('Positive association with age\nCD8 Naive T cells') + 
  theme(axis.title = element_text(size = 18))

ggsave(p,filename = paste0(figpath, 'cd8n.varexp.age.positive.pdf'), width = 4, height = 4)




# CD8 CD161 
# logcpm matrix
cd161.genes = unique(unlist(li$CD8_CD161_Tcell,use.names = FALSE))
mdat = pb$CD8_CD161_Tcell
mdat = edgeR::cpm(mdat, log = TRUE)
d2 = as.data.frame(as.matrix(t(mdat[cd161.genes, ])))
d2 = d2 %>%
  rownames_to_column('sid') %>%
  separate(sid, into = c('sample', 'celltype'), sep = '~') 
d2 = full_join(d2, samplemd, by = 'sample')
d2$sampleid = str_sub(d2$sampleid, -3,-1)

# scale age 
scale.simple = function(x) { (x - mean(x)) / sd(x)}
d2$age = scale.simple(d2$age)

# fit models 
dmat = d2 %>% select(age, all_of(cd161.genes))
age.scaled = d2$age

age.coef = apply(dmat[ ,2:ncol(dmat)],MARGIN = 2, function(x) { 
  y = lm(x ~ 0 + age.scaled)
  return(y)
} )

age.res = lapply(age.coef, broom::tidy)
age.res = bind_rows(age.res,.id = 'gene')

age.pos = age.res %>%  filter(estimate > 0) %$% gene
names(vl)
age.var = data.frame(age.var = agelist[[9]][age.pos]) %>%  rownames_to_column('gene')
p = 
  ggplot(age.var, aes(y = reorder(gene, age.var) , x = age.var*100, label = gene )) +
  mtheme + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab('') + 
  xlab('% variance explained by age') +
  xlim(c(-5, 40)) +
  geom_point(color = "#E64B35FF") +
  ggrepel::geom_text_repel(data = age.var %>% 
                             filter(age.var > 0.08) %>% 
                             filter(gene %in% c('HLA-DQA1', 'KRT1', 'IFNG', 'CCL4', 'KLRD1', 
                                                "CCL5", "CFP", "CD74", "HLA-DRA",  "IL2RB",  
                                                "IL17RA",  "TNFRSF8", "MAP3K8",   "SERPINB6", "CD38",   
                                                "TYK2", "CD70",     "PROK2" ,   "RGS1"    )),
                           size = 3,
                           force        = 1,
                           nudge_x      = -20,
                           direction    = "y",
                           hjust        = 1,
                           segment.size = 0.001) + 
  ggtitle('Positive association with age\nCD8 CD161+ T cells') + 
  theme(axis.title = element_text(size = 18))
ggsave(p,filename = paste0(figpath, 'cd161T.varexp.age.positive.pdf'), width = 4, height = 4)



# write lists 
plot(age.coef)
pos.age = age.coef[age.coef > 0]
neg.age = age.coef[age.coef < 0]
data.table::fwrite(list(names(pos.age)),file = paste0(datapath,'pos.age.cd8n.txt'))
data.table::fwrite(list(names(neg.age)),file = paste0(datapath,'neg.age.cd8n.txt'))
```

### Fig 2 & Fig S2. mixed effects timed vaccination response model – unadjuvanted cohort. <a name="fig2.1"></a>

Mixed effects covariate adjusted model of vaccination effects across
donors within each protein based subset – unadjuvanted cohort. Day 1 and
day 7 post vaccination effects estimates. Derive the effect size of
vaccination effect to rank genes and run enrichment.

mid_res/1_H1N1_pseudobulk_DE/1_h1_mixed_effect_workflow_V4.r

``` r
# H1N1 differential expression testing and gene set enrichment analysis
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(variancePartition))
suppressMessages(library(scglmmr))
source("functions/analysis_functions.R")

# make output directories 
datapath = here("mid_res/1_H1N1_pseudobulk_DE/dataV4/")
dir.create(datapath, recursive = TRUE)
figpath = here("mid_res/1_H1N1_pseudobulk_DE/figuresV4/")
dir.create(figpath, recursive = TRUE)

# parallel options 
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# read processed pseudobulk data 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))

# subset to unadjuvanted cohort and remove cell type string from sample names 
pb = lapply(pb, function(x) x = x[ ,1:40])
cnames = gsub("~.*","",colnames(pb[[1]]))
pb = lapply(pb, function(x) {
  x %>% 
    as.data.frame() %>% 
    setNames(nm = cnames) %>% 
    as.matrix()
})


# sample metadata 
samplemd = readRDS(file = here('data/samplemd.rds')) %>% filter(! adjmfc.group %in% 'AS03')
samplemd$scaledage = as.vector(scale(samplemd$age))
names(samplemd)[names(samplemd) == 'sampleid'] <- 'subjectid'
names(samplemd)[names(samplemd) == 'adjmfc.group'] <- 'group'
samplemd$group = str_replace_all(string = samplemd$group,pattern = ' ', replacement = '')

# format 
samplemd = samplemd %>% 
  mutate(time.group = paste(timepoint, group,sep = "_")) %>% 
  remove_rownames() %>% 
  column_to_rownames('sample')

# relevel combined factor 
samplemd$time.group = factor(samplemd$time.group, 
                             levels = c('d0_high', 'd1_high', 'd7_high', 
                                        'd0_low', 'd1_low', 'd7_low'))

samplemd$timepoint = factor(samplemd$timepoint, levels = c('d0', 'd1', 'd7'))

# create separate model metadata for the separate cohorts being tested 
d1 = samplemd[samplemd$time_cohort == 'd1', ] %>% droplevels()
d7 = samplemd[samplemd$time_cohort == 'd7', ] %>% droplevels()
d0 = samplemd[samplemd$timepoint == 'd0', ] %>% droplevels()


# subset the bulk lists for each time cohort 
d1d = lapply(pb, function(x){ x = x[ , rownames(d1)]})
d7d = lapply(pb, function(x){ x = x[ , rownames(d7)]})
d0d = lapply(pb, function(x){ x = x[ , rownames(d0)]})


################################
# fit day 1 model 
f1 <- ~ 0 + timepoint + batch + gender + age + (1|subjectid) 

# set up contrast matrix (based on first element of list) 
d = edgeR::DGEList(counts = d1d[[1]], samples = d1)
cmat = getContrast(exprObj = d, formula = f1, data = d1, coefficient = c( 'timepointd1', 'timepointd0'))
plotContrasts(cmat)

# run on each subset 
fit1 = v1 = list()
for (i in 1:length(d1d)) {
  # init data 
  meta = d1 
  form = f1 
  contrast_matrix = cmat
  counts = d1d[[i]]
  
  # dge list 
  d = edgeR::DGEList(counts = counts, samples = meta)
  
  # filter cell type specific lowly expressed genes and calc norm factors 
  gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, design = as.factor(d$samples$timepoint))
  print(names(d1d)[i]);print(table(gtable))
  d = d[gtable, keep.lib.sizes=FALSE]
  d = edgeR::calcNormFactors(object = d)
  
  # get voom observation level weights 
  v = voomWithDreamWeights(counts = d, 
                           formula = form,
                           data = meta, 
                           BPPARAM = pparam, 
                           plot = TRUE, save.plot = TRUE)
  # fit contrast mixed model 
  fitmm = dream(exprObj = v, formula = form, data = meta,
                L = contrast_matrix, useWeights = TRUE,
                BPPARAM = pparam, REML = TRUE)
  
  fitmm = variancePartition::eBayes(fit = fitmm)
  # save results 
  v1[[i]] = v
  fit1[[i]] = fitmm
}
names(v1) = names(fit1) = names(d1d)


################################
# fit day 7 model (uses same formula)
f1 <- ~ 0 + timepoint + batch + gender + age + (1|subjectid) 

# set up contrast matrix (based on first element of list) 
d = edgeR::DGEList(counts = d7d[[1]], samples = d7)
cmat = getContrast(exprObj = d, formula = f1, data = d7, coefficient = c( 'timepointd7', 'timepointd0'))
plotContrasts(cmat)

# run on each subset 
fit7 = v7 = list()
for (i in 1:length(d7d)) {
  # init data 
  meta = d7 
  form = f1 
  contrast_matrix = cmat
  counts = d7d[[i]]
  
  # dge list 
  d = edgeR::DGEList(counts = counts, samples = meta)
  
  # filter cell type specific lowly expressed genes 
  gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, design = as.factor(d$samples$timepoint))
  table(gtable)
  d = d[gtable, keep.lib.sizes=FALSE]
  d = edgeR::calcNormFactors(object = d)
  
  # get voom observation level weights 
  v = voomWithDreamWeights(counts = d, 
                           formula = form, 
                           data = meta, 
                           BPPARAM = pparam, 
                           plot = TRUE, save.plot = TRUE)
  # fit contrast mixed model 
  fitmm = dream(exprObj = v, formula = form, data = meta,
                L = contrast_matrix, useWeights = TRUE, 
                BPPARAM = pparam, REML = TRUE)
  
  fitmm = variancePartition::eBayes(fit = fitmm)
  # save results 
  v7[[i]] = v
  fit7[[i]] = fitmm
}
names(v7) = names(fit7) = names(d7d)

# run baseline model using limma 
# set up fixed effects model to run with limma 
mod0 <- model.matrix(~ 0 + group + batch + gender + age, data = d0)
colnames(mod0) = c("high", "low", 'batch2', "genderM", "age")
c0 = makeContrasts(adjmfc = high - low, levels = colnames(mod0))

fit0 = v0 = cont0 = list()
for (i in 1:length(d0d)) {
  # init data 
  meta = d0
  # form = f1 
  contrast_matrix = c0
  counts = d0d[[i]]
  
  # dge list 
  d = edgeR::DGEList(counts = counts, samples = meta)
  
  # filter cell type specific lowly expressed genes ** Change grouping factor for filter by expression to group
  gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, design = as.factor(d$samples$group))
  table(gtable)
  d = d[gtable, keep.lib.sizes=FALSE]
  d = edgeR::calcNormFactors(object = d)
  
  # get voom observation level weights
  v = voom(counts = d, design = mod0, save.plot = TRUE, plot = TRUE)
  #v = voomWithDreamWeights(counts = d, formula = form, data = meta, BPPARAM = pparam, plot = TRUE, save.plot = TRUE)
  # fit contrast mixed model 
  fit = limma::lmFit(object = v,design = mod0)
  cfit = contrasts.fit(fit = fit, contrasts = c0)
  eb = limma::eBayes(fit = cfit)
  # save results 
  v0[[i]] = v
  fit0[[i]] = fit
  cont0[[i]] = eb
}
names(v0) = names(fit0) = names(cont0) = names(d0d)



# save model fitting data 
saveRDS(object = samplemd, file = paste0(datapath, 'samplemd.rds'))
saveRDS(object = d1, file = paste0(datapath, 'd1.rds'))
saveRDS(object = d7, file = paste0(datapath, 'd7.rds'))
saveRDS(object = d0, file = paste0(datapath, 'd0.rds'))
saveRDS(object = d1d, file = paste0(datapath, 'd1d.rds'))
saveRDS(object = d7d, file = paste0(datapath, 'd7d.rds'))
saveRDS(object = d0d, file = paste0(datapath, 'd0d.rds'))

# save model fits 
# d0
saveRDS(object = fit0, file = paste0(datapath, 'fit0.rds'))
saveRDS(object = cont0, file = paste0(datapath, 'cont0.rds'))
saveRDS(object = v0, file = paste0(datapath, 'v0.rds'))
# day 1 
saveRDS(object = fit1, file = paste0(datapath, 'fit1.rds'))
saveRDS(object = v1, file = paste0(datapath, 'v1.rds'))
# day 7
saveRDS(object = fit7, file = paste0(datapath, 'fit7.rds'))
saveRDS(object = v7, file = paste0(datapath, 'v7.rds'))

# sessioninfo
sessionInfo()
```

Enrichment of curated gene signature pathways based on genes ranked by
vaccination effects above.  
mid_res/1_H1N1_pseudobulk_DE/2_rungsea_day1_day7_V4.r

``` r
# R version 4.0.5
# H1N1 differential expression testing and gene set enrichment analysis
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))

# parallel options for FseaList
BiocParallel::register(BiocParallel::SnowParam(4))
pparam = BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# set data path 
datapath = here("mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/")
dir.create(datapath)

# load pathways to be tested. 
sig_test = readRDS(file = here('signature_curation/combined_sig_sub.rds'))
core_sigs = readRDS(file = here('signature_curation/sig_test_sub.rds'))


# load each time statistical contrast model result extract contrast and rank genes by t statistic 
fit1 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/fit1.rds'))
fit7 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/fit7.rds'))
r1 = ExtractResult(model.fit.list = fit1, what = 'lmer.z.ranks', coefficient.number = 1, coef.name = 'L1')
r7 = ExtractResult(model.fit.list = fit7, what = 'lmer.z.ranks', coefficient.number = 1, coef.name = 'L1')

# run gene set enrichment Day 1 models 
# run unbiased modules and core signatures from past flu studies
g1c = FgseaList(rank.list.celltype = r1, pathways = core_sigs, BPPARAM = pparam)
g1f = FgseaList(rank.list.celltype = r1, pathways = sig_test, BPPARAM = pparam)

# day 7 
g7f = FgseaList(rank.list.celltype = r7, pathways = sig_test, BPPARAM = pparam)

# save 
saveRDS(object = g1c, file = paste0(datapath,'g1c.rds'))
saveRDS(object = g1f, file = paste0(datapath,'g1f.rds'))
saveRDS(object = g7f, file = paste0(datapath,'g7f.rds'))


sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
# [1] viridis_0.5.1            viridisLite_0.3.0        scglmmr_0.1.0            variancePartition_1.25.6
# [5] BiocParallel_1.24.1      limma_3.46.0             magrittr_2.0.1           here_1.0.1
# [9] forcats_0.5.1            stringr_1.4.0            dplyr_1.0.4              purrr_0.3.4
# [13] readr_1.4.0              tidyr_1.1.2              tibble_3.0.6             ggplot2_3.3.3
# [17] tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
# [1] tidyselect_1.1.0            lme4_1.1-26                 RSQLite_2.2.7               AnnotationDbi_1.52.0
# [5] grid_4.0.5                  scatterpie_0.1.7            munsell_0.5.0               codetools_0.2-18
# [9] statmod_1.4.35              withr_2.4.3                 colorspace_2.0-0            GOSemSim_2.16.1
# [13] Biobase_2.50.0              rstudioapi_0.13             stats4_4.0.5                ggsignif_0.6.0
# [17] DOSE_3.16.0                 labeling_0.4.2              MatrixGenerics_1.2.1        Rdpack_2.1.1
# [21] emmeans_1.5.4               GenomeInfoDbData_1.2.4      polyclip_1.10-0             pheatmap_1.0.12
# [25] bit64_4.0.5                 farver_2.0.3                rprojroot_2.0.2             downloader_0.4
# [29] coda_0.19-4                 vctrs_0.4.1                 generics_0.1.2              TH.data_1.0-10
# [33] R6_2.5.0                    doParallel_1.0.16           GenomeInfoDb_1.26.7         graphlayouts_0.7.2
# [37] locfit_1.5-9.4              bitops_1.0-6                cachem_1.0.4                fgsea_1.16.0
# [41] DelayedArray_0.16.3         assertthat_0.2.1            scales_1.1.1                multcomp_1.4-16
# [45] ggraph_2.0.5                enrichplot_1.10.2           gtable_0.3.0                egg_0.4.5
# [49] tidygraph_1.2.0             sandwich_3.0-0              rlang_1.0.2                 slanter_0.2-0
# [53] splines_4.0.5               rstatix_0.7.0               broom_0.7.5                 abind_1.4-5
# [57] BiocManager_1.30.10         reshape2_1.4.4              modelr_0.1.8                backports_1.2.1
# [61] qvalue_2.22.0               clusterProfiler_3.18.1      tools_4.0.5                 ellipsis_0.3.2
# [65] gplots_3.1.1                RColorBrewer_1.1-2          BiocGenerics_0.36.1         Rcpp_1.0.6
# [69] plyr_1.8.6                  progress_1.2.2              zlibbioc_1.36.0             RCurl_1.98-1.3
# [73] prettyunits_1.1.1           ggpubr_0.4.0                cowplot_1.1.1               S4Vectors_0.28.1
# [77] zoo_1.8-8                   SummarizedExperiment_1.20.0 haven_2.3.1                 ggrepel_0.9.1
# [81] fs_1.5.0                    data.table_1.14.0           DO.db_2.9                   openxlsx_4.2.3
# [85] reprex_1.0.0                mvtnorm_1.1-1               packrat_0.7.0               matrixStats_0.58.0
# [89] hms_1.0.0                   GSVA_1.38.2                 xtable_1.8-4                pbkrtest_0.5-0.1
# [93] RhpcBLASctl_0.21-247.1      XML_3.99-0.6                rio_0.5.16                  readxl_1.3.1
# [97] IRanges_2.24.1              gridExtra_2.3               compiler_4.0.5              KernSmooth_2.23-18
# [101] crayon_1.4.1                shadowtext_0.0.9            minqa_1.2.4                 ggfun_0.0.4
# [105] snow_0.4-3                  lubridate_1.7.9.2           DBI_1.1.1                   tweenr_1.0.2
# [109] dbplyr_2.1.0                MASS_7.3-53.1               boot_1.3-27                 Matrix_1.3-2
# [113] car_3.0-10                  cli_3.3.0                   rbibutils_2.0               parallel_4.0.5
# [117] igraph_1.2.6                GenomicRanges_1.42.0        pkgconfig_2.0.3             rvcheck_0.1.8
# [121] foreign_0.8-81              xml2_1.3.2                  foreach_1.5.1               annotate_1.68.0
# [125] XVector_0.30.0              GeneOverlap_1.26.0          estimability_1.3            rvest_0.3.6
# [129] digest_0.6.27               graph_1.68.0                cellranger_1.1.0            fastmatch_1.1-0
# [133] edgeR_3.32.1                GSEABase_1.52.1             curl_4.3                    gtools_3.8.2
# [137] nloptr_1.2.2.2              lifecycle_1.0.0             nlme_3.1-152                jsonlite_1.7.2
# [141] aod_1.3.1                   carData_3.0-4               pillar_1.4.7                lattice_0.20-41
# [145] fastmap_1.1.0               httr_1.4.2                  survival_3.2-10             GO.db_3.12.1
# [149] glue_1.6.2                  zip_2.1.1                   iterators_1.0.13            bit_4.0.4
# [153] ggforce_0.3.3               stringi_1.5.3               blob_1.2.1                  org.Hs.eg.db_3.12.0
# [157] caTools_1.18.1              memoise_2.0.0
```

Generate figures from model results above for day 1 vaccination effects
in unadjuvanted cohort. Derive dell type specific vaccination genes and
shared core interferon signature.  
mid_res/1_H1N1_pseudobulk_DE/3_V4_figures.r

``` r
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
```

Generate figures from model results above for day 7 vaccination effects
in unadjuvanted cohort.  
mid_res/1_H1N1_pseudobulk_DE/4_V4_figures_d7.R

``` r
# R version 4.0.5
# H1N1 differential expression testing and gene set enrichment analysis
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))

# set fig path 
figpath = here("mid_res/1_H1N1_pseudobulk_DE/figuresV4/")



# set theme for subset plots 
mtheme1 = list(
  theme_bw(base_size = 10.5), 
  theme(text = element_text(color = 'black')),
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold",family = "Helvetica"), 
        axis.text.y =  element_text(size = 12, color = 'black'))
  )

# load day 7 mixed model gene set enrichment results 
g7f = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/g7f.rds'))

# filter subset of signals to visualize 
g7f = lapply(g7f, function(x) 
  x %>%  
    filter(!str_sub(pathway, 1,5) == 'REACT' ) %>% 
    filter(NES > 0)
   )

# save global 
p = PlotFgsea(gsea_result_list = g7f, p.threshold = 0.01, NES_filter = 0.1)
ggsave(p, filename = paste0(figpath, 'gsea.pos.d7.pdf'), width = 8, height = 9)

# cell type specific signals to highlight.
d7d = p$data
d7d$pathway  = as.character(d7d$pathway)
# define B cell signals 
bsub  = c(
  'LI.S2 B cell surface signature',
  'LI.M47.0 enriched in B cells (I)',
  'CD40_ACT',
  'LI.M5.0 regulation of antigen presentation and immune response',
  'LI.S8 Naive B cell surface signature',
  'KEGG_AMINOACYL_TRNA_BIOSYNTHESIS',
  'KEGG_OXIDATIVE_PHOSPHORYLATION',
  'LI.M212 purine nucleotide biosynthesis',
  'LI.M234 transcription elongation, RNA polymerase II',
  'LI.M32.0 platelet activation (I)',
  'LI.M227 translation initiation',
  'LI.M37.0 immune activation - generic cluster'
)

bsub.plot = d7d %>% filter(celltype == 'BC_Naive' & pathway %in% bsub)

# give shorter names 
bsub.plot$pathway[bsub.plot$pathway == 'LI.M5.0 regulation of antigen presentation and immune response'] <-
  'LI.M5.0 regulation of antigen presentation'
bsub.plot$pathway[bsub.plot$pathway == 'KEGG_AMINOACYL_TRNA_BIOSYNTHESIS'] <-
  'kegg aminoacyl tRNA biosynthesis'
bsub.plot$pathway[bsub.plot$pathway == 'KEGG_OXIDATIVE_PHOSPHORYLATION'] <- 
  'kegg oxidatie phosphorylation'

# save plot 
p = ggplot(bsub.plot, aes(x = NES, y = reorder(pathway, NES), size = -log10(padj)) ) + 
  mtheme1 +
  ylab("") +
  xlab('Normalized Enrichment Score') + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  xlim(-1,3) +
  geom_point(shape = 21 , fill ='red' ) + 
  scale_size_area() + 
  ggtitle('Day 7 induced: Naive B cells')
ggsave(p,filename = paste0(figpath,'bnaive.d7.pdf'), width = 7, height = 3)


# t cell (EM )
tsub.plot = d7d %>%  filter(celltype == 'CD4_Efct_Mem_Tcell')
# give shorter names 
tsub.plot$pathway = as.character(tsub.plot$pathway)
tsub.plot$pathway[tsub.plot$pathway == 'KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION'] <- 
  'kegg valine leucine isoleucine degratation'
tsub.plot$pathway[tsub.plot$pathway == 'KEGG_PEROXISOME'] <- 'kegg peroxisome'
tsub.plot$pathway[tsub.plot$pathway == 'KEGG_FATTY_ACID_METABOLISM'] <- 
  'kegg fatty acid metabolism'
tsub.plot$pathway[tsub.plot$pathway == 'KEGG_PRIMARY_IMMUNODEFICIENCY'] <- 
  'kegg primary immunodeficiency'

# save plot 
p = ggplot(tsub.plot, aes(x = NES, y = reorder(pathway, NES), size = -log10(padj)) ) + 
  mtheme1 +
  ylab("") +
  xlab('Normalized Enrichment Score') + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  xlim(-1,3) +
  geom_point(shape = 21 , fill ='red' ) + 
  scale_size_area() + 
  ggtitle('Day 7 induced: CD4 effector memory T cells')
p
ggsave(p,filename = paste0(figpath,'cd4mem.d7.pdf'), width = 7, height = 3)
```

Visualize core shared across subsets interferon signature defined in
script 3.  
mid_res/1_H1N1_pseudobulk_DE/5_shared_core_fin_state.r

``` r
# R version 4.0.5
# H1N1 differential expression testing and gene set enrichment analysis
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
source('functions/MattPMutils.r')
# set fig path 
figpath = here("mid_res/1_H1N1_pseudobulk_DE/figuresV4/")

# load model fitting data and pseudobulk data saved in mixed model workflow (scipt 1)
d1d = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/d1d.rds'))

# log CPM the pseudobulk data 
lcpm = lapply(d1d, edgeR::cpm, log = TRUE)  

# read Core signature
ms1.names = readRDS(file = here("mid_res/1_H1N1_pseudobulk_DE/dataV4/ms1.names.rds"))

av2 = lapply(lcpm, function(x){ x[ms1.names, ] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    mutate(timepoint = str_sub(sample, -2, -1))
})

# for (i in 1:length(av2)) {
#   av2[[i]]$celltype = names(av2)[i]
# }
av2 = av2 %>% bind_rows(.id = 'celltype')
av3 = av2 %>% gather(gene, average, ADAR:ZBP1)
av4 = av3 %>% group_by(sample, celltype, timepoint) %>%
  summarize(core_ifn_score = mean(average)) 

av4 = av4 %>% mutate(subject = str_sub(sample, 1,3))
av4$celltype = str_replace_all(string = av4$celltype,pattern = '_', replacement = ' ')

# colors 
cu1 = unname(sapply(c("grey", "orange"), col.alpha, 0.5))
cu2 = c("grey", "orange")
p=
  ggplot(av4, aes(x = celltype, y = core_ifn_score ,  color = timepoint, fill = timepoint )) +
  theme_bw() +
  geom_boxplot() +
  #scale_x_discrete(position = "top") +
  ylab("shared day 1 induced IFN state") +
  xlab("")+
  scale_color_manual(values = cu2 ) +
  scale_fill_manual(values = cu1 ) +
  theme(axis.title.y  = element_text(size = 12, color = 'black' )) +
  theme(axis.text.x = element_text(size = 11, color = "black") )  +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, color = 'black'))+ 
  theme(legend.position = 'top')
p
ggsave(p, filename = paste0(figpath, "core_interferon_state.pdf"), width = 6, height = 5)



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
```

### Fig S2 visualization of day 7 post vaccination phenotypes and predictive signature deconvolution <a name="fig2.2"></a>

Comparison of day 7 signatures predictive of antibody response in
microarray and aggregated CITE-seq data  
mid_res/array_bulk_comparison/array_bulk_comparison.R

``` r
# R 4.0.5 
suppressMessages(library(tidyverse))
suppressMessages(library(scglmmr))
suppressMessages(library(here))
suppressMessages(library(BiocParallel))
suppressMessages(library(edgeR))
suppressMessages(library(variancePartition))

# set paths 
figpath = here('mid_res/array_bulk_comparison/figures/')
dir.create(figpath)
datapath = here('mid_res/array_bulk_comparison/generated_data/')
dir.create(datapath)

# set parallel options 
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

#pb = readRDS(file = here('mid_res/pb.ds'))
#fit7 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/fit7.rds'))

core7 = readRDS(file = here('signature_curation/core_d7.rds'))
core7$`LI.M156.0 plasma cell b cell Ig`


# read processed pseudobulk data 
# subset to unadjuvanted cohort and remove cell type string from sample names 
s = readRDS(here('data/h1h5_annotated_with_meta.rds'))


# day1.cohort = c("200", "205","207", "236", "237", "250", "273" ,"279")
day7.cohort = c("201", "209", "212", "215", "229", "233",
                "234", "245", "256",  "261",  "268", "277")

md7 = s@meta.data %>%  
  filter(cohort == 'H1N1') %>% 
  filter(sampleid %in% day7.cohort) %>% 
  arrange(sampleid, timepoint)
umi7 = s@raw.data[ ,rownames(md7)]

# pseudobulk all cells
scell = lapply(X = split(md7, f = md7$sample), FUN = rownames)
csample = lapply(scell, function(x) Matrix::rowSums(umi7[ ,x]))
pbmat = as.data.frame(t(do.call(cbind, csample))) %>% t()


#define day 7 metadata 
met = md7 %>% 
  select(sample,timepoint , subjectid = sampleid) %>% 
  group_by(sample, subjectid, timepoint) %>% 
  distinct() %>%
  ungroup() %>% 
  column_to_rownames('sample')
met$timepoint = factor(met$timepoint, levels = c('d0', 'd7'))

# check order 
stopifnot(isTRUE(all.equal(colnames(pbmat), rownames(met))))

# filter features  
gene.keep = filterByExpr(y = pbmat, design = met$timepoint,min.count = 3)
pbmat = pbmat[gene.keep, ]

# fit model 
f1 = ~ 0 + timepoint + (1|subjectid)
L1 = makeContrastsDream(formula = f1, data =  met,
                        contrasts =  "timepointd7 - timepointd0")
v7 = voomWithDreamWeights(counts = pbmat,formula = f1,
                          BPPARAM = pparam,data = met)
result7 = dream(exprObj =  v7,formula = f1,data = met,
                L = L1, BPPARAM = pparam, useWeights = TRUE)
# save 
saveRDS(result7,file = paste0(datapath,'result7.rds'))


## Part II fit same model on microarray data 

# array data coefficient 
# read array data 
array = data.table::fread("data/CHI_H1N1_data/microarray/CHI_GE_matrix_gene.txt", data.table = F) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("gene") %>% 
  select(-matches("day70")) %>% 
  select(., matches("day0|day7")) %>% 
  select(-matches("pre")) 

# day 7 samples; no day 7 data for subject 209 
array7 =
  array %>% 
  select(which(substr(names(.),1,3) %in% day7.cohort)) %>% 
  select(-matches("209"))

# Metadata 
d7md  = 
  colnames(array7) %>% 
  base::as.data.frame() %>% 
  rename(sample = ".") %>% 
  mutate(timepoint = str_sub(sample, -4,-1)) %>% 
  mutate(subjectid = str_sub(sample, 1,3)) %>% 
  column_to_rownames("sample")

d7md$timepoint = factor(d7md$timepoint, levels = c('day0', 'day7'))

# check order 
stopifnot(isTRUE(all.equal(colnames(array7), rownames(d7md))))

# test same genes
gene.sub = rownames(pbmat)
array7 = array7[gene.sub, ]
gene.keep2 = !is.na(Matrix::rowSums(array7))
array7 = array7[gene.keep2, ]


# fit model
L1.1 = makeContrastsDream(formula = f1, data =  d7md,
                          contrasts =  "timepointday7 - timepointday0")
# no weights for normalized microarray data ; same formula
result7.1 = dream(exprObj =  array7, formula = f1,data = d7md,
                L = L1.1, BPPARAM = pparam, useWeights = FALSE)
# save 
saveRDS(result7.1,file = paste0(datapath,'result7.1.rds'))


# comparison 

# extract results for array data 
ra = ExtractResult(model.fit.list = list('array' = result7.1), 
                   coefficient.number = 1,
                   coef.name = 'timepointday7 - timepointday0')
ra$array$logFC.array = ra$array$logFC
ra = ra$array %>%  
  select(gene,  logFC.array)
# extract results for CITE-seq bulk data 
rc = ExtractResult(model.fit.list = list('CITE-seq Bulk' = result7),
                   coefficient.number = 1, 
                   coef.name = 'timepointd7 - timepointd0')
rc$`CITE-seq Bulk`$logFC.CITEseq = rc$`CITE-seq Bulk`$logFC
rc = rc$`CITE-seq Bulk` %>% 
  select(gene,  logFC.CITEseq)


# visualize correlation between signals 
d = full_join(ra, rc)
dsub = d %>% filter(gene %in% core7$`LI.M156.0 plasma cell b cell Ig`)
p = 
  ggplot(dsub, aes(x = logFC.array , y = logFC.CITEseq)) + 
  theme_bw() + 
  geom_smooth(method = 'lm', color = 'black') + 
  geom_point(shape = 21, color = 'white', fill = 'black') + 
  ggrepel::geom_text_repel(data = dsub, mapping = aes(label = gene), 
                           segment.size = 0.1, box.padding = 0.1, 
                           max.overlaps = 4,  size = 2.6) + 
  ggpubr::stat_cor(method = "pearson") + 
  ylab("LI.M156 Module \n CITE-seq Day 7 log2 FC") + 
  xlab("LI.M156 Module \n Microarray Day 7 log2 FC") 
ggsave(p, filename = paste0(figpath, "/m156_gene_correlation.pdf"),width = 3, height = 3)
```

Single cell deconvolution of predictive signatures across protein based
subsets.  
mid_res/d7_predictive_deconvolution/1.d7predictive.score.deconvolution.r  
*This script uses R 3.5.1*

``` r
# Day 7 signature deconvolution direct with module score.  
# R 3.5.1 
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(magrittr))
suppressMessages(library(here))
source("functions/analysis_functions.R")
source('functions/MattPMutils.r')

# save path 
figpath = here("mid_res/d7_predictive_deconvolution/figures/")
dir.create(figpath)

# day 7 core signatures. 
sig7 = readRDS("signature_curation/core_d7.rds")

# h1 data baseline and day 7 cells. 
h1 = ReadCohort(joint_object_dir = "data/h1h5_annotated_with_meta.rds", cohort = "H1N1")
h1 = SetAllIdent(h1, id = "time_cohort") %>% 
  SubsetData(ident.use = "d7")

# add module score 
h1 = AddModuleScore(h1, genes.list = sig7, seed.use = 1, enrich.name = names(sig7))

# get long for mfor visualization of module score distribtion. 
df_sig = h1@meta.data %>% select(
  celltype_joint,
  timepoint,
  LI.M156_Plasma_Cell = `LI.M156.0 plasma cell b cell Ig6`,
  CHI_4 = 'CHI 4 d710',
  CHI_5 = `CHI 5 d711`,
  CHI_d7_Response = `CHI d7 Response9`
) %>%
  mutate(celltype_b = if_else(
    celltype_joint %in% c("BC_Mem", "BC_Naive",  "CD38_Bcell", "pDC"),
    true = celltype_joint,
    false = "other"
  )) %>%
  gather(module, module_score, LI.M156_Plasma_Cell:CHI_d7_Response) %>%
  mutate(timepoint = factor(timepoint, levels = c("d0", "d7"))) %>%
  mutate(celltype_b = str_replace_all(
    string = celltype_b,
    pattern = "_",
    replacement = " "
  )) %>%
  mutate(celltype_joint = str_replace_all(
    string = celltype_joint,
    pattern = "_",
    replacement = " "
  )) %>%
  mutate(celltype_b = factor(
    celltype_b,
    levels = c("CD38 Bcell", "pDC", "BC Naive", "BC Mem", "other")
  ))

# plot CHI predictive sig from array 
subplot = df_sig %>% filter(timepoint == "d7" &
                              module %in% c("CHI_d7_Response", "CHI_4", "LI.M156_Plasma_Cell"))
subplot = subplot %>% mutate(
  module = dplyr::recode(
    module,
    "LI.M156_Plasma_Cell" = "M156",
    "CHI_d7_Response" = "Antibody sig",
    "CHI_4" = "CHI 4"
  )
)

grey1 = col.alpha(acol = 'grey',alpha = 0.2)

p = ggplot(subplot,   aes(
  x = reorder(celltype_joint, module_score),
  y = module_score,
  fill = celltype_joint
)) +
  theme_bw() +
  geom_violin(
    show.legend = FALSE,
    scale = "width", 
    size = 0.2, 
    draw_quantiles = c(0.5)
  ) +
  facet_wrap( ~ module,
              as.table = TRUE,
              scales = "free_x",
              ncol = 4) +
  geom_hline(yintercept  = 0,
             size = 0.5,
             linetype = "dashed") +
  coord_flip() +
  scale_fill_manual(values = c(rep(grey1, 5), "red", rep(grey1, 17)))  +
  theme(
    panel.spacing = unit(0.2, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 7, color = "black"),
    panel.border = element_blank()
  ) +
  xlab("") + ylab("single cell score distribution \n day 7 bulk predictive signatures") +
  ggtitle("Day 7 post-vaccination")
ggsave(p, filename = paste0(figpath,"d7core_Tirosh_mod_score_celltypes_sub.pdf"), width = 5, height = 5)
```

Single cell raw UMI count deconvolution of TNFRSF17 gene.  
mid_res/d7_predictive_deconvolution/2.TNFRSF17.deconvolution.r4.0.5.r

``` r
# R 4.0.5
suppressMessages(library(tidyverse))
suppressMessages(library(scglmmr))
suppressMessages(library(here))
suppressMessages(library(magrittr))
source('functions/analysis_functions.R')

# save path 
figpath = here('mid_res/d7_predictive_deconvolution/figures/')

# single cell composition of m156 
core7 = readRDS(file = here('signature_curation/core_d7.rds'))
h1 = readRDS(file = here('data/h1h5_annotated_with_meta.rds'))

m156_sub = intersect(rownames(h1@data), core7$`LI.M156.0 plasma cell b cell Ig`)
gene = as.data.frame(as.matrix(t(as.matrix(h1@data[m156_sub, ]))))
#pdf = cbind(gene, h1@meta.data)
#pdf = pdf %>% gather(gene, normcount, m156_sub[1]:m156_sub[length(m156_sub)])

# composisiton of TNFRSF17
gene = as.data.frame(as.matrix(t(as.matrix(h1@raw.data[m156_sub, ]))))
pdf = cbind(gene, h1@meta.data)
pdf = pdf %>% gather(gene, count, m156_sub[1]:m156_sub[length(m156_sub)])
tnf = pdf %>% 
  filter(cohort == 'H1N1') %>% 
  filter(gene == "TNFRSF17" & timepoint %in% c("d0", "d7") ) %>% 
  group_by(celltype_joint, timepoint, gene) %>% 
  summarise(n = sum(count)) %>%
  ungroup() 


# visualize distribution
tnf$celltype_joint = str_replace_all(tnf$celltype_joint, pattern = '_',replacement = ' ')
p = ggplot(tnf, aes(x = reorder(celltype_joint, n), y = n, fill = timepoint)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  theme_bw() +
  coord_flip() + 
  theme(axis.text.x =  element_text(size = 9, color = 'black'), axis.title.x = element_text(color = 'black')) + 
  theme(axis.text.y =  element_text(size = 9, color = 'black'), axis.title.y = element_text(color = 'black')) + 
  ylab("UMI counts") + 
  scale_fill_manual(values = c("grey48", "red")) + 
  theme(axis.title.y  =  element_blank()) + 
  theme(legend.position = c(0.7, 0.2) ) + 
  ggtitle("TNFRSF17 Gene") 
p
ggsave(p, filename = paste0(figpath, "TNFRSF17_composition.pdf"),width = 3, height = 4.3)

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] magrittr_2.0.3    here_1.0.1        scglmmr_0.1.0     forcats_0.5.1     stringr_1.4.0     dplyr_1.0.4       purrr_0.3.4      
# [8] readr_1.4.0       tidyr_1.1.2       tibble_3.1.8      ggplot2_3.3.3     tidyverse_1.3.0   viridis_0.5.1     viridisLite_0.3.0
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.2                  tidyselect_1.2.0            lme4_1.1-26                 RSQLite_2.2.7              
# [5] AnnotationDbi_1.52.0        grid_4.0.5                  BiocParallel_1.24.1         scatterpie_0.1.7           
# [9] munsell_0.5.0               codetools_0.2-18            statmod_1.4.35              withr_2.4.3                
# [13] colorspace_2.0-0            GOSemSim_2.16.1             Biobase_2.50.0              rstudioapi_0.13            
# [17] ROCR_1.0-11                 stats4_4.0.5                ggsignif_0.6.0              DOSE_3.16.0                
# [21] labeling_0.4.2              MatrixGenerics_1.2.1        Rdpack_2.1.1                emmeans_1.5.4              
# [25] GenomeInfoDbData_1.2.4      polyclip_1.10-0             bit64_4.0.5                 farver_2.0.3               
# [29] pheatmap_1.0.12             rprojroot_2.0.2             downloader_0.4              coda_0.19-4                
# [33] vctrs_0.5.1                 generics_0.1.2              TH.data_1.0-10              R6_2.5.0                   
# [37] doParallel_1.0.16           GenomeInfoDb_1.26.7         graphlayouts_0.7.2          locfit_1.5-9.4             
# [41] bitops_1.0-6                cachem_1.0.4                fgsea_1.16.0                DelayedArray_0.16.3        
# [45] assertthat_0.2.1            scales_1.1.1                multcomp_1.4-16             ggraph_2.0.5               
# [49] enrichplot_1.10.2           gtable_0.3.0                egg_0.4.5                   tidygraph_1.2.0            
# [53] sandwich_3.0-0              rlang_1.0.6                 slanter_0.2-0               splines_4.0.5              
# [57] rstatix_0.7.0               broom_0.7.5                 BiocManager_1.30.10         reshape2_1.4.4             
# [61] abind_1.4-5                 modelr_0.1.8                backports_1.2.1             qvalue_2.22.0              
# [65] clusterProfiler_3.18.1      tools_4.0.5                 ellipsis_0.3.2              gplots_3.1.1               
# [69] RColorBrewer_1.1-2          BiocGenerics_0.36.1         Rcpp_1.0.9                  plyr_1.8.6                 
# [73] progress_1.2.2              zlibbioc_1.36.0             RCurl_1.98-1.3              prettyunits_1.1.1          
# [77] ggpubr_0.4.0                cowplot_1.1.1               S4Vectors_0.28.1            zoo_1.8-8                  
# [81] SummarizedExperiment_1.20.0 haven_2.4.3                 ggrepel_0.9.1               fs_1.5.0                   
# [85] variancePartition_1.25.6    data.table_1.14.0           DO.db_2.9                   openxlsx_4.2.3             
# [89] RANN_2.6.1                  reprex_1.0.0                mvtnorm_1.1-1               packrat_0.7.0              
# [93] matrixStats_0.58.0          hms_1.0.0                   GSVA_1.38.2                 xtable_1.8-4               
# [97] pbkrtest_0.5-0.1            RhpcBLASctl_0.21-247.1      XML_3.99-0.6                rio_0.5.16                 
# [101] readxl_1.3.1                IRanges_2.24.1              gridExtra_2.3               compiler_4.0.5             
# [105] KernSmooth_2.23-18          crayon_1.4.1                shadowtext_0.0.9            minqa_1.2.4                
# [109] ggfun_0.0.4                 lubridate_1.8.0             DBI_1.1.1                   tweenr_1.0.2               
# [113] dbplyr_2.1.0                MASS_7.3-53.1               boot_1.3-27                 Matrix_1.4-1               
# [117] car_3.0-10                  cli_3.4.1                   rbibutils_2.0               parallel_4.0.5             
# [121] igraph_1.2.6                GenomicRanges_1.42.0        pkgconfig_2.0.3             rvcheck_0.1.8              
# [125] foreign_0.8-81              xml2_1.3.2                  foreach_1.5.1               annotate_1.68.0            
# [129] XVector_0.30.0              GeneOverlap_1.26.0          estimability_1.3            rvest_0.3.6                
# [133] digest_0.6.27               graph_1.68.0                cellranger_1.1.0            fastmatch_1.1-0            
# [137] edgeR_3.32.1                GSEABase_1.52.1             curl_4.3                    gtools_3.8.2               
# [141] nloptr_1.2.2.2              lifecycle_1.0.3             nlme_3.1-152                jsonlite_1.7.2             
# [145] aod_1.3.1                   carData_3.0-4               limma_3.46.0                fansi_0.4.2                
# [149] pillar_1.8.1                lattice_0.20-41             fastmap_1.1.0               httr_1.4.2                 
# [153] survival_3.2-10             GO.db_3.12.1                glue_1.6.2                  zip_2.1.1                  
# [157] iterators_1.0.13            bit_4.0.4                   ggforce_0.3.3               stringi_1.5.3              
# [161] blob_1.2.1                  org.Hs.eg.db_3.12.0         caTools_1.18.1              memoise_2.0.0
```

### Fig 2. bottom up single cell reconstruction of single cell monocyte pseudotime <a name="fig2.3"></a>

*This section is run with R 3.5.1*

Construct mRNA-based monocyte single cell latent space with DDRTree.
Infer pseudotime using monocle.

mid_res/monocyte_map/1_monocyte_map.r

``` r
# R 3.5 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
suppressMessages(library(here))
#source("de_workflow-master/downstream_analysis_functions.r")
source("functions/analysis_functions.R")
btm = readRDS("signature_curation/BTM_li.rds")
datapath = here('mid_res/monocyte_map/generated_data/'); dir.create(datapath)

# select celltype on which to to run pseudotime analysis 
celltype_use = c("CD14_Mono", "CD16_Mono")

# subset time cohort 
sub = ReadCohort(joint_object_dir = "data/h1h5_annotated_with_meta.rds", cohort = "H1N1") %>% 
  SetAllIdent(id = "time_cohort") %>% 
  SubsetData(ident.use  = "d1", subset.raw = TRUE) %>% 
  SetAllIdent(id = "celltype_joint") %>% 
  SubsetData(ident.use = celltype_use, subset.raw = TRUE)

# addd proteins as meta data 
prot_dat = as.data.frame(t(sub@assay$CITE@data))
sub = AddMetaData(sub, metadata = prot_dat)

#### use DDR tree to calculate trajectory 
library(monocle)
sm = monocle::importCDS(sub)
sm = BiocGenerics::estimateSizeFactors(sm)
sm <- detectGenes(sm, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(sm),num_cells_expressed >= 15))

## Select genes 
time1_genes = differentialGeneTest(sm[expressed_genes,],fullModelFormulaStr = "~timepoint",cores = 4)
rpgene =  grep(pattern = "RPL|RPS|MT-|RP11", x = expressed_genes, value = TRUE)
t1genes = time1_genes %>% 
  rownames_to_column("gene") %>% 
  filter(qval < 0.15) %>% 
  arrange(qval) %>% 
  filter(!gene %in% rpgene)
t1gene = t1genes %$% gene

# set ordering filter and reduce dimensions by the ddr tree algorithm 
msm = setOrderingFilter(sm, ordering_genes = t1gene)
sm = reduceDimension(sm, max_components = 2, 
                     reduction_method = "DDRTree", 
                     residualModelFormulaStr = "~ sampleid")
sm = orderCells(sm, reverse = TRUE)
saveRDS(sm, file = paste0(datapath, "sm_cd14_cd16_d1_monocle_object.rds"))
```

Integrate pseudotime, surface protein levels and time relative to
vaccination with mRNA based pseudotime calculated above to interpret 3
branches.  
mid_res/monocyte_map/2_monocyte_mat_visualization.I.r

``` r
# R 3.5 
# make visualization of main monocyte pseudotime axis. 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
suppressMessages(library(here))
suppressMessages(library(monocle))
source("functions/analysis_functions.R")
source('functions/MattPMutils.r')
btm = readRDS("signature_curation/BTM_li.rds")
figpath = here("mid_res/monocyte_map/figures/"); dir.create(figpath, recursive = TRUE)
sm = readRDS(file = here("mid_res/monocyte_map/generated_data/sm_cd14_cd16_d1_monocle_object.rds"))

######
# visualization 
p = plot_cell_trajectory(sm, color_by = "Pseudotime")
df = ggplot_build(p)[["plot"]][["data"]] %>% 
  rename(component_1 = data_dim_1  ,component_2 = data_dim_2)
df$adjmfc.time = factor(df$adjmfc.time, levels = c("d0 low", "d1 low", "d0 high", "d1 high"))
library(cowplot) 

# set plot theme 
theme.set = list(theme_bw(), 
                 theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(), 
                       axis.ticks.y = element_blank(),
                       axis.ticks.x = element_blank()))
# reverse order to make the root node go left to right 
df$component_1 = -1* df$component_1

# create main plot
p1 = ggplot(df, aes(x = component_1, y = component_2, fill = Pseudotime)) +
  theme.set + 
  geom_point(shape = 21, size = 2.2, color = "grey" ,stroke = 0.1) +
  theme(legend.position = c(0.1, 0.7)) + 
  theme(strip.background = element_blank()) + 
  ylab("mRNA trajectory component 2") + 
  xlab("mRNA trajectory component 1") + 
  scale_fill_viridis_c(option = "B") + 
  scale_color_manual(values = c("grey", "black"))
ggsave(p1, filename = paste0(figpath, "mono_trajectory_only.png"), width = 6, height = 5)

# make a background plot on which to add the canvas marginal plots
pnull = ggplot(df, aes(x = component_1, y = component_2, fill = Pseudotime)) +
  theme.set + 
  theme(legend.position = c(0.1, 0.7)) + 
  theme(strip.background = element_blank()) + 
  theme(axis.title.x = element_text(size = 17)) + 
  theme(axis.title.y = element_text(size = 17)) +
  ylab("mRNA trajectory component 2") + 
  xlab("mRNA trajectory component 1") + 
  scale_fill_viridis_c(option = "plasma") + 
  scale_color_manual(values = c("grey", "black"))

# real time for top margin 
time.col = c( 
  col.alpha(acol = 'black', 0.1), 
  col.alpha(acol = ggsci::pal_jama()(2)[2], 0.4) 
  ) 
timecol2 = c(
  'black',
  ggsci::pal_jama()(2)[2] 
)
xd = axis_canvas(pnull, axis = "x") +
  geom_density(data = df, aes(x = component_1, color = timepoint), size = 1) + 
  scale_color_manual(values = time.col)

## test 
xd = axis_canvas(pnull, axis = "x") +
  geom_density(data = df, aes(x = component_1, fill = timepoint, color = timepoint), size = 1) + 
  scale_fill_manual(values = time.col) + 
  scale_color_manual(values = timecol2)
xd
## 


# CD14 vs CD16 protein for bottom margin 
x16 = axis_canvas(pnull, axis = "x") +
  geom_smooth(data = df, 
              aes(x = component_1, y = CD16_PROT), 
              method = "loess", se = TRUE,
              color = 'black') + 
  geom_smooth(data = df,
              aes(x = component_1, y = CD14_PROT),
              method = "loess", se = TRUE, 
              color = 'grey') 
# add to plot   
p2 <- insert_xaxis_grob(pnull, xd, grid::unit(.2, "null"), position = "top")
p4 =  insert_xaxis_grob(p2, x16, grid::unit(.3, "null"), position = "bottom")
p6 = ggdraw(p4)
p6
# save 
ggsave(p6, filename = paste0(figpath, "monocyte_merged_plot.pdf"), width = 6, height = 8)
ggsave(p6, filename = paste0(figpath, "monocyte_merged_plot.png"), width = 6, height = 8)
```

Integrate “bottom up” single cell monocyte pseudotime reconstruction
with “top down” mixed effects vaccine perturbation phenotypes. Within
the genes defined by the pseudobulk mixed effects models (the section
above) and in the leading edge of curated pathway enrichments, calculate
branch dependent differential expression using BEAM. Define categories
of genes based on their behavior across single cell pseudotime.
mid_res/monocyte_map/3_monocyte_perturbation_integration.r

``` r
# R 3.5 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
suppressMessages(library(here))
library(monocle)
source("functions/analysis_functions.R")
source("functions/MattPMutils.r")
figpath = here("mid_res/monocyte_map/figures/"); dir.create(figpath, recursive = TRUE)
datapath = here("mid_res/monocyte_map/generated_data/"); dir.create(datapath, recursive = TRUE)

# load monocle object 
sm = readRDS(file = here("mid_res/monocyte_map/generated_data/sm_cd14_cd16_d1_monocle_object.rds"))

# load day 1 enrichment from monocytes 
g1c = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/g1c.rds'))
g1 = g1c$CD14_Mono %>% 
  filter(NES > 0) %>% 
  filter(padj < 0.05)
monole = g1$leadingEdge
names(monole) = g1$pathway
cgene2 = unique(unlist(monole))
saveRDS(cgene2,file = paste0(datapath, 'cgene2.rds'))

# define branch dependent genes
de.branch <- BEAM(sm[cgene2, ], branch_point = 1, cores = 4)
saveRDS(de.branch,file = paste0(datapath,'de_branch.rds'))
de.branch = readRDS(file = here('mid_res/monocyte_map/generated_data/de_branch.rds'))
de.branch.sub = de.branch %>% filter(qval < 0.05)
branch.genes = as.character(de.branch.sub$gene_short_name)


# Visualization of branch genes 
led = exprs(sm)[branch.genes, ] %>% 
  t() %>% 
  as.matrix() %>% 
  as.data.frame()
d = cbind(led, sm@phenoData@data)

dd = d %>% select(OASL, CCL2, IFITM2, FCER1G, TNFSF10, FCGR1B)
dd$Pseudotime = d$Pseudotime
dd$timepoint = d$timepoint

# mono act 
dd = dd %>%  filter(Pseudotime > 5)
dd = dd %>% gather(gene, value, OASL:FCGR1B)

# vis theme 
mtheme =  list(
  theme_bw(), 
  geom_smooth(size = 2), 
  theme(axis.title = element_text(size = 18)), 
  ggsci::scale_color_jama(alpha = 0.8),
  theme(legend.position = c(0.2, 0.8)),
  xlab('Pseudotime') 
)

# example category 1 gene 
p1 = ggplot(dd %>% filter(gene %in% c('CCL2')), 
            aes(x = Pseudotime, y = value,  color = timepoint)) + 
  ylab('CCL2 Expression') +
  mtheme 
ggsave(p1, filename = paste0(figpath,'Category1_1_CCL2.pdf'), width = 3.5, height = 3.5)

# example category 2 gene 
p1 = ggplot(dd %>% filter(gene %in% c('TNFSF10')), 
            aes(x = Pseudotime, y = value,  color = timepoint)) + 
  ylab('TNFSF10 Expression') +
  mtheme 
ggsave(p1, filename = paste0(figpath,'Category1_2_TNFSF10.pdf'), width = 3.5, height = 3.5)

# example category 2 gene
p1 = ggplot(dd %>% filter(gene %in% c('FCER1G')), 
            aes(x = Pseudotime, y = value,  color = timepoint)) + 
  ylab('FCER1G Expression') +
  mtheme 
ggsave(p1, filename = paste0(figpath,'Category2_1_FCER1G.pdf'), width = 3.5, height = 3.5)

# example category 2 gene
p1 = ggplot(dd %>% filter(gene %in% c('IFITM2')), 
            aes(x = Pseudotime, y = value,  color = timepoint)) + 
  mtheme + 
  ylab('IFITM2 Expression')
ggsave(p1, filename = paste0(figpath,'Category2__1_IFITM2.pdf'), width = 3.5, height = 3.5)



#########################
# leading edge signatures analysis 

# Interferon 
sig = intersect(branch.genes, monole$`reactome interferon signaling`)
pdf(file = paste0(figpath, 'IFN_branch_de.pdf'),width = 5, height = 5)
plot_genes_branched_heatmap(sm[sig, ],
                            branch_point = 1,
                            cores = 1,
                            num_clusters = 3,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()

# set categ 
cat2 = c('IFITM2', 'PTPN1', 'EIF4E2', 'IFITM3', 'HLA-C')
cat1 = sig[!sig %in% cat2]
# get data 
led = exprs(sm)[sig, ] %>% 
  t() %>% 
  as.matrix() %>% 
  as.data.frame()
d = cbind(led, sm@phenoData@data)
dd = d %>% select(sig)
dd = apply( dd, 2, scale.simple) %>%  as.data.frame()
dd$Pseudotime = d$Pseudotime
dd$timepoint = d$timepoint
index1 = sig[1]
index2 = sig[length(sig)]
d3 = dd %>% gather(gene, value, index1:index2 ) 
d3$cat = ifelse(d3$gene %in% cat1, '1', '2')



p1 = ggplot(data = d3 %>%  
              filter(Pseudotime > 5 & cat ==1 )) + 
  geom_smooth(data = d3 %>% 
                filter(Pseudotime > 5 & timepoint == 'd0' & cat ==1 ), 
              mapping = aes(x  = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(1), se = FALSE) + 
  geom_smooth(data = d3 %>% filter(Pseudotime > 5 & timepoint == 'd1' & cat ==1 ),
              mapping =  aes(x  = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(2)[2],
              se = FALSE, size = 2) + 
  theme_bw() + 
  geom_vline(xintercept = 9.5, linetype = 'dashed') + 
  theme(strip.background = element_blank()) + 
  ylab('reactome interferon \n Category 1 genes ') + 
  theme(axis.title = element_text(size = 18))
p1
ggsave(p1,filename = paste0(figpath, 'IFNcat1.pdf'), width = 3.7, height = 3.5)


p2 = ggplot(data = d3 %>%  
              filter(Pseudotime > 5 & cat ==2 )) + 
  geom_smooth(data = d3 %>% 
                filter(Pseudotime > 5 & timepoint == 'd0' & cat ==2 ), 
              mapping = aes(x  = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(1), se = FALSE) + 
  geom_smooth(data = d3 %>% 
                filter(Pseudotime > 5 & timepoint == 'd1' & cat ==2),
              mapping =  aes(x  = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(2)[2], 
              se = FALSE, size = 2) + 
  theme_bw() + 
  geom_vline(xintercept = 9.5, linetype = 'dashed') + 
  theme(strip.background = element_blank()) + 
  ylab('reactome interferon \n Category 2 genes ') + 
  theme(axis.title = element_text(size = 18))
p2
ggsave(p2,filename = paste0(figpath, 'IFNcat2.pdf'), width = 3.7, height = 3.5)


# mtor hypoxia 
sig = intersect(branch.genes,
                c(monole$`HALLMARK hypoxia`, monole$`HALLMARK MTORC1 signaling`)
                )

pdf(file = paste0(figpath, 'mtorhypoxia_branch_de.pdf'),width = 5, height = 5)
plot_genes_branched_heatmap(sm[sig, ],
                            branch_point = 1,
                            cores = 1,
                            num_clusters = 3,
                            use_gene_short_name = T,
                            show_rownames = T)

dev.off()

#set categ 
cat2 = c('CTSC', 'PFKL', 'ACTR3', 'CITED2', 'PGK1', 'INSIG1', 'CHST2')
cat1 = sig[!sig %in% cat2]

# get data 
led = exprs(sm)[sig, ] %>% 
  t() %>% 
  as.matrix() %>% 
  as.data.frame()
d = cbind(led, sm@phenoData@data)
dd = d %>% select(sig)
dd = apply( dd, 2, scale.simple) %>%  as.data.frame()
dd$Pseudotime = d$Pseudotime
dd$timepoint = d$timepoint
index1 = sig[1]
index2 = sig[length(sig)]
d3 = dd %>% gather(gene, value, index1:index2 ) 
d3$cat = ifelse(d3$gene %in% cat1, '1', '2')


p1 = ggplot(data = d3 %>% 
              filter(Pseudotime > 5 & cat ==1 )) + 
  geom_smooth(data = d3 %>% 
                filter(Pseudotime > 5 & timepoint == 'd0' & cat ==1 ), 
              mapping = aes(x  = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(1), se = FALSE) + 
  geom_smooth(data = d3 %>% 
                filter(Pseudotime > 5 & timepoint == 'd1' & cat ==1 ),
              mapping =  aes(x = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(2)[2], 
              se = FALSE, size = 2) + 
  theme_bw() + 
  geom_vline(xintercept = 9.5, linetype = 'dashed') + 
  theme(strip.background = element_blank()) + 
  ylab('MTORC1 and Hypoxia\n Category 1 genes ') + 
  theme(axis.title = element_text(size = 18))
p1
ggsave(p1,filename = paste0(figpath, 'mtorcat1.pdf'), width = 3.7, height = 3.5)


p2 = ggplot(data = d3 %>% 
              filter(Pseudotime > 5 & cat ==2 )) + 
  geom_smooth(data = d3 %>% 
                filter(Pseudotime > 5 & timepoint == 'd0' & cat ==2 ), 
              mapping = aes(x  = Pseudotime, y = value, group = gene),
              color = ggsci::pal_jama(alpha = 0.2)(1), se = FALSE) + 
  geom_smooth(data = d3 %>% filter(Pseudotime > 5 & timepoint == 'd1' & cat ==2),
              mapping =  aes(x  = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(2)[2], 
              se = FALSE, size = 2) + 
  theme_bw() + 
  geom_vline(xintercept = 9.5, linetype = 'dashed') + 
  theme(strip.background = element_blank()) + 
  ylab('MTORC1 and Hypoxia\n Category 2 genes ') + 
  theme(axis.title = element_text(size = 18))
ggsave(p2,filename = paste0(figpath, 'mtorcat2.pdf'), width = 3.7, height = 3.5)
```

Pathway enrichment within the integrated pseudotime and perturbation
based gene categories defined above.  
note due to issues the enrichr R package servers, the html web based
enrichr server was used here. This script contains links to those
results and the genes used as input. This script uses the
LeadingEdgeIndexed function from scglmmr requiring R 4.0.5.
mid_res/monocyte_map/4_genecat_enrichr.r

``` r
# R 4.0.5
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
library(enrichR)
figpath = here("mid_res/monocyte_map/figures/"); dir.create(figpath, recursive = TRUE)
datapath = here("mid_res/monocyte_map/generated_data/"); dir.create(datapath, recursive = TRUE)


# monocyte leading edge d1 
g1c = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/g1c.rds'))
g1c = lapply(g1c, function(x) x %>%  filter(NES > 0))
mono.le = LeadingEdgeIndexed(gsea.result.list = g1c,padj.threshold = 0.05)


# branch ependent genes 
de.branch = readRDS(file = here('mid_res/monocyte_map/data/de_branch.rds'))
de.branch.sub = de.branch %>% filter(qval < 0.05)
branch.genes = as.character(de.branch.sub$gene_short_name)


# enrichr
dbs <- c("GO_Molecular_Function_2015",
         "GO_Cellular_Component_2015",
         "GO_Biological_Process_2015")

cat2.ifn =  c('IFITM2', 'PTPN1', 'EIF4E2', 'IFITM3', 'HLA-C')
cat1.ifn = intersect(mono.le$CD14_Mono$`reactome interferon signaling`, branch.genes)
cat1.ifn = setdiff(cat1.ifn, cat2.ifn)


# mtor 
cat2.mtor  =c('CTSC', 'PFKL', 'ACTR3', 'CITED2', 'PGK1', 'INSIG1', 'CHST2')
cat1.mtor = intersect(mono.le$CD14_Mono$`HALLMARK MTORC1 signaling`, branch.genes)
cat1.mtor = setdiff(cat1.mtor, cat2.mtor)

# Cat 1 Mtor 
# https://maayanlab.cloud/Enrichr/enrich?dataset=eb3028063e4833deb6212e24235dffee

# Cat 2 mtor 
# https://maayanlab.cloud/Enrichr/enrich?dataset=ace0e94a6d3c1fa0037be65c0f677aa2

# Cat 1 IFN 
# https://maayanlab.cloud/Enrichr/enrich?dataset=537c52fe1a6621033587e9048bf20e98

# Cat 2 ifn 
# https://maayanlab.cloud/Enrichr/enrich?dataset=e3ef3178aeb695f05f90ae855b80da21
```

### Fig 3. & FigS3 mixed effects timed vaccination response model – AS03 CITE-seq cohort <a name="fig3.1"></a>

Fit combined mixed model with unadjuvanted and adjuvanted subjects and
apply contrast to define difference in 24h fold change post vaccination
within protein subsets adjusted for age and sex.  
mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/1_V4_AS03_contrastmodel.r

``` r
# R version 4.0.5 
# H5 vs H1 day 1 DE contrast pseudobulk model 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(variancePartition))
suppressMessages(library(scglmmr))
source("functions/analysis_functions.R")

# make output directories 
datapath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/"); 
dir.create(datapath, recursive = TRUE)
figpath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/figuresV4/"); 
dir.create(figpath, recursive = TRUE)

# parallel options for dream lme4 fits 
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# read metadata to exract sample names in analysis 
samplemd = readRDS(file = here('data/samplemd.rds')) 
d1sx = samplemd %>% filter(time_cohort =='d1') %$% sample

# read processed pseudobulk data 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))

# remove cell type string from sample names and subset to day 1 cohort
cnames = gsub("~.*","",colnames(pb[[1]]))
pb = lapply(pb, function(x){ 
  x %>% as.data.frame() %>% 
    setNames(nm = cnames) %>% 
    select(all_of(d1sx)) %>% 
    as.matrix() 
  })

# sample metadata for contrast model 
samplemd =readRDS(file = here('data/samplemd.rds')) 
samplemd = 
  samplemd %>% 
  filter(! time_cohort == 'd7') %>% 
  mutate(scaledage = (age - mean(age)) / sd(age)) %>% 
  rename('subjectid' = sampleid) %>% 
  rename('group' = adjmfc.group) %>%  
  mutate(group = ifelse(group %in% c('high', 'low'), "NOAS03", "AS03")) %>% 
  mutate(time.group = paste(timepoint, group, sep = "_"))
# get rid of space
samplemd$time.group = str_replace_all(
  string = samplemd$time.group,
  pattern = ' ',
  replacement = ''
)
# re-level time.group into ordered combined factor 
samplemd$time.group = factor(
  samplemd$time.group,
  levels = c('d0_AS03', 'd1_AS03', 'd0_NOAS03', 'd1_NOAS03')
  )
# format metadata 
samplemd = samplemd %>% 
  remove_rownames() %>% 
  column_to_rownames('sample')

# designmat 
met = samplemd[ ,c('gender', 'scaledage', 'time.group')]
mat = model.matrix( ~ 0 + time.group +  gender + scaledage, data = met)

################################
# specify random intercept model formula with combined factor 
f1 <- ~ 0 + time.group + gender + scaledage + (1|subjectid) 

# specify contrast matrix to test the fold change difference 
# based on levels of time.group this should be cmat = c(-1, 1, 1, -1, 0, 0)
L2 = makeContrastsDream(formula = f1, data = samplemd,
  contrasts = c(delta = "(time.groupd1_AS03 - time.groupd0_AS03) - (time.groupd1_NOAS03 - time.groupd0_NOAS03)")
  )
plotContrasts(L2) + ggsave(filename = paste0(figpath,'contrastmodel.pdf'), width = 7, height = 4)

# fit model on each subset 
# init store 
fit1 = v1 = list()
for (i in 1:length(pb)) {
  
  # init data 
  meta = samplemd
  form = f1 
  contrast_matrix = L2
  counts = pb[[i]]
  
  # dge list 
  d = edgeR::DGEList(counts = counts, samples = meta)
  
  # filter cell type specific lowly expressed genes and calc norm factors 
  gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, design = mat)
  print(names(pb)[i]);print(table(gtable))
  d = d[gtable, keep.lib.sizes=FALSE]
  d = edgeR::calcNormFactors(object = d)
  
  
  # get voom observation level weights 
  v = voomWithDreamWeights(counts = d,
                           formula = form,
                           data = meta,
                           BPPARAM = pparam,
                           plot = TRUE, save.plot = TRUE)
  
  # fit contrast mixed model 
  fitmm = dream(exprObj = v, 
                formula = form,
                data = meta,
                L = contrast_matrix,
                BPPARAM = pparam, 
                useWeights = TRUE, REML = TRUE)
  
  # save results 
  v1[[i]] = v
  fit1[[i]] = fitmm
}
names(v1) = names(fit1) = names(pb)

# save day 1 contrast fit 
saveRDS(object = fit1, file = paste0(datapath, 'fit1.rds'))
saveRDS(object = v1, file = paste0(datapath, 'v1.rds'))


# apply eBayes to model fit 
fit1 = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/fit1.rds'))
fit1e  = lapply(fit1, variancePartition::eBayes) 
saveRDS(object = fit1e, file = paste0(datapath, 'fit1e.rds'))


# save model fitting data 
saveRDS(object = samplemd, file = paste0(datapath, 'samplemd.rds'))
saveRDS(object = pb, file = paste0(datapath, 'pb.rds'))
saveRDS(object = L2, file = paste0(datapath, 'L2.rds'))

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] viridis_0.5.1            viridisLite_0.3.0        scglmmr_0.1.0            variancePartition_1.25.6 BiocParallel_1.24.1      limma_3.46.0            
# [7] magrittr_2.0.1           SeuratObject_4.0.0       Seurat_4.0.1             here_1.0.1               forcats_0.5.1            stringr_1.4.0           
# [13] dplyr_1.0.4              purrr_0.3.4              readr_1.4.0              tidyr_1.1.2              tibble_3.0.6             ggplot2_3.3.3           
# [19] tidyverse_1.3.0         
# 
# loaded via a namespace (and not attached):
# [1] estimability_1.3            scattermore_0.7             coda_0.19-4                 bit64_4.0.5                 irlba_2.3.3                
# [6] multcomp_1.4-16             DelayedArray_0.16.3         data.table_1.14.0           rpart_4.1-15                RCurl_1.98-1.3             
# [11] doParallel_1.0.16           generics_0.1.0              snow_0.4-3                  BiocGenerics_0.36.1         RhpcBLASctl_0.21-247.1     
# [16] cowplot_1.1.1               TH.data_1.0-10              RSQLite_2.2.7               shadowtext_0.0.9            RANN_2.6.1                 
# [21] future_1.21.0               bit_4.0.4                   enrichplot_1.10.2           spatstat.data_2.1-0         xml2_1.3.2                 
# [26] lubridate_1.7.9.2           httpuv_1.5.5                SummarizedExperiment_1.20.0 assertthat_0.2.1            hms_1.0.0                  
# [31] promises_1.2.0.1            progress_1.2.2              caTools_1.18.1              dbplyr_2.1.0                readxl_1.3.1               
# [36] igraph_1.2.6                DBI_1.1.1                   htmlwidgets_1.5.3           spatstat.geom_2.0-1         stats4_4.0.5               
# [41] ellipsis_0.3.1              ggpubr_0.4.0                backports_1.2.1             annotate_1.68.0             aod_1.3.1                  
# [46] deldir_0.2-10               MatrixGenerics_1.2.1        vctrs_0.3.6                 Biobase_2.50.0              ROCR_1.0-11                
# [51] abind_1.4-5                 cachem_1.0.4                withr_2.4.1                 ggforce_0.3.3               emmeans_1.5.4              
# [56] sctransform_0.3.2           prettyunits_1.1.1           goftest_1.2-2               cluster_2.1.2               DOSE_3.16.0                
# [61] lazyeval_0.2.2              crayon_1.4.1                labeling_0.4.2              edgeR_3.32.1                pkgconfig_2.0.3            
# [66] tweenr_1.0.2                GenomeInfoDb_1.26.7         nlme_3.1-152                rlang_0.4.10                globals_0.14.0             
# [71] lifecycle_1.0.0             miniUI_0.1.1.1              sandwich_3.0-0              downloader_0.4              modelr_0.1.8               
# [76] cellranger_1.1.0            rprojroot_2.0.2             polyclip_1.10-0             GSVA_1.38.2                 matrixStats_0.58.0         
# [81] lmtest_0.9-38               graph_1.68.0                Matrix_1.3-2                carData_3.0-4               boot_1.3-27                
# [86] zoo_1.8-8                   reprex_1.0.0                pheatmap_1.0.12             ggridges_0.5.3              png_0.1-7                  
# [91] bitops_1.0-6                KernSmooth_2.23-18          blob_1.2.1                  qvalue_2.22.0               parallelly_1.23.0          
# [96] rstatix_0.7.0               S4Vectors_0.28.1            ggsignif_0.6.0              scales_1.1.1                memoise_2.0.0              
# [101] GSEABase_1.52.1             plyr_1.8.6                  ica_1.0-2                   gplots_3.1.1                zlibbioc_1.36.0            
# [106] compiler_4.0.5              scatterpie_0.1.7            RColorBrewer_1.1-2          lme4_1.1-26                 fitdistrplus_1.1-3         
# [111] cli_2.5.0                   XVector_0.30.0              lmerTest_3.1-3              listenv_0.8.0               patchwork_1.1.1            
# [116] pbapply_1.4-3               MASS_7.3-53.1               mgcv_1.8-34                 tidyselect_1.1.0            stringi_1.5.3              
# [121] GOSemSim_2.16.1             locfit_1.5-9.4              ggrepel_0.9.1               grid_4.0.5                  fastmatch_1.1-0            
# [126] tools_4.0.5                 rio_0.5.16                  future.apply_1.7.0          parallel_4.0.5              rstudioapi_0.13            
# [131] foreign_0.8-81              foreach_1.5.1               gridExtra_2.3               farver_2.0.3                Rtsne_0.15                 
# [136] ggraph_2.0.5                digest_0.6.27               rvcheck_0.1.8               BiocManager_1.30.10         shiny_1.6.0                
# [141] Rcpp_1.0.6                  car_3.0-10                  GenomicRanges_1.42.0        broom_0.7.5                 egg_0.4.5                  
# [146] later_1.1.0.1               RcppAnnoy_0.0.18            org.Hs.eg.db_3.12.0         httr_1.4.2                  AnnotationDbi_1.52.0       
# [151] Rdpack_2.1.1                colorspace_2.0-0            rvest_0.3.6                 XML_3.99-0.6                fs_1.5.0                   
# [156] tensor_1.5                  reticulate_1.18             IRanges_2.24.1              splines_4.0.5               uwot_0.1.10                
# [161] statmod_1.4.35              spatstat.utils_2.1-0        graphlayouts_0.7.2          plotly_4.9.3                xtable_1.8-4               
# [166] jsonlite_1.7.2              nloptr_1.2.2.2              tidygraph_1.2.0             ggfun_0.0.4                 R6_2.5.0                   
# [171] pillar_1.4.7                htmltools_0.5.1.1           mime_0.10                   glue_1.4.2                  fastmap_1.1.0              
# [176] minqa_1.2.4                 clusterProfiler_3.18.1      codetools_0.2-18            fgsea_1.16.0                mvtnorm_1.1-1              
# [181] lattice_0.20-41             spatstat.sparse_2.0-0       numDeriv_2016.8-1.1         pbkrtest_0.5-0.1            curl_4.3                   
# [186] leiden_0.3.7                gtools_3.8.2                zip_2.1.1                   openxlsx_4.2.3              GO.db_3.12.1               
# [191] survival_3.2-10             munsell_0.5.0               DO.db_2.9                   GenomeInfoDbData_1.2.4      iterators_1.0.13           
# [196] haven_2.3.1                 reshape2_1.4.4              gtable_0.3.0                rbibutils_2.0               spatstat.core_2.0-0
```

Gene set enrichment and defining As03 specific cell type specific
leading edge phenotypes.  
mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/2_V4_AS03_contrastmodel_enrichment.r

``` r
# R version 4.0.5 
# H5 vs H1 day 1 DE contrast pseudobulk model Enrichment
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(scglmmr))
suppressMessages(library(variancePartition))
suppressMessages(library(magrittr))

# specify output directories and init parallel opts
datapath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/")
dir.create(datapath)
figpath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/figuresV4/gsea/")
dir.create(figpath, recursive = TRUE)

# parallel options 
BiocParallel::register(BiocParallel::SnowParam(4))
pparam = BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# load combined modules -- rm the RP gene outliers 
cmod = readRDS(file = here('signature_curation/combined_sig_sub.rds'))

# load contrast fit results 
fit1e = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/fit1e.rds'))
toprank = ExtractResult(
  model.fit.list = fit1e,
  what = 'lmer.z.ranks',
  coefficient.number = 1,
  coef.name = 'delta'
)

# gsea on combined modules
gc = FgseaList(rank.list.celltype = toprank,pathways = cmod, BPPARAM = pparam)
saveRDS(gc,file = paste0(datapath, 'gc.rds'))
gc = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/gc.rds'))

# index leading edge genes 
li = LeadingEdgeIndexed(gsea.result.list = gc,padj.threshold = 0.1)
li = Filter(li,f =  length)
saveRDS(li,file = paste0(datapath, 'li.rds'))

# jaccard enrichment
gsub = lapply(gc, function(x) x %>%  filter(padj < 0.1))
gsub = Filter(gsub, f = nrow)
ji = EnrichmentJaccard(gsealist = gsub, 
                       indexedgenes = li, 
                       saveplot = FALSE, 
                       figpath = figpath, 
                       returnJaccardMtx = TRUE)
saveRDS(ji$sortedgsea,file = paste0(datapath, 'sortedgsea.rds'))

# subet of leading edge genes upregulated in contrast model only. 
g.up = lapply(gc, function(x)
  x %>%  filter(padj < 0.1 & NES > 0)
  )
li.up = LeadingEdgeIndexed(gsea.result.list = g.up, padj.threshold = 0.1)
li.up = Filter(li.up, f = length)
saveRDS(li.up,file = paste0(datapath, 'li.up.rds'))

# add top gnees not included in gsea pathways 
res = ExtractResult(model.fit.list = fit1e, coefficient.number = 1, coef.name = 'delta')
topgene = lapply(res, function(x)
  x %>% 
    filter(logFC > 0.25 & P.Value < 0.03) %$% 
    gene
  )
saveRDS(topgene, file = paste0(datapath, 'topgene.rds'))

# 20% fdr subset jaccard enrichment 
li2 = LeadingEdgeIndexed(gsea.result.list = gc,padj.threshold = 0.2)
saveRDS(li2,file = paste0(datapath,'li2.rds'))

gc.up = lapply(gc, function(x) x %>% filter(padj < 0.2 & NES > 0))
li2.up = LeadingEdgeIndexed(gsea.result.list = gc.up,padj.threshold = 0.2)
li2.up = Filter(li2.up, f = length)
saveRDS(li2.up, file = paste0(datapath,'li2.up.rds'))
```

Calculate log cpm for gene distribution visualization in next script.  
mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/3_V4_calc_logcpm_tidyaveragedata.r

``` r
# average gene distributions - part 1 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
source(here('functions/MattPMutils.r'))

# calculate logcpm 
pb = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/pb.rds'))
samplemd = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/samplemd.rds'))
time.group = factor(samplemd$time.group)
av = list()
for (i in 1:length(pb)) {
  dge =  edgeR::DGEList( pb[[i]] ) 
  gtable = edgeR::filterByExpr(y = dge$counts, min.count = 3,  design = time.group)
  dge = dge[gtable, ]
  av[[i]]  = edgeR::cpm(dge, log = TRUE)
}
names(av) = names(pb)

# create tidy data for gene visualization and calculation of average signature scores 
# get tidy summary data 
av_tidy = list()
for (i in 1:length(av)) {
  ct = names(av)[i]
  gs = rownames(av[[i]])
  av_tidy[[i]] = GetTidySummary(av.exprs.list = av, 
                                celltype.index = i,
                                genes.use = gs)  %>% 
    mutate(cohort = if_else( str_sub(sample, 1,2) == "H5", "H5N1", "H1N1")) %>% 
    mutate(group = paste(str_sub(sample, -2,-1), cohort))
  av_tidy[[i]]$group= factor(av_tidy[[i]]$group, 
                             levels = 
                               c("d0 H1N1" ,"d1 H1N1" ,"d0 H5N1","d1 H5N1"))
  av_tidy[[i]]$group = plyr::revalue(av_tidy[[i]]$group, 
                                     c("d0 H1N1" = "d0 No-AS03",
                                       "d1 H1N1" = "d1 No-AS03",
                                       "d0 H5N1" = "d0 AS03",
                                       "d1 H5N1" = "d1 AS03"))
}
names(av_tidy) = names(av)
saveRDS(av_tidy, file = paste0(datapath,'av_tidy.rds'))
```

Visualize gene distributions of AS03 specific perturbation effects
across time between unadjuvanted and AS03 adjuvanted cohorts.  
mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/5_V4_AS03model_gene_distributions.r

``` r
# R version 4.0.5 
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(scglmmr))
source(here('functions/MattPMutils.r'))

# save path 
datapath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/as03fig/")
figpath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/figuresV4/as03fig/")

# load aggregated tidy data 
av_tidy = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gene_dist/av_tidy.rds'))

# load mixed model fit res v4
fit1e = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/fit1e.rds'))
fitres = ExtractResult(model.fit.list = fit1e, coefficient.number = 1,coef.name = 'delta')

# load gsea results 
gc = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/gc.rds'))

# load leading edge indexed 
li = readRDS(file = 'mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li.rds')
li2 = readRDS(file = 'mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li2.rds')

# load top genes 
topgene = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/topgene.rds'))


# specify specs
cu = c("grey48", "grey", "grey48",  "deepskyblue3")
cu.alpha = sapply(cu, col.alpha, alpha = 0.4) %>% unname()

mtheme = list(
  theme_bw(base_size = 10.5), 
  theme(axis.title.x = element_blank()), 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)), 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold",family = "Helvetica"), 
        axis.text.y =  element_text(size = 6),
        # new line 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = 'black')
  )
)
box_gg = list(
  ylab("log CPM"),
  geom_boxplot(show.legend = FALSE, outlier.shape = NA), 
  mtheme 
)


# CD14 Monocytes 
for (i in 1:length(li$CD14_Mono)) {
  names(li$CD14_Mono)[i] %>%  print()
  dplyr::intersect(topgene$CD14_Mono, unlist(li$CD14_Mono[i],use.names = FALSE) ) %>% print()
}

# Format monocyte subset plot 
gene_highlight = c('MB21D1', 'FPR2', 'P2RY13', 'TLR4')
mplt = av_tidy$CD14_Mono %>% filter(gene %in% gene_highlight)
mplt$gene[mplt$gene == "MB21D1"] = "CGAS"
mplt$gene = factor(mplt$gene, levels = c("CGAS",'FPR2', 'P2RY13', 'TLR4'))

# monocyte subset 
p = ggplot(mplt, aes(x =group, y = count, fill = group , color = group)) +
  mtheme +
  box_gg + 
  facet_wrap(~gene, scales = "free", nrow = 1) + 
  scale_color_manual(values = cu) +
  scale_fill_manual(values = cu.alpha) 
p
ggsave(p, filename = paste0(figpath, 'mono.subset.pdf'), width = 3.5, height = 2)


# mDC 
for (i in 1:length(li2$mDC)) {
  names(li2$mDC)[i] %>%  print()
  dplyr::intersect(topgene$mDC, unlist(li2$mDC[i],use.names = FALSE) ) %>% print()
}

gene_highlight = c('FPR1', 'CCR1', 'P2RY13', 'TLR4')
mplt = av_tidy$mDC %>% filter(gene %in% gene_highlight)
mplt$gene = factor(mplt$gene , levels =c('FPR1', 'CCR1', 'P2RY13', 'TLR4'))
# mDC subset 
p = ggplot(mplt, aes(x =group, y = count, fill = group , color = group)) +
  mtheme +
  box_gg + 
  facet_wrap(~gene, scales = "free", nrow = 1) + 
  scale_color_manual(values = cu) +
  scale_fill_manual(values = cu.alpha) 
p
ggsave(p, filename = paste0(figpath, 'mdc.subset.pdf'), width = 3.5, height = 2)


# B cell Naive
for (i in 1:length(li$BC_Naive)) {
  names(li$BC_Naive)[i] %>%  print()
  dplyr::intersect(topgene$BC_Naive, unlist(li$BC_Naive[i],use.names = FALSE) ) %>%
    print()
}
bplt = av_tidy$BC_Naive %>% filter(gene == 'PMAIP1')
bplt$gene[bplt$gene == "PMAIP1"] = "NOXA (PMAIP1)"

# B cell subset 
p = ggplot(bplt, aes(x =group, y = count, fill = group , color = group)) +
  mtheme +
  box_gg + 
  facet_wrap(~gene, scales = "free", nrow = 1) + 
  scale_color_manual(values = cu) +
  scale_fill_manual(values = cu.alpha) 
p
ggsave(p, filename = paste0(figpath, 'Bcell.PMAIP1.pdf'), width = 1.4, height = 2.5)
```

### Fig 3. & FigS3 B cell AS03 phenotype analysis <a name="fig3.2"></a>

Further analysis of B cell phenotypes. Single cell model of apoptosis
signature and correlation of apoptosis signature with B CD40 Activation
signature.  
mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/6_V4_bsignal.r

``` r
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
```

Distribution of B cell apoptosis genes from all fitted mixed model z
statistics from fold change contrast within naive B cells.  
mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/7_bgenes.r

``` r
# Save results table 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
source(here('functions/MattPMutils.r'))
datapath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/bsig_data/")
figpath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/figuresV4/bsig/")

# read gsea and mixed model results 
fit1e = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/fit1e.rds'))
res = ExtractResult(model.fit.list = fit1e, coefficient.number = 1, coef.name = 'delta')
gc = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/gc.rds'))
m160 = gc$BC_Naive$leadingEdge[[6]]
apop = c('PMAIP1 (NOXA)', 'BCL2', 'BTG2', 'BTG1')
res$BC_Naive$gene[res$BC_Naive$gene == 'PMAIP1'] <- 'PMAIP1 (NOXA)'
p = 
  ggplot(res$BC_Naive, aes(x = logFC, y = z.std)) + 
  theme_bw() + 
  geom_bin2d(bins = 400, show.legend = FALSE, fill = 'black') +
  ylab('Mixed model contrast \n standardized z statistic') + 
  xlab('Difference in day 1 log fold change\nAS03 vs unadjuvanted') +
  geom_point(data = res$BC_Naive %>%  filter(gene %in% apop | gene %in% m160), aes(x = logFC, y = z.std),
             color = 'deepskyblue3', size = 2) +
  ggrepel::geom_text_repel(data = res$BC_Naive %>%  filter(gene %in% apop), 
                           aes(x = logFC, y = z.std, label = gene),
                           color = 'red',size = 5,
                           nudge_x = -0.2,
                           nudge_y = -1, box.padding = 1,
                           segment.size = 0.1) + 
  theme(title = element_text(size = 18)) + 
  ggtitle('Naive B cells')
p
ggsave(p,filename = paste0(figpath,'bcn_contrast_genes.pdf'), width = 5, height = 5)
```

Analysis of surface plasmon resonance data to correlate strain specific
with non strain antibody avidity in AS03 adjuvanted donors.
mid_res/ru/RU_binding.r

``` r
suppressMessages(library(tidyverse))
suppressMessages(library(here))

figpath =here("mid_res/ru/figures/"); dir.create(figpath)
datapath =here("mid_res/ru/generated_data/"); dir.create(datapath, recursive = TRUE)


# define adjuvant subjects 
met = data.table::fread(here("data/CHI_H5N1_data/clinical_info_adj.txt"))

# define adjuvant subjects 
adj.subjects = met %>%
  select(`Subject ID`, Adjuvant) %>%
  filter(Adjuvant == 'Adj') %>% 
  select(subjectid = `Subject ID`, Adjuvant) 

# load SPR data 
ru = read_delim(file = "data/CHI_H5N1_data/MN_abbinding/MPMEDIT_CHI_H5N1_AS03_SPR_data_2017_Khurana_SK.txt",delim = '\t')
CITE = c("H5N1-011", "H5N1-017", "H5N1-021", "H5N1-031", "H5N1-038", "H5N1-043")
ru.sub = 
  ru %>%
  separate(Sera,into = c('sx','time'),sep = '-') %>% 
  filter(time == ' D42') %>% 
  filter(subjectid %in% adj.subjects$subjectid) %>% 
  mutate(cite = ifelse(subjectid %in% CITE, '1', '0')) %>% 
saveRDS(ru.sub,file = paste0(datapath,'ru.sub.rds'))


p = 
  ggplot(data  = ru.sub, aes(x = Indo_RU_HA1, y = Viet_RU_HA1 )) + 
  theme_bw() +
  geom_smooth(method = "lm", color = "black") + 
  geom_point(shape = 21 , size = 3, fill = "deepskyblue3") + 
  theme(legend.position = "top") + 
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10, color = 'black')) +
  # ggrepel::geom_text_repel(data = ru.sub %>% filter( cite == 1), size = 2.8, aes(label = subjectid)) +
  xlab("Day 42 Antibody Binding (RU) \n Heterologous H5N1 strain (Vietnam)") + 
  ylab("Day 42 Antibody Binding (RU) \n vaccine H5N1 strain (Indonesia)") + 
  ggpubr::stat_cor(method = "pearson")
ggsave(p, filename = paste0(figpath,"RU_plot.pdf" ), width = 3.5, height = 3.5)
```

### Fig 3. & FigS3 mixed effects timed vaccination response model – AS03 Validatiton cohort <a name="fig3.3"></a>

Fit mixed model with unadjuvanted and AS03 adjuvanted subjects, apply
contrast to define difference in 24h fold change post vaccination on
FACS sorted subsets from validation cohort. Format data:  
mid_res/vand/1_V4_vand_d1_adjvsnon_DATAFORMAT.r

``` r
# R version 4.0.5 
# H5 vs H1 day 1 DE contrast pseudobulk model 
# reanalysis of Howard et. al AS03 vs PBS cotrol with sorted lineages rnaseq data
# data from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0167488#pone.0167488.s009
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

###### save paths  
datapath = here("mid_res/vand/generated_data/"); dir.create(datapath)

 # load log CPM data from supp table 2 list of celltypes 
vand_datapath = here("data/vand/data/")
ctd = list.files(path = vand_datapath, full.names = T)
ctd = ctd[-1]
e = lapply(ctd, function(x){read_delim(x, delim = ",")})

# get celltypes to name list 
celltypes = lapply(e,function(x){ str_sub(colnames(x), 3,5) %>% unique}) %>% unlist
celltypes = setdiff(celltypes, "SEM")
names(e) = celltypes

# setup data for limma / dream lme4
e = lapply(e, function(x){
  x = x %>%
    select(-ENSEMBL63_GENE_ID) %>% 
    select(matches('D000|D001|ENSEMBL63_GENE_NAME')) %>% 
    column_to_rownames("ENSEMBL63_GENE_NAME") 
})


# create metadata 
md = lapply(e, function(x){
  colnames(x) %>%
    as.data.frame() %>% 
    mutate(group = str_sub(., -4,-1)) %>% 
    mutate(group = if_else(group == "_PBS", "xPBS", group )) %>% 
    mutate(subjectid = str_sub(., 1,1)) %>% 
    mutate(timepoint = str_sub(., 7,10)) %>% 
    mutate(timepoint = ifelse(timepoint == 'D000', 'd0', 'd1')) %>% 
    mutate(time.group = paste(timepoint, group,sep = "_")) %>% 
    column_to_rownames(".")
})

# rm samples with missing data.
lapply(md, function(x) { x$subjectid})
### 
# md 1 = P 
# md 5 = P 
# md 6 = C J 

## remove missing data from sample data 
md[[1]] = md[[1]] %>% 
  rownames_to_column("sample") %>% 
  filter(!subjectid %in% "P") %>%
  column_to_rownames("sample")
md[[5]] = md[[5]] %>%
  rownames_to_column("sample") %>% 
  filter(!subjectid %in% "P") %>%
  column_to_rownames("sample")
md[[6]] = md[[6]] %>% 
  rownames_to_column("sample") %>% 
  filter(!subjectid %in% c("C", "J")) %>%
  column_to_rownames("sample")

#  remove missing data from RNAseq data 
e[[1]] = e[[1]][ ,rownames(md[[1]])]
e[[5]] = e[[5]][ ,rownames(md[[5]])]
e[[6]] = e[[6]][ ,rownames(md[[6]])]

#############################################
## create model matrix and check model rank 
d1m = lapply(md, function(x){
  x =  x %>% mutate_if(is.character, as.factor) %$% time.group ; 
  x = model.matrix(~0 + x) ; 
  colnames(x) = str_sub(colnames(x), start = -9, end = -1)
  return(x)
})

# re QC model 
for (i in 1:length(d1m)) {
  model = d1m[[i]] ; print(i)
  stopifnot(Matrix::rankMatrix(model) == ncol(model)) ; stopifnot(any(colSums(model) == 0) == FALSE)
}
stopifnot(all.equal(
  lapply(md, rownames), lapply(e,colnames)
))
  
#  Confirm data has been normalized: 
# lapply(e, boxplot)
saveRDS(md, file = paste0(datapath, 'md.rds'))
saveRDS(e, file = paste0(datapath, 'e.rds'))

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] magrittr_2.0.1  forcats_0.5.1   stringr_1.4.0   dplyr_1.0.4     purrr_0.3.4     readr_1.4.0     tidyr_1.1.2     tibble_3.0.6   
# [9] ggplot2_3.3.3   tidyverse_1.3.0 here_1.0.1     
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.6        cellranger_1.1.0  pillar_1.4.7      compiler_4.0.5    dbplyr_2.1.0      tools_4.0.5       lattice_0.20-41  
# [8] jsonlite_1.7.2    lubridate_1.7.9.2 lifecycle_1.0.0   gtable_0.3.0      pkgconfig_2.0.3   rlang_0.4.10      Matrix_1.3-2     
# [15] reprex_1.0.0      cli_2.5.0         rstudioapi_0.13   DBI_1.1.1         haven_2.3.1       withr_2.4.1       xml2_1.3.2       
# [22] httr_1.4.2        fs_1.5.0          generics_0.1.0    vctrs_0.3.6       hms_1.0.0         rprojroot_2.0.2   grid_4.0.5       
# [29] tidyselect_1.1.0  glue_1.4.2        R6_2.5.0          readxl_1.3.1      modelr_0.1.8      backports_1.2.1   scales_1.1.1     
# [36] ellipsis_0.3.1    rvest_0.3.6       assertthat_0.2.1  colorspace_2.0-0  stringi_1.5.3     munsell_0.5.0     broom_0.7.5      
# [43] crayon_1.4.1     
```

Fit model:  
mid_res/vand/2_V4_vand_d1_adjvsnon_MODEL.r

``` r
# R version 4.0.5 
# vanderilt validation cohort of contrast model results from CITE-seq 
# analysis within analagous main immune cell lineage. 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(variancePartition))
suppressMessages(library(scglmmr))

###### save paths  
datapath = here("mid_res/vand/generated_data/")

# parallel options 
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# load data 
e = readRDS(file = here('mid_res/vand/generated_data/e.rds'))
md = readRDS(file = here('mid_res/vand/generated_data/md.rds'))

# process metadata for contrast model 
# order time.group factor 
md_lmer = lapply(md, function(x){
  x = x %>% mutate(time.group = factor(time.group, 
      levels = c("d0_AS03", "d1_AS03", "d0_xPBS", "d1_xPBS")))
  })

# specify formula for lme4 / dream and set up contrasts 
f1 <- ~ 0 + time.group + (1|subjectid)  

# specify contrast matrix to test the fold change difference 
# based on levels of time.group this should be cmat = c(-1, 1, 1, -1)
L2 = makeContrastsDream(
  formula = f1,
  data = md_lmer[[1]],
  contrasts = c(delta = "(time.groupd1_AS03 - time.groupd0_AS03) - (time.groupd1_xPBS - time.groupd0_xPBS)")
)
plotContrasts(L2) + ggsave(filename = paste0(figpath,'contrastmodel.pdf'), width = 7, height = 4)


# fit model on each subset 
fit1 = fitne = list()
for (i in 1:length(e)) {
  
  # init data 
  norm_dat = e[[i]]
  meta = md_lmer[[i]]
  form = f1 
  contrast_matrix = L2

    # fit contrast mixed model on prenormalized values 
  fitmm = dream(exprObj = norm_dat, 
                formula = form,
                data = meta,
                L = contrast_matrix,
                BPPARAM = pparam, 
                useWeights = FALSE, 
                REML = TRUE)
  # save results 
  fitne[[i]] = fitmm
  fit1[[i]] = variancePartition::eBayes(fitmm) 
}
names(fit1) = names(fitne) = names(e)

# Save 
# day 1 contrast fit 
saveRDS(object = fit1, file = paste0(datapath, 'fit1.rds'))
saveRDS(object = fitne, file = paste0(datapath, 'fit1.rds'))
# model fitting data 
saveRDS(object = md_lmer, file = paste0(datapath, 'md_lmer.rds'))
saveRDS(object = L2, file = paste0(datapath, 'L2.rds'))

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scglmmr_0.1.0            variancePartition_1.25.6 BiocParallel_1.24.1      limma_3.46.0             magrittr_2.0.1          
# [6] here_1.0.1               forcats_0.5.1            stringr_1.4.0            dplyr_1.0.4              purrr_0.3.4             
# [11] readr_1.4.0              tidyr_1.1.2              tibble_3.0.6             ggplot2_3.3.3            tidyverse_1.3.0         
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.1.0            lme4_1.1-26                 RSQLite_2.2.7               AnnotationDbi_1.52.0       
# [5] grid_4.0.5                  scatterpie_0.1.7            munsell_0.5.0               codetools_0.2-18           
# [9] statmod_1.4.35              withr_2.4.1                 colorspace_2.0-0            GOSemSim_2.16.1            
# [13] Biobase_2.50.0              rstudioapi_0.13             stats4_4.0.5                ggsignif_0.6.0             
# [17] DOSE_3.16.0                 labeling_0.4.2              MatrixGenerics_1.2.1        Rdpack_2.1.1               
# [21] emmeans_1.5.4               GenomeInfoDbData_1.2.4      polyclip_1.10-0             pheatmap_1.0.12            
# [25] bit64_4.0.5                 farver_2.0.3                rprojroot_2.0.2             downloader_0.4             
# [29] coda_0.19-4                 vctrs_0.3.6                 generics_0.1.0              TH.data_1.0-10             
# [33] R6_2.5.0                    doParallel_1.0.16           GenomeInfoDb_1.26.7         graphlayouts_0.7.2         
# [37] locfit_1.5-9.4              bitops_1.0-6                cachem_1.0.4                fgsea_1.16.0               
# [41] DelayedArray_0.16.3         assertthat_0.2.1            scales_1.1.1                multcomp_1.4-16            
# [45] ggraph_2.0.5                enrichplot_1.10.2           gtable_0.3.0                egg_0.4.5                  
# [49] tidygraph_1.2.0             sandwich_3.0-0              rlang_0.4.10                splines_4.0.5              
# [53] rstatix_0.7.0               broom_0.7.5                 BiocManager_1.30.10         reshape2_1.4.4             
# [57] abind_1.4-5                 modelr_0.1.8                backports_1.2.1             qvalue_2.22.0              
# [61] clusterProfiler_3.18.1      tools_4.0.5                 ellipsis_0.3.1              gplots_3.1.1               
# [65] RColorBrewer_1.1-2          BiocGenerics_0.36.1         Rcpp_1.0.6                  plyr_1.8.6                 
# [69] progress_1.2.2              zlibbioc_1.36.0             RCurl_1.98-1.3              prettyunits_1.1.1          
# [73] ggpubr_0.4.0                viridis_0.5.1               cowplot_1.1.1               S4Vectors_0.28.1           
# [77] zoo_1.8-8                   SummarizedExperiment_1.20.0 haven_2.3.1                 ggrepel_0.9.1              
# [81] fs_1.5.0                    data.table_1.14.0           lmerTest_3.1-3              DO.db_2.9                  
# [85] openxlsx_4.2.3              reprex_1.0.0                mvtnorm_1.1-1               matrixStats_0.58.0         
# [89] hms_1.0.0                   GSVA_1.38.2                 xtable_1.8-4                pbkrtest_0.5-0.1           
# [93] RhpcBLASctl_0.21-247.1      XML_3.99-0.6                rio_0.5.16                  readxl_1.3.1               
# [97] IRanges_2.24.1              gridExtra_2.3               compiler_4.0.5              KernSmooth_2.23-18         
# [101] crayon_1.4.1                shadowtext_0.0.9            minqa_1.2.4                 ggfun_0.0.4                
# [105] snow_0.4-3                  lubridate_1.7.9.2           DBI_1.1.1                   tweenr_1.0.2               
# [109] dbplyr_2.1.0                MASS_7.3-53.1               boot_1.3-27                 Matrix_1.3-2               
# [113] car_3.0-10                  cli_2.5.0                   rbibutils_2.0               parallel_4.0.5             
# [117] igraph_1.2.6                GenomicRanges_1.42.0        pkgconfig_2.0.3             rvcheck_0.1.8              
# [121] numDeriv_2016.8-1.1         foreign_0.8-81              xml2_1.3.2                  foreach_1.5.1              
# [125] annotate_1.68.0             XVector_0.30.0              estimability_1.3            rvest_0.3.6                
# [129] digest_0.6.27               graph_1.68.0                cellranger_1.1.0            fastmatch_1.1-0            
# [133] edgeR_3.32.1                GSEABase_1.52.1             curl_4.3                    gtools_3.8.2               
# [137] nloptr_1.2.2.2              lifecycle_1.0.0             nlme_3.1-152                jsonlite_1.7.2             
# [141] aod_1.3.1                   carData_3.0-4               viridisLite_0.3.0           pillar_1.4.7               
# [145] lattice_0.20-41             fastmap_1.1.0               httr_1.4.2                  survival_3.2-10            
# [149] GO.db_3.12.1                glue_1.4.2                  zip_2.1.1                   iterators_1.0.13           
# [153] bit_4.0.4                   ggforce_0.3.3               stringi_1.5.3               blob_1.2.1                 
# [157] org.Hs.eg.db_3.12.0         caTools_1.18.1              memoise_2.0.0     
```

Test enrichment of CITE-seq derived AS03 specific cell phenotypes in the
validation cohort.  
mid_res/vand/3_V4_vand_enrCITE_sgnals.r

``` r
# gene set nenrichment signals in vand cohort 
# R 4.0.5 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(variancePartition))
suppressMessages(library(scglmmr))

###### save paths  
datapath = here("mid_res/vand/generated_data/")

# parallel opts
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# load signatures 
li.up = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li.up.rds'))
li2.up = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li2.up.rds'))
li = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li.rds'))
#li.full = unlist(li.up,use.names = FALSE, recursive = TRUE) %>%  unique 
mono.combined  = li.up$CD14_Mono %>%  unlist() %>% unique()
mdc.combined = li2.up$mDC %>% unlist() %>% unique()
nb.combined = li.up$BC_Naive %>% unlist() %>% unique()
t.combined = list(
  'Tcell.combined' = c(
    li.up$CD4_CD161_Mem_Tcell,
    li.up$CD4_CD25_Tcell,
    li.up$CD4_Efct_Mem_Tcell,
    li.up$CD4Naive_Tcell,
    li.up$CD8_Mem_Tcell,
    li.up$CD8_Naive_Tcell
  ) %>%
    unlist() %>%
    unique()
)
  

# add additional b cell signatures from apoptosis hypothesis 
sig.test = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/bsig_data/sig.test.rds'))

# add combined signals 
li.up$CD14_Mono$combined.signature = mono.combined
li2.up$mDC$combined.signature = mdc.combined
li$BC_Naive$combined.signature = nb.combined
li$BC_Naive = c(li$BC_Naive, sig.test)

# load vand fits and extract ranks 
fit1 = readRDS(file = here('mid_res/vand/generated_data/fit1.rds'))
vand.rank = ExtractResult(model.fit.list = fit1,
                          what = 'lmer.z.ranks', 
                          coefficient.number = 1, 
                          coef.name = 'delta')

# CD14 monocyte test in total sorted monocyte
mv = FgseaList(
  rank.list.celltype = list('MNC' = vand.rank$MNC),
  pathways = li.up$CD14_Mono,
  BPPARAM = pparam
)

# mDC test in sorted DC
dcv = FgseaList(
  rank.list.celltype = list('DNC' = vand.rank$DNC),
  pathways = li2.up$mDC,
  BPPARAM = pparam
)

# naive BC test in sorted total B
bcv = FgseaList(
  rank.list.celltype = list('BCL' = vand.rank$BCL),
  pathways = li$BC_Naive,
  BPPARAM = pparam
)


# T cell combined in sorted T celsl 
tcv = FgseaList(
  rank.list.celltype = list('TCL' = vand.rank$TCL),
  pathways = c(
    li.up$CD4_CD161_Mem_Tcell,
    li.up$CD4_CD25_Tcell,
    li.up$CD4_Efct_Mem_Tcell,
    li.up$CD4Naive_Tcell,
    li.up$CD8_Mem_Tcell,
    li.up$CD8_Naive_Tcell,
    t.combined
  ),
  BPPARAM = pparam
)

# save objects
saveRDS(object = mv,file = paste0(datapath, 'mv.rds'))
saveRDS(object = dcv,file = paste0(datapath, 'dcv.rds'))
saveRDS(object = bcv,file = paste0(datapath, 'bcv.rds'))
saveRDS(object = tcv,file = paste0(datapath, 'tcv.rds'))
```

### Fig 3. & FigS3 AS03 specific cell phenotypes validation comparison figures <a name="fig3.4"></a>

Combined figures for the CITE-seq and validation cohort.  
mid_res/combined_contrast/1_combined_contrast_vand_citeseq.r

``` r
# average gene distributions - part 2
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
source(here('functions/MattPMutils.r'))
# output directories 
figpath = here("mid_res/combined_contrast/figures/")
dir.create(figpath)

# load CITE results 
gc = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/gc.rds'))

#######
# CD14 monocytes cite-seq 
######
mo = gc$CD14_Mono %>%
  as.data.frame() %>% 
  filter(padj < 0.1 & NES > 0)
mo$pathway[mo$pathway == "REACTOME_PEPTIDE_LIGAND_BINDING_RECEPTORS" ] <- 'Reactome peptide ligand binding receptors'
mo$pathway[mo$pathway == "LI.M4.3 myeloid cell enriched receptors and transporters" ] <- 'LI.M4.3 myeloid receptors and transporters'
mo$pathway[mo$pathway == "LI.M37.0 immune activation - generic cluster" ] <- 'LI.M37.0 immune activation'

###### 
# mDCs cite-seq
######
dc = gc$mDC %>%
  as.data.frame() %>% 
  filter(padj<0.2 & NES > 0)
dc$pathway[dc$pathway == "REACTOME_PEPTIDE_LIGAND_BINDING_RECEPTORS" ] <- 'Reactome peptide ligand binding receptors'
dc$pathway[dc$pathway == "REACTOME_CLASS_A1_RHODOPSIN_LIKE_RECEPTORS" ] <- 'Reactome rhodopsin-like receptors'
dc$pathway[dc$pathway == "REACTOME_NFKB_AND_MAP_KINASES_ACTIVATION_MEDIATED_BY_TLR4_SIGNALING_REPERTOIRE" ] <- 'Reactome NFKB activation via TLR4'
dc$pathway[dc$pathway == "REACTOME_TAK1_ACTIVATES_NFKB_BY_PHOSPHORYLATION_AND_ACTIVATION_OF_IKKS_COMPLEX" ] <- 'Reactome TAK1 activates NFKB'
dc$pathway[dc$pathway == "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION" ] <- 'Reactome ligand receptor interaction'
dc$pathway[dc$pathway == "REACTOME_DNA_REPLICATION" ] <- 'Reactome DNA replication'
dc$pathway[dc$pathway == "REACTOME_ACTIVATED_TLR4_SIGNALLING" ] <- 'Reactome activated TLR4 signaling'
dc$pathway[dc$pathway == "REACTOME_G_ALPHA_I_SIGNALLING_EVENTS" ] <- 'Reactome G alpha signaling'


#####
# B cells cite-seq 
#####
bn = gc$BC_Naive %>%
  as.data.frame() %>% 
  filter(padj < 0.1) 

# add string for cohort 
mo$cohort = 'CITE-seq'
dc$cohort = 'CITE-seq'
bn$cohort = 'CITE-seq'

#################
# validation cohort
#################
mv = readRDS(file = here('mid_res/vand/generated_data/mv.rds')) 
dcv = readRDS(file = here('mid_res/vand/generated_data/dcv.rds')) 
bcv = readRDS(file = here('mid_res/vand/generated_data/bcv.rds')) 
mv = as.data.frame(mv$MNC)
dcv = as.data.frame(dcv$DNC)
bcv = as.data.frame(bcv$BCL)

# plot exactly as in the CITE-eq subset 
mv$pathway[mv$pathway == "REACTOME_PEPTIDE_LIGAND_BINDING_RECEPTORS" ] <- 'Reactome peptide ligand binding receptors'
mv$pathway[mv$pathway == "LI.M4.3 myeloid cell enriched receptors and transporters" ] <- 'LI.M4.3 myeloid receptors and transporters'
mv$pathway[mv$pathway == "LI.M37.0 immune activation - generic cluster" ] <- 'LI.M37.0 immune activation'


# DCS 
# plot exactly as in the CITE-eq subset 
dcv$pathway[dcv$pathway == "REACTOME_PEPTIDE_LIGAND_BINDING_RECEPTORS" ] <- 'Reactome peptide ligand binding receptors'
dcv$pathway[dcv$pathway == "REACTOME_CLASS_A1_RHODOPSIN_LIKE_RECEPTORS" ] <- 'Reactome rhodopsin-like receptors'
dcv$pathway[dcv$pathway == "REACTOME_NFKB_AND_MAP_KINASES_ACTIVATION_MEDIATED_BY_TLR4_SIGNALING_REPERTOIRE" ] <- 'Reactome NFKB activation via TLR4'
dcv$pathway[dcv$pathway == "REACTOME_TAK1_ACTIVATES_NFKB_BY_PHOSPHORYLATION_AND_ACTIVATION_OF_IKKS_COMPLEX" ] <- 'Reactome TAK1 activates NFKB'
dcv$pathway[dcv$pathway == "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION" ] <- 'Reactome ligand receptor interaction'
dcv$pathway[dcv$pathway == "REACTOME_DNA_REPLICATION" ] <- 'Reactome DNA replication'
dcv$pathway[dcv$pathway == "REACTOME_ACTIVATED_TLR4_SIGNALLING" ] <- 'Reactome activated TLR4 signaling'
dcv$pathway[dcv$pathway == "REACTOME_G_ALPHA_I_SIGNALLING_EVENTS" ] <- 'Reactome G alpha signaling'

# append with cohort 
dcv$cohort = 'validation'
dcv$celltype = 'sorted DC'

bcv$cohort = 'validation'
bcv$celltype = 'sorted B cells'

mv$cohort = 'validation'
mv$celltype = 'sorted monocytes'

# combine
col.keep = c('pathway', 'pval', 'padj', 'NES', 'celltype', 'cohort') 
r.list = list(dcv, bcv, mv, mo, dc, bn)
r.list = lapply(r.list, function(x) x %>% select(all_of(col.keep)))
d = bind_rows(r.list)

# group
d$main = ifelse(d$celltype %in% c('mDC', 'sorted DC'), yes = 'DC', no = d$celltype)
d$main = ifelse(d$celltype %in% c('CD14_Mono', 'sorted monocytes'), yes = 'Mono', no = d$main)
d$main = ifelse(d$celltype %in% c('BC_Naive', 'sorted B cells'), yes = 'BC', no = d$main)

d2 = d %>% filter(!celltype %in% c('BC_Naive', 'sorted B cells'))
d2$main = factor(d2$main, levels = c('Mono', 'DC'))


# add asterisk 
d2 = d2 %>%  filter(!pathway == 'combined.signature')
d3 = 
  d2 %>% 
  mutate(padj.validation = ifelse(cohort == 'validation', padj, no = Inf)) %>% 
  mutate(padj.citeseq = ifelse(cohort == 'CITE-seq', padj,no = Inf)) %>% 
  mutate(pathway.new = ifelse( padj.validation < 0.01, yes = paste0(' * ', pathway), no = pathway))
d3 %>% filter(cohort == 'validation')
d2$pathway = plyr::mapvalues(d2$pathway,from = d3$pathway,to = d3$pathway.new)

p = 
  ggplot(d2, 
       aes(x = NES, y = reorder(pathway, NES), size = -log10(padj), fill=cohort )) +
  xlim(c(-1,3)) +
  theme_bw() +
  geom_point(shape = 21) +
  scale_fill_manual(values = c(col.alpha('deepskyblue3', 0.5), col.alpha('#90C983',0.7))) + 
  ylab("") +
  xlab('Normalized Enrichment Score') + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  facet_grid(vars(main), scales = 'free', space = 'free') +
  theme_bw(base_size = 9) + 
  theme(axis.text = element_text(color = 'black')) +
  guides(fill = guide_legend(override.aes = list(size = 5))) 
p
ggsave(p, filename = paste0(figpath,'combined_as03_model.pdf'), width = 5, height = 3)
  


# B cells 
d3 = d %>% filter(celltype %in% c('BC_Naive', 'sorted B cells'))
d3$pathway = factor(d3$pathway, levels = c(
  "CD40_ACT", 
  "LI.S2 B cell surface signature",                     
  "LI.M47.0 enriched in B cells (I)",
  "LI.M69 enriched in B cells (VI)",                         
  "combined.signature", 
  "apoptosis.signature",                                     
  "LI.M160 leukocyte differentiation", 
  "LI.M165 enriched in activated dendritic cells (II)",    
  "LI.M43.0 myeloid, dendritic cell activation via NFkB (I)"
))

p =
  ggplot(d3 %>%filter(!pathway == 'combined.signature') %>%filter(cohort == 'validation'),
    aes(x = NES,y = reorder(pathway, NES),size = -log10(padj),fill = cohort)) +
  theme_bw() +
  geom_point(shape = 21) +
  scale_fill_manual(values = col.alpha('#90C983',0.7))  +
  ylab("") +
  xlab('Normalized Enrichment Score') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(color = 'black')) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme(legend.key.height = unit(0.3, units = 'cm'))
p
ggsave(p, filename = paste0(figpath,'validation_bc_as03_model.pdf'), width = 5, height = 2.3)
```

### Fig 4. Define high responder baseline cell phenotypes from multivariate model with enrichment <a name="fig4.1"></a>

Use model output from 1_h1_mixed_effect_workflow_V4.r to rank genes
based moderated t test statistics of high vs low responder effect at
baseline adjusted for age sex & batch.  
mid_res/baseline_response/1_baseline_gseaV3.r

``` r
# R version 4.0.5
# H1N1 Unadjuvanted group baseline gene set enrichment analysis
# gsea based on the ranks of the contrast high vs low responder pre vaccinaiton
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))

# set save paths 
datapath = here("mid_res/baseline_response/dataV3/"); dir.create(datapath)
figpath = here("mid_res/baseline_response/figuresV3/"); dir.create(figpath)

# parallel options for FseaList
BiocParallel::register(BiocParallel::SnowParam(4))
pparam = BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# load combined pathways 
mods = readRDS(file = here('signature_curation/combined_signatures.rds'))

# load baseline contrast, rank genes run gsea 
cont0 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV3/cont0.rds'))
r0 = ExtractResult(model.fit.list = cont0, what = 'gene.t.ranks',coefficient.number = 1, coef.name = 'adjmfc')

# run fgea on each cell type
g0 = FgseaList(rank.list.celltype = r0, pathways = mods, BPPARAM = pparam)
saveRDS(object = g0, file = paste0(datapath, 'g0.rds'))
```

These results read in a set of shortened module / pathway names for
visualization that is in the starting data folder Curate enrichment
results pt 1. mid_res/baseline_response/2_curate_gseaV3.r

``` r
# R version 4.0.5
# Curate H1N1 Unadjuvanted group baseline gene set enrichment analysis
# calculate pairwise jaccard index and reduce enrichments to major signals with low mutual information
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))

# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")

# read baseline enrichemnt results 
mrm = readRDS(file = here('signature_curation/module_rmlist.rds'))
g0 = readRDS(file = here('mid_res/baseline_response/dataV3/g0.rds'))
g0 = lapply(g0, function(x) x %>% filter(!pathway %in% mrm))
filtered_g0 = lapply(g0, function(x) x %>% filter(padj < 0.05))

# compute jaccard index of leadingedge genes within celltype  
li = LeadingEdgeIndexed(gsea.result.list = filtered_g0, padj.threshold = 0.05)

jres = EnrichmentJaccard(gsealist = filtered_g0, indexedgenes = li, 
                         saveplot = FALSE,
                         figpath = figpath,
                         returnJaccardMtx = TRUE, 
                         fontsize_row = 7.5, fontsize_col = 7.5)
d = jres$sortedgsea %>% 
  mutate(leadingEdge = map_chr(leadingEdge, toString)) %>% 
  select(celltype, av_jaccard,everything())
write_delim(d,file = paste0(datapath, 'g0jaccard.csv'),delim = ',')

# save the jaccard matrices 
jmats = jres$jaccard_matrix_list
saveRDS(jmats ,file = paste0(datapath, 'jmats.rds'))

sessionInfo()
```

Curate enrichment pt 2  
mid_res/baseline_response/3_gsea.vis.r

``` r
# R version 4.0.5
# visualization of gene set enrichment results. 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))
# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")

# res text gsea curated 
d = data.table::fread(here('mid_res/baseline_response/dataV3/g0jaccard.curated.txt')) %>% 
  filter(include ==1) %>% 
  mutate(signal = paste(celltype, pathway, sep = '~')) 

# gsea res raw 
g0 = readRDS(file = here('mid_res/baseline_response/dataV3/g0.rds'))

# filter to the gene sets from curated results. 
g0.sub = list()
for (i in 1:length(g0)) {
  g0.sub[[i]] =
    g0[[i]] %>% 
    mutate(signal = paste(celltype, pathway, sep = '~')) %>% 
    filter(signal %in% d$signal) %>% 
    mutate(celltype = str_replace_all(celltype,pattern = '_',replacement = ' '))
}
names(g0.sub) = names(g0)
# plot
p = PlotFgsea(gsea_result_list = g0.sub, padj_filter = 0.01)
ggsave(p,filename = paste0(figpath, 'gsea.g0sub.baseline.pdf'), width = 9.5, height = 6)
# save object test curated fgsea formatted results. 
saveRDS(g0.sub, file = paste0(datapath, 'g0.sub.rds'))
```

Within each subset calculate log cpm of cell type specific leading edge
enrichment phenotypes based on high vs low responder model.  
mid_res/baseline_response/4_baseline_exprs_amz_score.r

``` r
# create baseline expression leading edge module correlation matrix 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))
# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")

# read pb data, subset to day 0 non adj, subset out day 0 metadata. 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))
pb = lapply(pb, function(x) x = x[ ,1:40])
cnames = gsub("~.*","",colnames(pb[[1]]))
pb = lapply(pb, function(x){
  x %>% as.data.frame() %>% setNames(nm = cnames) %>% as.matrix() 
  })
d0 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV3/d0.rds'))
d0d = lapply(pb, function(x){ x = x[ , rownames(d0)]})


# convert pb data to log counts per million
d.norm = list()
for (i in 1:length(d0d)) {
  d = edgeR::DGEList(counts = d0d[[i]], samples = d0)
  gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, 
                               design = as.factor(d$samples$group))
  d = d[gtable, keep.lib.sizes=FALSE]
  d.norm[[i]] = edgeR::cpm(y = d, log = TRUE, prior.count = 1)
}
names(d.norm) = names(d0d)

# get leading edge genes from cur. baseline mods 
g0.sub = readRDS(file = here('mid_res/baseline_response/dataV3/g0.sub.rds'))
li.g0 = LeadingEdgeIndexed(gsea.result.list = g0.sub, padj.threshold = 0.05)
li.g0 = base::Filter(length, li.g0)

# subset normalized expression to subsets with baseline enrichments 
d.norm = d.norm[names(li.g0)]

res = list()
for (i in 1:length(d.norm)) {
  stopifnot(all.equal( names(d.norm[i]), names(li.g0[i]) ))
  zscore = scglmmr::calc_avg_module_zscore(
    module.list = li.g0[[i]], average.data.frame = d.norm[[i]]
  )
  rownames(zscore) = paste(rownames(zscore), names(d.norm[i]), sep = '~')
  res[[i]] = zscore
}

ds = do.call(rbind, res) %>% t()
saveRDS(ds, file = paste0(datapath, 'ds.rds'))

sessionInfo()
```

### Fig 4. Correlate expression of baseline high responder phenotypes with plasmablast response <a name="fig4.2"></a>

For each subject calculate the day 7 fold change of the predictive
antibody response signature (array data).
mid_res/baseline_response/4b_d7FC_response_sig.r

``` r
# baseline nettwork correlation 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))

# Day 7 array data and modules in day 7 data  # load day 7 gene sets 
sig7 = readRDS("signature_curation/core_d7.rds")
module.list = list("CHI_Day7_Response" = sig7$`CHI d7 Response`)

datapath = here("mid_res/baseline_response/dataV3/")

# read CHI array data 
subjects = data.table::fread(here('data/full_metadata/full_sample_metadata.txt')) %>% 
  filter(CITEdata == 1 & vaccine_cohort == 'H1N1') %$% 
  subjectid

array7 = 
  data.table::fread(here("data/CHI_H1N1_data/microarray/CHI_GE_matrix_gene.txt"), data.table = F) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("gene") %>% 
  select(which(substr(names(.),1,3) %in% subjects)) %>% 
  select(., matches("day0|day1|day7")) %>% 
  select(-matches(c("day70|pre|day1") )) %>% 
  select(-c('207_day0', '209_day0', '279_day0'))

# calculate microarray fold changge (data is already in log space)
t0 =  array7[ ,seq(from=1, to = ncol(array7), by = 2)]
t1 = array7[ ,seq(from=2, to = ncol(array7), by = 2)]
stopifnot(str_sub(colnames(t1), 1,3) == str_sub(colnames(t0), 1,3))
stopifnot(dim(t0) == dim(t1))
fc7 = t1 - t0
fc7 = as.data.frame(fc7)

# calculate the average z score of the day 7 fold change values across samples
d7res = calc_avg_module_zscore(module.list = module.list, average.data.frame = fc7)
names(d7res) = str_sub(names(d7res), 1,3)

# save results 
saveRDS(d7res, file = paste0(datapath, 'd7res.rds'))

sessionInfo()
```

Correlate expression of cell type specific baseline states with the day7
response signature.  
mid_res/baseline_response/4c_d7sigFC_vs_baseline_correlation.r

``` r
# correlate module expression avz with d7 FC in antibody predictive signature. 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))
source(here('functions/MattPMutils.r'))
# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/d7cor/");
dir.create(figpath, recursive = TRUE)

# day 7 response signature fc 
d7res = readRDS(file = here('mid_res/baseline_response/dataV3/d7res.rds'))

# baseline expression correlation 
ds = readRDS(file = here('mid_res/baseline_response/dataV3/ds.rds'))

# created shorter names 
new.names = data.table::fread(
  file = here('mid_res/baseline_response/dataV3/baseline.module.name.shortened.txt'),
  sep = '\t'
)
new.names$cname2 = paste(new.names$celltype, new.names$shortname, sep = ' :: ')
colnames(ds) = new.names$cname2

# format 
d7form = ds %>%
  as.data.frame() %>% 
  rownames_to_column('subject') %>% 
  mutate(subject = str_sub(subject, 1,3)) %>% 
  filter(subject %in% colnames(d7res)) %>% 
  column_to_rownames('subject')
saveRDS(d7form, file = paste0(datapath,'d7form.rds'))
d7form = readRDS(file = here('mid_res/baseline_response/dataV3/d7form.rds'))

# pairwise correlation with d7 response 
dd = cbind(t(d7res), as.data.frame(d7form))
saveRDS(dd, file = paste0(datapath,'dd.rds'))
d7.cor = Hmisc::rcorr(as.matrix(dd),type = 'spearman')
saveRDS(d7.cor, file = paste0(datapath,'d7.cor.rds'))
d7.cor = readRDS(file = here('mid_res/baseline_response/dataV3/d7.cor.rds'))

# aes set 
plotattr = list(
  theme_bw(),
  geom_point(shape = 21, size = 3.5, stroke = 0.8), 
  #geom_text(nudge_y = 0.05, size = 3), 
  ylab('Day 7 FC Antibody response signature'),
  theme(axis.title.y = element_text(size = 11)), 
  theme(axis.title.x = element_text(size = 8)), 
  scale_fill_manual(values = c(col.alpha("red",0.7), col.alpha("dodgerblue", 0.7))),
  theme(legend.position = 'none', legend.key.size = unit(0.29, units = 'cm')), 
  theme(aspect.ratio = 1) 
)

# specify response in vis.  
high.responders = c("205","207","209","212","215","234","237","245","250","256")

# calculate and vis correlation between day7 response signature fold change (bulk)
# versus log cpm of the day 7 signatures associated with adjMFC group 
for (i in 1:length(colnames(d7form))) {
  #i = 1 
  mod.names = c(colnames(d7form)[i], colnames(t(d7res)))
  cplot = cbind(as.data.frame(d7form[, i]), t(d7res))
  cplot$subject = rownames(d7form)
  cplot$response = ifelse(cplot$subject %in% high.responders, 'high', 'low')
  colnames(cplot)[1:2] = mod.names
  
  # -- for correlation -- 
  dsub = as.data.frame(cbind(v1 = cplot[ ,1], v2 = cplot[ ,2]))
  
  # calculate and vis. correlation across all subjects; color by response. 
  p = ggpubr::ggscatter(cplot, x = mod.names[1], y = mod.names[2],
                        color = col.alpha('white','0.01'),
                        add.params = list(color = "black", fill = "grey")
                        )  +
    plotattr +
    aes(fill = response) + 
    ggpubr::stat_cor(data = dsub, aes(x = v1, y = v2), method = 'spearman',
                     inherit.aes = FALSE, label.x.npc = "left",label.y.npc = "top", cor.coef.name =  "rho",) 
  # save name by cell type first 
  module = sub("^[^::]*::", "", mod.names[1])
  celltype = gsub("::.*", "", mod.names[1])
  ggsave(p, filename = paste0(figpath, celltype, module, 'd7cor.pdf'), width = 3.1, height = 3.1)
}
```

### Fig 4. & FigS4 Construct and visualize high responder multicellular network phenotypes<a name="fig4.3"></a>

Create baseline cell phenotype correlation network and ‘calculate shared
latent information’ for intracellular correlations.  
mid_res/baseline_response/5_baseline_sli_correlation_network.r

``` r
# baseline nettwork correlation 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
source('functions/MattPMutils.r')

# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")

################################
# Part I correct intracellular correlations for overlapping gene content
################################

# pairwise jaccard index of all leading edge genes from modules 
library(scglmmr)
g0.sub = readRDS(file = here("mid_res/baseline_response/dataV3/g0.sub.rds"))
li.index = scglmmr::LeadingEdgeIndexed(gsea.result.list = g0.sub,padj.threshold = 0.05)

# remove lists with no enrichments indexed by cell type 
g0.sub = base::Filter(g0.sub, f = nrow)
li.index = base::Filter(li.index, f = length)

# reorder enrichment to match leading edge index
g0.sub = g0.sub[names(li.index)]
stopifnot(isTRUE(all.equal(names(g0.sub), names(li.index))))

ss = li.index[c('CD14_Mono', 'CD16_Mono', 'mDC', 'MAIT_Like')]
ss = lapply(ss, function(x) x %>%  unlist() %>% unname() %>% unique())


VennDiagram::venn.diagram(ss, 
                          # color 
                          fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"), 
                          # text
                          cat.cex = 0.4,cat.default.pos = "inner",cat.fontfamily = "Helvetica",
                          # file 
                          imagetype="png", height = 1080, width = 1080, resolution = 300,compression = "lzw",
                          filename = paste0(figpath,"gpvenn.png")
)
go = VennDiagram::calculate.overlap(ss)
go <- unlist(
  lapply(
    1:length(ss), function(j) {
           combn(names(ss), j, simplify = FALSE)
           }),
  recursive = FALSE
  )
names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
str(combs)


# get jaccard index matrix of intercellular leading edge genes 
# from different enrichments 
jmat = EnrichmentJaccard(gsealist = g0.sub, 
                         indexedgenes = li.index, 
                         returnJaccardMtx = TRUE)
jmat = jmat$jaccard_matrix_list

# load of baseline module expression across donors only of 
# leading edge genes from baseline enrichments 
ds = readRDS(here('mid_res/baseline_response/dataV3/ds.rds'))
data.table::fwrite(ds,file = paste0(here('git_ignore/ds.csv')),sep = ',')

# created shorter names 
new.names = data.table::fread(
  file = here('mid_res/baseline_response/dataV3/baseline.module.name.shortened.txt'),
  sep = '\t'
  )
new.names$cname2 = paste(new.names$celltype, new.names$shortname, sep = ' :: ')
colnames(ds) = new.names$cname2

# split module expression by cell types to subtract jacard similariry index 
# from intracellular spearman correlation coefficient caluclated below 
# also add shorter names names 
ds2 = ds %>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column('cname2') %>% 
  full_join(new.names, by = "cname2") %>% 
  select(-c('module', 'cname', 'shortname')) %>% 
  select(celltype, everything())
# split
ds2 = split(ds2,f = ds2$celltype)

# calculate Spearman correlation for intracellular correlations 
ds.cor = lapply(ds2, function(x) { 
  mtr =  x %>%  
    select(-c('celltype')) %>% 
    remove_rownames() %>% 
    column_to_rownames('cname2') %>% 
    t() 
  return(Hmisc::rcorr(mtr, type = 'spearman')$r)
  })

# subset to celltypes with multiple enrichments 
ds.cor = ds.cor[names(jmat)]

# append the names of the jaccard matrix with cell type so that format is 
# celltype :: module name to match correlation matrix
for (i in 1:length(jmat)) {
  mt = jmat[[i]]
  colnames(mt) = rownames(mt) = plyr::mapvalues(
    rownames(mt), 
    from = new.names$module, 
    to = new.names$shortname
    )
  colnames(mt) = rownames(mt) = paste(names(jmat)[i], colnames(mt),sep = ' :: ')
  jmat[[i]] = mt
}

# calculate shared latent informaiton of intracellular correlations
# subtract jaccard similariry from the spearman correlation coefficient 
# diagonal corrects to 0 
sli = list()
for (i in 1:length(jmat)) {
  stopifnot(isTRUE(all.equal(rownames(ds.cor[[i]]), rownames(jmat[[i]]))) )
  sli[[i]] = ds.cor[[i]] - jmat[[i]]
}

# now calculate the full spearman correlation matrix without subtracting jmat 
# replace the matrix values from intracellular correlations with the SLI values 
# which are stored in sli list by celltype 
spearmanmat = Hmisc::rcorr(ds, type = 'spearman')
saveRDS(spearmanmat, file = paste0(datapath, 'spearmanmat.rds'))
spearmanmat = readRDS(file = here('mid_res/baseline_response/dataV3/spearmanmat.rds'))
rhomat = spearmanmat$r
mat = rhomat
for (i in 1:length(sli)) {
  # ged index of cols and rows 
  row.replace = which(rownames(mat) %in% rownames(sli[[i]]))
  col.replace = which(colnames(mat) %in% colnames(sli[[i]]))
  
  # we need to be replacing values along the square diagonal of the matrix 
  # confirm these are the same 
  stopifnot(isTRUE(all.equal(row.replace, col.replace)))
  
  # check structure 
  stopifnot(isTRUE(all.equal(
    # full spearman matrix subset by rows of celltype i 
    mat[row.replace, col.replace],
    # original spearman correlation matrix for celltype i 
    ds.cor[[i]]
    )))
  
  # replace iteratively 
  mat[row.replace,col.replace] = sli[[i]]
  
  # this should be going down with each iteration 
  print(sum(mat))

}

# save SLI corrected matrix 
saveRDS(mat,file = paste0(datapath,'mat.rds'))
mat = readRDS(file = here('mid_res/baseline_response/dataV3/mat.rds'))

################################
# Part II visualization
################################

# plot pre and post adjustment correlation map 
# spearmancorrelation coefficient rhomat
cu = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(11)
range <- max(abs(rhomat))
# without clustering
dev.off()
pdf(file = paste0(figpath,'preadjusted.cormat.baseline.pdf'),width = 6, height = 5)
pheatmap::pheatmap(rhomat, color = cu, 
                   cluster_rows = FALSE, cluster_cols = FALSE, 
                   border_color = NA,
                   breaks = seq(-range, range, length.out = 11),
                   fontsize_row = 5, fontsize_col = 0.01)
dev.off()

# plot sli corrected matrix 
pdf(file = paste0(figpath,'post.SLIadjusted.cormat.baseline.pdf'),width = 6, height = 5)
pheatmap::pheatmap(mat, color = cu, 
                   cluster_rows = FALSE, cluster_cols = FALSE, 
                   border_color = NA,
                   breaks = seq(-range, range, length.out = 11),
                   fontsize_row = 5, fontsize_col = 0.01)
dev.off()

# figure  -- full network visualization as matrix (non pruned)
# cluster the square matrix 
#diag(mat)[diag(mat) > 0] = 0 
range <- max(abs(mat))
cu3 = BuenColors::jdb_palette('solar_flare', type = 'continuous') %>%  as.vector
cu3 = cu3[seq(from = 0 , to = 1000,length.out = 20)]
rownames(mat) = str_replace_all(string = rownames(mat),pattern = '_',replacement = ' ')
pdf(file = paste0(figpath,'post.clustered.SLIadjusted.cormat.baseline.pdf'),width = 8, height = 6.5)
pheatmap::pheatmap(mat, 
                   color = cu3, 
                   #cluster_rows = FALSE, cluster_cols = FALSE, 
                   border_color = NA,
                   treeheight_row = 10, 
                   treeheight_col = 20,
                   breaks = seq(-range, range, length.out = 20),
                   fontsize_row = 5.5,
                   fontsize_col = 0.01)
dev.off()


# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scglmmr_0.1.0   magrittr_2.0.1  here_1.0.1      forcats_0.5.1   stringr_1.4.0   dplyr_1.0.4     purrr_0.3.4    
# [8] readr_1.4.0     tidyr_1.1.2     tibble_3.0.6    ggplot2_3.3.3   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] lme4_1.1-26                 tidyselect_1.1.0            RSQLite_2.2.7               AnnotationDbi_1.52.0       
# [5] htmlwidgets_1.5.3           grid_4.0.5                  BiocParallel_1.24.1         scatterpie_0.1.7           
# [9] munsell_0.5.0               codetools_0.2-18            statmod_1.4.35              withr_2.4.3                
# [13] colorspace_2.0-0            GOSemSim_2.16.1             Biobase_2.50.0              knitr_1.39                 
# [17] rstudioapi_0.13             stats4_4.0.5                ggsignif_0.6.0              DOSE_3.16.0                
# [21] Rdpack_2.1.1                MatrixGenerics_1.2.1        emmeans_1.5.4               GenomeInfoDbData_1.2.4     
# [25] polyclip_1.10-0             bit64_4.0.5                 farver_2.0.3                pheatmap_1.0.12            
# [29] rprojroot_2.0.2             downloader_0.4              coda_0.19-4                 vctrs_0.4.1                
# [33] generics_0.1.2              TH.data_1.0-10              xfun_0.30                   doParallel_1.0.16          
# [37] R6_2.5.0                    GenomeInfoDb_1.26.7         graphlayouts_0.7.2          locfit_1.5-9.4             
# [41] pals_1.7                    bitops_1.0-6                cachem_1.0.4                fgsea_1.16.0               
# [45] DelayedArray_0.16.3         assertthat_0.2.1            scales_1.1.1                multcomp_1.4-16            
# [49] ggraph_2.0.5                nnet_7.3-15                 enrichplot_1.10.2           gtable_0.3.0               
# [53] egg_0.4.5                   tidygraph_1.2.0             sandwich_3.0-0              rlang_1.0.2                
# [57] slanter_0.2-0               splines_4.0.5               rstatix_0.7.0               dichromat_2.0-0            
# [61] broom_0.7.5                 checkmate_2.0.0             BiocManager_1.30.10         reshape2_1.4.4             
# [65] abind_1.4-5                 modelr_0.1.8                backports_1.2.1             qvalue_2.22.0              
# [69] Hmisc_4.5-0                 clusterProfiler_3.18.1      tools_4.0.5                 ellipsis_0.3.2             
# [73] gplots_3.1.1                RColorBrewer_1.1-2          BiocGenerics_0.36.1         Rcpp_1.0.6                 
# [77] plyr_1.8.6                  progress_1.2.2              base64enc_0.1-3             zlibbioc_1.36.0            
# [81] RCurl_1.98-1.3              prettyunits_1.1.1           ggpubr_0.4.0                rpart_4.1-15               
# [85] viridis_0.5.1               cowplot_1.1.1               S4Vectors_0.28.1            zoo_1.8-8                  
# [89] SummarizedExperiment_1.20.0 haven_2.3.1                 ggrepel_0.9.1               cluster_2.1.2              
# [93] fs_1.5.0                    variancePartition_1.25.6    data.table_1.14.0           DO.db_2.9                  
# [97] openxlsx_4.2.3              reprex_1.0.0                mvtnorm_1.1-1               packrat_0.7.0              
# [101] matrixStats_0.58.0          hms_1.0.0                   GSVA_1.38.2                 xtable_1.8-4               
# [105] pbkrtest_0.5-0.1            RhpcBLASctl_0.21-247.1      XML_3.99-0.6                rio_0.5.16                 
# [109] jpeg_0.1-8.1                BuenColors_0.5.6            readxl_1.3.1                IRanges_2.24.1             
# [113] gridExtra_2.3               compiler_4.0.5              maps_3.4.0                  KernSmooth_2.23-18         
# [117] crayon_1.4.1                shadowtext_0.0.9            minqa_1.2.4                 htmltools_0.5.2            
# [121] ggfun_0.0.4                 Formula_1.2-4               lubridate_1.7.9.2           DBI_1.1.1                  
# [125] tweenr_1.0.2                dbplyr_2.1.0                MASS_7.3-53.1               boot_1.3-27                
# [129] Matrix_1.3-2                car_3.0-10                  cli_3.3.0                   rbibutils_2.0              
# [133] parallel_4.0.5              igraph_1.2.6                GenomicRanges_1.42.0        pkgconfig_2.0.3            
# [137] rvcheck_0.1.8               foreign_0.8-81              foreach_1.5.1               xml2_1.3.2                 
# [141] annotate_1.68.0             XVector_0.30.0              GeneOverlap_1.26.0          estimability_1.3           
# [145] rvest_0.3.6                 digest_0.6.27               graph_1.68.0                cellranger_1.1.0           
# [149] fastmatch_1.1-0             htmlTable_2.1.0             edgeR_3.32.1                GSEABase_1.52.1            
# [153] curl_4.3                    gtools_3.8.2                nloptr_1.2.2.2              nlme_3.1-152               
# [157] lifecycle_1.0.0             jsonlite_1.7.2              aod_1.3.1                   carData_3.0-4              
# [161] mapproj_1.2.8               viridisLite_0.3.0           limma_3.46.0                pillar_1.4.7               
# [165] lattice_0.20-41             fastmap_1.1.0               httr_1.4.2                  survival_3.2-10            
# [169] GO.db_3.12.1                glue_1.6.2                  zip_2.1.1                   iterators_1.0.13           
# [173] png_0.1-7                   bit_4.0.4                   ggforce_0.3.3               stringi_1.5.3              
# [177] blob_1.2.1                  org.Hs.eg.db_3.12.0         latticeExtra_0.6-29         caTools_1.18.1             
# [181] memoise_2.0.0 
```

Visualize results of baseline high responder networks, part 1.  
mid_res/baseline_response/5b_network_construction_and_visualization.r

``` r
# baseline nettwork correlation 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
source('functions/MattPMutils.r')

# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")

# load SLI corrected matrix 
mat = readRDS(file = here('mid_res/baseline_response/dataV3/mat.rds'))

# Create mDC and innate sub-network 
dn = data.frame(mods = rownames(mat)) %>%
  mutate(name = mods) %>% 
  separate(name, into = c('celltype', 'module'),sep = ' :: ')
innate.sub = dn %>% 
  filter(celltype %in% c('CD14_Mono', 'CD16_Mono', 'MAIT_Like', 'mDC', 'BC_Naive'))
ms = innate.sub$mods

# rm QCd modules 
#m.rm = readRDS(here('mid_res/baseline_response/dataV3/m.rm.rds'))
m.rm = c(
  "CD14_Mono :: M111.1 viral sensing IRF2",
  "MAIT_Like :: Kegg Ag Presentation",           
  "MAIT_Like :: reactome interferon alpha beta"
)
ms = ms[!ms %in% m.rm]

# calculate adjusted p values across the spearman correlation matrix 
# load uncorrected matrix Hmisc obejct containing correlation p values 
spearmanmat = readRDS(file = here('mid_res/baseline_response/dataV3/spearmanmat.rds'))
padj = p.adjust.cormat(hmisc.cor =  spearmanmat, method = 'fdr')
saveRDS(padj, file = paste0(datapath,'padj.rds'))

# filter to the innate subnetwork 
# filter edges in SLI adjusted network
# only include correlations with adjusted p < 0.05
mat2 = mat
mat2 = mat2[ms, ms]
padj = padj[ms, ms]
stopifnot(isTRUE(all.equal(colnames(padj), colnames(mat2))))
# filter based on adjusted p 
mat2[padj > 0.05] <- 0

# save the pruned mat2 
saveRDS(mat2,file = paste0(datapath,'mat2.rds'))
mat2 = readRDS(file = here('mid_res/baseline_response/dataV3/mat2.rds'))


# make graph of the strongly linked edges pruned above
net <- graph_from_adjacency_matrix(
  mat2, weighted = TRUE,
  mode = 'undirected',
  diag = FALSE
  )


# prune the graph further to retain links above the median weight
med.weight <- median(E(net)$weight)
mat3 = mat2
mat3[mat3 < med.weight] <- 0
saveRDS(mat3,file = paste0(datapath,'mat3.rds'))
mat3 = readRDS(file = here('mid_res/baseline_response/dataV3/mat3.rds'))

# make a subhraph with stonger connections above prev median weight.
net <- graph_from_adjacency_matrix(
  mat3,
  weighted = TRUE,
  mode = 'undirected',
  diag = FALSE
  )

# create network annotations frame for vertices 
d = data.frame(signal = colnames(mat3)) %>%
  mutate(s = signal) %>%
  separate(s,into = c('celltype', 'module'),sep = ' :: ')

# specify vertex attributes 
V(net)$celltype = d$celltype
V(net)$module = d$module

# calculate network degree and hubs / authority scores (same for undirected)
V(net)$degree <- degree(net)                        
V(net)$hubs <- hub.score(net)$vector                
V(net)$authorities <- authority.score(net)$vector   

# add vertex property information to `d``
d$degree = V(net)$degree
d$hubscore = V(net)$hubs

# add d7 correlation to d 
d7.cor = readRDS(file = here('mid_res/baseline_response/dataV3/d7.cor.rds'))
d7.cor.p = d7.cor$P[1, -1][ms]
d7.cor.rho = d7.cor$r[1, -1][ms]
# check orders correct 
stopifnot(isTRUE(all.equal(names(d7.cor.p), d$signal)))
d$d7cor.p = d7.cor.p
d$d7cor.rho = d7.cor.rho

# add this information to the network 
V(net)$d7cor.p = d7.cor.p
V(net)$d7cor.rho = d7.cor.rho

# specify edge width as the weight 
E(net)$width = E(net)$weight 

# save network 
saveRDS(net,file = paste0(datapath,'net.rds'))
net = readRDS(file = here('mid_res/baseline_response/dataV3/net.rds'))

############################
# plot hubs 
signal.highlight = d %>% filter(celltype == 'CD14_Mono') %$% signal
signal.highlight2 = d %>% filter(celltype == 'CD16_Mono') %$% signal

cu3 = c('#FFD38F', '#F4A69B', '#A7DDEA', '#8ACFC3', '#9FABC4')

p = 
  ggplot(d, aes(y = reorder(module, hubscore), x = hubscore ,  fill = celltype, label = signal )) + 
  theme_bw() + 
  geom_point(shape =21, size = 3.5) + 
  #ggsci::scale_fill_npg(alpha = 0.8) + 
  scale_fill_manual(values = cu3) + 
  ylab('') + 
  theme(legend.position = c(0.8,0.15), legend.key.size = unit(0.2,units = 'cm')) +
  theme(axis.text = element_text(color = 'black')) +
  xlab('Hub Score') + 
  theme(aspect.ratio = 1.1)  +
  ggrepel::geom_text_repel(data = d %>% filter(signal %in% signal.highlight & hubscore > 0.75), 
                           size = 2.5, nudge_y = 0, nudge_x = -0.3, seed = 2, segment.size = 0.1,
                           force = 40,
                           max.overlaps = 10) + 
  ggrepel::geom_text_repel(data = d %>% filter(signal %in% signal.highlight2 & hubscore > 0.85), 
                           size = 2.5, nudge_y = 0, nudge_x = -0.3, box.padding = 0.4, seed = 1, 
                           max.overlaps =10,force = 40,
                           segment.size = 0.1) 
p
ggsave(p, filename = paste0(figpath, 'innate.subnetwork.hub.pdf'), width = 6.5, height = 6.5)

# specify colors for nodes 
cu = c( col.alpha('orange', alpha = 0.5), ggsci::pal_npg(alpha = 0.5)(4))
c.celltype = cu[factor(V(net)$celltype)]
# layout network.
lay <- layout_in_circle(net)

# specify celltypes for leend 
cts = str_replace_all(string = levels(factor(V(net)$celltype)),pattern = '_',replacement = ' ')


# version with vertiices highlighted 
# specify sve path for subgraph plots 
figpath3 = here('mid_res/baseline_response/figuresV3/subgraphsLABELED/'); dir.create(figpath3, recursive = TRUE)
# plot the subgraphs 
for (i in 1:length(unique(V(net)))) {
  
  # highlight edges 
  # specify subset highlighted 
  edge.highlight.t = incident(net, v = V(net)[i], mode="all")
  # for savig 
  signal = str_replace_all( names(V(net)[i]), pattern = ' :: ', replacement = '  ')
  
  # make a new network 
  net.sp = net 
  # set size for highlighted edge 
  E(net.sp)$width = 1.1
  # remove the other edges for visualization 
  ot = E(net.sp)[!E(net.sp) %in% edge.highlight.t]
  net.sp <- delete_edges(net.sp, edges = ot)
  E(net.sp)$color = col.alpha('black',0.7)
  # plot network 
  pdf(file = paste0(figpath3,signal,'.subnetwork.pdf'),width = 10, height = 10)
  plot(net.sp, 
       layout = lay, 
       vertex.label = names(V(net.sp)),
       vertex.label.cex=0.3, 
       edge.size = E(net.sp)$weight,
       vertex.size = log(V(net.sp)$degree+ 1)*4,
       edge.curved = 0.1,
       vertex.color = c.celltype,
       vertex.size=degree(net.sp)) 
  legend(x=0.75, y=1.2, 
         cts,
         pch=21, 
         pt.bg=cu,
         pt.cex=1,  cex=.5, bty="n",ncol=1)
  dev.off()
  
}

#### This not run in published workflow 
### Commented out for published workflow bc redundant w code below -- this used to make figs. 
# specify sve path for subgraph plots 
# figpath2 = here('mid_res/baseline_response/figuresV3/subgraphs/'); dir.create(figpath2, recursive = TRUE)
# # plot the subgraphs 
# for (i in 1:length(unique(V(net)))) {
#   
#   # highlight edges 
#   # specify subset highlighted 
#   edge.highlight.t = incident(net, v = V(net)[i], mode="all")
#   # for savig 
#   signal = str_replace_all( names(V(net)[i]), pattern = ' :: ', replacement = '  ')
# 
#   # make a new network 
#   net.sp = net 
#   # set size for highlighted edge 
#   E(net.sp)$width = 1.1
#   # remove the other edges for visualization 
#   ot = E(net.sp)[!E(net.sp) %in% edge.highlight.t]
#   net.sp <- delete_edges(net.sp, edges = ot)
#   E(net.sp)$color = col.alpha('black',0.7)
#   # plot network 
#   pdf(file = paste0(figpath2,signal,'.subnetwork.pdf'),width = 4.5, height = 4.5)
#   plot(net.sp, 
#        layout = lay, 
#        vertex.label = NA,
#        vertex.label.cex=0.3, 
#        edge.size = E(net.sp)$weight,
#        vertex.size = log(V(net.sp)$degree+ 1)*4,
#        edge.curved = 0.1,
#        vertex.color = c.celltype,
#        vertex.size=degree(net.sp)) 
#   legend(x=0.75, y=1.2, 
#          cts,
#          pch=21, 
#          pt.bg=cu,
#          pt.cex=1,  cex=.5, bty="n",ncol=1)
#   dev.off()
# }
# 
```

Visualize results of baseline high responder networks, part 2.  
mid_res/baseline_response/5c_network_correlations.r

``` r
# baseline nettwork correlation 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(igraph))
source('functions/MattPMutils.r')

# set paths 
datapath = here("mid_res/baseline_response/dataV3")
figpath = here("mid_res/baseline_response/figuresV3/network_correlations/");
dir.create(figpath,recursive = TRUE)


# load of baseline module expression across donors only of 
# leading edge genes from baseline enrichments 
ds = readRDS(here('mid_res/baseline_response/dataV3/ds.rds'))
#data.table::fwrite(ds,file = paste0(here('git_ignore/ds.csv')),sep = ',')
# created shorter names -- read these in 
new.names = data.table::fread(
  file = here('mid_res/baseline_response/dataV3/baseline.module.name.shortened.txt'),
  sep = '\t'
)
new.names$cname2 = paste(new.names$celltype, new.names$shortname, sep = ' :: ')
colnames(ds) = new.names$cname2

# fix subject names 
dp = ds %>%
  as.data.frame() %>% 
  rownames_to_column('subject') %>% 
  mutate(subject = str_sub(subject, 1,3)) %>% 
  column_to_rownames('subject')

# read matrix and network 
net = readRDS(file = here('mid_res/baseline_response/dataV3/net.rds'))
edf = as_long_data_frame(net)

# readadj p vas for comparison 
padj = readRDS(file = here('mid_res/baseline_response/dataV3/padj.rds'))

# aes set 
plotattr = list(
  theme_bw(),
  geom_point(shape = 21, size = 3.5, stroke = 0.8), 
  theme(axis.title.y = element_text(size = 8)), 
  theme(axis.title.x = element_text(size = 8)), 
  scale_fill_manual(values = c(col.alpha("red",0.7), col.alpha("dodgerblue", 0.7))),
  theme(legend.position = 'none', legend.key.size = unit(0.29, units = 'cm')), 
  theme(aspect.ratio = 1) 
)

# specify response in vis.  
high.responders = c("205","207","209","212", "215",
                    "234","237","245","250","256")

for (i in 1:nrow(edf)) {
 
  # get edge to plot from the data framed network  
  mod.names = c(edf[i, ]$from_name, edf[i, ]$to_name)
  
  cplot = dp %>% 
    select(all_of(mod.names)) %>% 
    rownames_to_column('subject') %>% 
    mutate(response = ifelse(subject %in% high.responders, 'high', 'low'))
  
  
  ctp = cor.test(cplot[ , 2], cplot[ ,3], method = 'spearman', exact = FALSE)$p.value
  adjusted.p = padj[mod.names[1], mod.names[2]]
  print(ctp < adjusted.p)

  # calculate and vis. correlation across all subjects; color by response. 
  p = ggpubr::ggscatter(cplot, 
                        x = mod.names[1], 
                        y = mod.names[2],
                        conf.int = FALSE, 
                        color = col.alpha('white','0.01'),
                        add.params = list(color = "black", fill = "grey"))  +
    plotattr +
    aes(fill = response) 
  
  # save name by cell type first 
  modsave = str_replace_all(mod.names,pattern = ' :: ', replacement = '..')
  ggsave(p, filename = paste0(figpath, modsave[1], '___', modsave[2], 'cor.pdf'), width = 3.1, height = 3.1)
  
}

# draw a legend 
fp2 = here('mid_res/baseline_response/figuresV3/')
p2 = p + theme(legend.position = 'top')
legend <- cowplot::get_legend(p2)
pdf(file = paste0(fp2, 'LEGEND.pdf'),width = 2, height = 1)
grid::grid.draw(legend)
dev.off()
```

### Fig 4. Early kinetics of baseline states <a name="fig4.4"></a>

Fit single cell mixed model to test day 1 post vaccination kinetics of
baseline cell phenotypes defined above.  
mid_res/baseline_response/6_sc.kinetic.singlecellmodel.r

``` r
# Early kinetics of baseline states 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))
# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")


# load single cell data 
h1 = readRDS(file = here('data/h1h5_annotated_with_meta.rds'))
md = h1@meta.data %>% 
  filter(cohort == 'H1N1') %>% 
  filter(time_cohort == 'd1') %>% 
  mutate(group_id = factor(adjmfc.time, levels = c("d0 low", "d1 low", "d0 high", "d1 high"))) %>% 
  mutate(subjectid = factor(sampleid)) %>% 
  mutate(sex = factor(gender)) %>% 
  mutate(scaled.age = as.numeric(scale(age))) %>% 
  mutate(celltype = celltype_joint) 
  
# add a covariate for number of cells per sample  
ncell_met = md %>% group_by(sample) %>% summarize(ncells = n())
md$ncell = plyr::mapvalues(x = md$sample, from = ncell_met$sample, to = ncell_met$ncells)
md$ncell = as.numeric(md$ncell)
md$log10ncell = log10(md$ncell) 

# subset normalized RNA data 
norm.rna = h1@data[ ,rownames(md)]


# get leading edge genes from cur. baseline mods 
g0.sub = readRDS(file = here('mid_res/baseline_response/dataV3/g0.sub.rds'))
li.g0 = LeadingEdgeIndexed(gsea.result.list = g0.sub, padj.threshold = 0.05)
li.g0 = base::Filter(length, li.g0)

# metadata by cell type 
cts = names(li.g0)
md = md %>% filter( celltype %in% cts )
ct.md = split( md, f = md$celltype )

# module score for each cell type of specific baseline enriched leading edge genes.
mod_scores = list()
for (i in 1:length(ct.md)) {
  
  # init data for subset i 
  rna = norm.rna[ ,rownames(ct.md[[i]])]
  mod.list = li.g0[[i]]
  
  # calculate single cell score for baseline-enriched module 
  mod_scores[[i]] = WeightedCellModuleScore(
    gene_matrix = rna, 
    module_list = mod.list, 
    threshold = 0, 
    cellwise_scaling = TRUE, 
    return_weighted = FALSE 
    )
  
  # add a "null" score of Gaussian noise as a reference 
  mod_scores[[i]]$null = rnorm(n = nrow(mod_scores[[i]]), mean = 0, sd = 1)
}

# specify save paths for marginal means plots.
plot_savepath1 = paste0(figpath, "/marginalmeans.m1/"); dir.create(plot_savepath1)
plot_savepath2 = paste0(figpath, "/marginalmeans.m2/"); dir.create(plot_savepath2)


# specify the 2 models 
f1 = 'modulescore ~ 0 + group_id + (1|subjectid)'
f2 = 'modulescore ~ 0 + group_id + log10ncell + scaled.age + sex + (1|subjectid)'


# fit sc mod mixed model on module scores. 
mm1 = mm2 = list()
for (i in 1:length(ct.md)) {
  
  stopifnot( nrow(ct.md[[i]]) == nrow(mod_scores[[i]]) )
  
  # formula 1
  mm1[[i]] = FitLmerContrast(module_data_frame = mod_scores[[i]], 
                              celltype_column = 'celltype', 
                              metadata = ct.md[[i]], 
                              lmer_formula = f1, 
                              plotdatqc = FALSE, 
                              fixed_effects = NULL,
                              figpath = plot_savepath1)
  
  # formula 2 
  mm2[[i]] = FitLmerContrast(module_data_frame = mod_scores[[i]], 
                             celltype_column = 'celltype', 
                             metadata = ct.md[[i]], 
                             lmer_formula = f2, 
                             plotdatqc = FALSE, 
                             fixed_effects = NULL,
                             figpath = plot_savepath2)  

}

mm1 = do.call(rbind, mm1)
mm2 = do.call(rbind, mm2)

saveRDS(mm1,file = paste0(datapath, 'mm1.rds'))
saveRDS(mm2,file = paste0(datapath, 'mm2.rds'))


sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] SeuratObject_4.0.0 Seurat_4.0.1       scglmmr_0.1.0      magrittr_2.0.1     here_1.0.1        
# [6] forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0       
# [11] tidyr_1.1.2        tibble_3.0.6       ggplot2_3.3.3      tidyverse_1.3.0   
# 
# loaded via a namespace (and not attached):
# [1] estimability_1.3            scattermore_0.7             coda_0.19-4                
# [4] bit64_4.0.5                 irlba_2.3.3                 multcomp_1.4-16            
# [7] DelayedArray_0.16.3         rpart_4.1-15                data.table_1.14.0          
# [10] RCurl_1.98-1.3              generics_0.1.0              BiocGenerics_0.36.1        
# [13] cowplot_1.1.1               TH.data_1.0-10              RSQLite_2.2.7              
# [16] shadowtext_0.0.9            RANN_2.6.1                  future_1.21.0              
# [19] bit_4.0.4                   enrichplot_1.10.2           spatstat.data_2.1-0        
# [22] xml2_1.3.2                  lubridate_1.7.9.2           httpuv_1.5.5               
# [25] SummarizedExperiment_1.20.0 assertthat_0.2.1            viridis_0.5.1              
# [28] hms_1.0.0                   promises_1.2.0.1            caTools_1.18.1             
# [31] dbplyr_2.1.0                readxl_1.3.1                igraph_1.2.6               
# [34] DBI_1.1.1                   htmlwidgets_1.5.3           spatstat.geom_2.0-1        
# [37] stats4_4.0.5                ellipsis_0.3.1              ggpubr_0.4.0               
# [40] backports_1.2.1             annotate_1.68.0             deldir_0.2-10              
# [43] MatrixGenerics_1.2.1        vctrs_0.3.6                 Biobase_2.50.0             
# [46] ROCR_1.0-11                 abind_1.4-5                 cachem_1.0.4               
# [49] withr_2.4.3                 ggforce_0.3.3               packrat_0.7.0              
# [52] emmeans_1.5.4               sctransform_0.3.2           goftest_1.2-2              
# [55] cluster_2.1.2               DOSE_3.16.0                 lazyeval_0.2.2             
# [58] crayon_1.4.1                edgeR_3.32.1                pkgconfig_2.0.3            
# [61] labeling_0.4.2              tweenr_1.0.2                GenomeInfoDb_1.26.7        
# [64] nlme_3.1-152                rlang_0.4.10                globals_0.14.0             
# [67] lifecycle_1.0.0             miniUI_0.1.1.1              sandwich_3.0-0             
# [70] downloader_0.4              modelr_0.1.8                cellranger_1.1.0           
# [73] rprojroot_2.0.2             polyclip_1.10-0             GSVA_1.38.2                
# [76] matrixStats_0.58.0          lmtest_0.9-38               graph_1.68.0               
# [79] Matrix_1.3-2                carData_3.0-4               boot_1.3-27                
# [82] zoo_1.8-8                   reprex_1.0.0                ggridges_0.5.3             
# [85] pheatmap_1.0.12             png_0.1-7                   viridisLite_0.3.0          
# [88] bitops_1.0-6                KernSmooth_2.23-18          blob_1.2.1                 
# [91] qvalue_2.22.0               parallelly_1.23.0           rstatix_0.7.0              
# [94] S4Vectors_0.28.1            ggsignif_0.6.0              scales_1.1.1               
# [97] memoise_2.0.0               GSEABase_1.52.1             plyr_1.8.6                 
# [100] ica_1.0-2                   gplots_3.1.1                zlibbioc_1.36.0            
# [103] compiler_4.0.5              scatterpie_0.1.7            RColorBrewer_1.1-2         
# [106] lme4_1.1-26                 fitdistrplus_1.1-3          cli_2.5.0                  
# [109] XVector_0.30.0              listenv_0.8.0               patchwork_1.1.1            
# [112] pbapply_1.4-3               mgcv_1.8-34                 MASS_7.3-53.1              
# [115] tidyselect_1.1.0            stringi_1.5.3               GOSemSim_2.16.1            
# [118] locfit_1.5-9.4              ggrepel_0.9.1               GeneOverlap_1.26.0         
# [121] grid_4.0.5                  fastmatch_1.1-0             tools_4.0.5                
# [124] future.apply_1.7.0          parallel_4.0.5              rio_0.5.16                 
# [127] rstudioapi_0.13             foreign_0.8-81              gridExtra_2.3              
# [130] farver_2.0.3                Rtsne_0.15                  ggraph_2.0.5               
# [133] digest_0.6.27               rvcheck_0.1.8               BiocManager_1.30.10        
# [136] shiny_1.6.0                 Rcpp_1.0.6                  GenomicRanges_1.42.0       
# [139] car_3.0-10                  broom_0.7.5                 egg_0.4.5                  
# [142] later_1.1.0.1               RcppAnnoy_0.0.18            org.Hs.eg.db_3.12.0        
# [145] httr_1.4.2                  AnnotationDbi_1.52.0        colorspace_2.0-0           
# [148] tensor_1.5                  rvest_0.3.6                 XML_3.99-0.6               
# [151] fs_1.5.0                    reticulate_1.18             IRanges_2.24.1             
# [154] splines_4.0.5               uwot_0.1.10                 statmod_1.4.35             
# [157] spatstat.utils_2.1-0        graphlayouts_0.7.2          plotly_4.9.3               
# [160] xtable_1.8-4                jsonlite_1.7.2              nloptr_1.2.2.2             
# [163] tidygraph_1.2.0             ggfun_0.0.4                 R6_2.5.0                   
# [166] pillar_1.4.7                htmltools_0.5.2             mime_0.10                  
# [169] glue_1.4.2                  fastmap_1.1.0               minqa_1.2.4                
# [172] clusterProfiler_3.18.1      BiocParallel_1.24.1         codetools_0.2-18           
# [175] fgsea_1.16.0                mvtnorm_1.1-1               spatstat.sparse_2.0-0      
# [178] lattice_0.20-41             slanter_0.2-0               curl_4.3                   
# [181] leiden_0.3.7                gtools_3.8.2                zip_2.1.1                  
# [184] GO.db_3.12.1                openxlsx_4.2.3              survival_3.2-10            
# [187] limma_3.46.0                munsell_0.5.0               DO.db_2.9                  
# [190] GenomeInfoDbData_1.2.4      haven_2.3.1                 reshape2_1.4.4             
# [193] gtable_0.3.0                spatstat.core_2.0-0        
```

Visualize results of kinetic analysis above.  
mid_res/baseline_response/7_sc.kinetic.singlecellmodel.figures.r

``` r
# Early kinetics of baseline states 
# visualize model results
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))

# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")
source('functions/MattPMutils.r')

mm2 = readRDS(file = here('mid_res/baseline_response/dataV3/mm2.rds')) %>% 
  filter(!module == 'null') %>% 
  filter(! singular_fit == 1)
mm2$cm = paste(mm2$celltype, mm2$module,sep = ' :: ')

# change to shorter names 
new.names = data.table::fread(file = here('mid_res/baseline_response/dataV3/baseline.module.name.shortened.txt'), sep = '\t')
new.names$cname2 = paste(new.names$celltype, new.names$shortname, sep = ' :: ')
new.names$cname1 = paste(new.names$celltype, new.names$module, sep = ' :: ')
mm2$cm = plyr::mapvalues(x = mm2$cm, from = new.names$cname1, to = new.names$cname2)

# assign to 'd'
d = mm2

# plot innate subset 
ds = d %>% filter(celltype  %in% c( 'CD14_Mono', 'CD16_Mono', 'mDC', "MAIT_Like" )) 

# filter the modules that did not have a effect in the single cell model 
m.rm = ds %>%
  filter(estimatetime0_group2vs1 > 0.1) %>%  
  mutate(m1s = estimatetime0_group2vs1 - std.errortime0_group2vs1) %>%
  filter(m1s > -0.01) %$% cm
saveRDS(m.rm, file = paste0(datapath, 'm.rm.rds'))
ds = ds %>% filter(cm %in% m.rm)

# remove cell label from module name to reduce clutter 
ds$cm = gsub(".*:","",ds$cm)
ds$celltype = str_replace_all(ds$celltype, pattern = '_',replacement =' ')


pl = list(
  # baseline 
  geom_point(shape = 23, size = 3, color = 'black', fill = col.alpha('red', 0.5)),
  geom_segment(aes(x = (estimatetime0_group2vs1 + -1*std.errortime0_group2vs1),
                   xend = estimatetime0_group2vs1 + 1*std.errortime0_group2vs1,
                   yend = cm),
               color = col.alpha('red', 0.5), 
               size = 2), 
  # day 1 
  geom_point(data = ds, aes(x = estimatetime1vs0, y = cm), size = 3, shape = 23, 
               color = 'black', fill = col.alpha('#e2a359', 0.5)),
  geom_segment(aes(x = (estimatetime1vs0 + -1*std.errortime1vs0),
                     xend = estimatetime1vs0 + 1*std.errortime1vs0,
                     yend = cm), 
                 color = col.alpha('#e2a359', 0.5),
                 size = 2), 
  theme(
    axis.text.y = element_text(color = 'black', size = 10),
    strip.background = element_blank(),
    axis.title.x = element_text(size = 14), 
    axis.text.x = element_text(size = 14), 
    strip.text = element_text(size = 12), 
    panel.spacing.x=unit(2, "lines")
    )
  )

# mdc, CD14 mono, CD16 mono, mait-like
p = 
  ggplot(ds, aes(x = estimatetime0_group2vs1, y = reorder(cm, estimatetime1vs0 ))) + 
  facet_grid(vars(celltype), scales = 'free', space = 'free') +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  pl + 
  scale_x_continuous( breaks= scales::pretty_breaks(), expand = c(0.15,0)) + 
  ylab('') + 
  xlab('contrast effect size ') 
p
ggsave(p, filename = paste0(figpath, 'mm2.innate.mait.pdf'), width = 5, height = 6.5)

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] magrittr_2.0.1  here_1.0.1      forcats_0.5.1   stringr_1.4.0   dplyr_1.0.4     purrr_0.3.4     readr_1.4.0    
# [8] tidyr_1.1.2     tibble_3.0.6    ggplot2_3.3.3   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.6        plyr_1.8.6        pillar_1.4.7      compiler_4.0.5    cellranger_1.1.0  dbplyr_2.1.0     
# [7] tools_4.0.5       digest_0.6.27     packrat_0.7.0     jsonlite_1.7.2    lubridate_1.7.9.2 lifecycle_1.0.0  
# [13] gtable_0.3.0      pkgconfig_2.0.3   rlang_0.4.10      reprex_1.0.0      cli_2.5.0         rstudioapi_0.13  
# [19] DBI_1.1.1         haven_2.3.1       withr_2.4.3       xml2_1.3.2        httr_1.4.2        fs_1.5.0         
# [25] generics_0.1.0    vctrs_0.3.6       hms_1.0.0         rprojroot_2.0.2   grid_4.0.5        tidyselect_1.1.0 
# [31] data.table_1.14.0 glue_1.4.2        R6_2.5.0          readxl_1.3.1      farver_2.0.3      modelr_0.1.8     
# [37] backports_1.2.1   scales_1.1.1      ellipsis_0.3.1    rvest_0.3.6       assertthat_0.2.1  colorspace_2.0-0 
# [43] labeling_0.4.2    stringi_1.5.3     munsell_0.5.0     broom_0.7.5       crayon_1.4.1     
```

### Fig 4. Analysis of mRNA vaccine data to define induction of high responder phenotypes <a name="fig4.5"></a>

Process data from GSE171964.  
mid_res/mrna/mrna_1\_setup.r

``` r
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

# load metadata 
p = read.delim(file = here('data/GSE171964/GSE171964_geo_pheno.csv'), sep = ',')

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
```

Renormalize raw ADT data with `dsb::ModelNegativeADTnorm()` and manually
gate cells.  
mid_res/mrna/mrna_2\_gatemono.r

``` r
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
d$pp3 = ifelse(d$tx==TRUE & d$ty==TRUE & d$CD3.ADT < tlm, yes = '1',no = '0')

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
```

Test baseline high responder cell phenotypes in the same subsets before
and after mRNA vaccination with mixed model.  
mid_res/mrna/mrna_3\_baseline.sig.test.mono.r

``` r
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
suppressMessages(library(emmeans))
source('functions/MattPMutils.r')

# set paths 
datapath = file.path(here('mid_res/mrna/generated_data/'))
figpath = file.path(here('mid_res/mrna/figures/'))

# load baseline monocyte leadingedge index unique genes 
gs0 = readRDS(file = here('mid_res/baseline_response/dataV3/g0.sub.rds'))
li0 = LeadingEdgeIndexed(gsea.result.list = gs0, padj.threshold = 0.05)
li0 = li0$CD14_Mono

# define ifn sigs 
li0.ifn = grepl('IFN|interferon',x = names(li0))

# subset signatures to ifn annotated and nonn ifn annotated 
li0_ifn = li0[li0.ifn]
li0_other = li0[!li0.ifn]
ifn.sig = li0_ifn %>%  unlist() %>%  unique()
non.ifn.sig = li0_other %>% unlist() %>% unique()

# further prune non ifn to not include overlappint genes with ifn sigs.
both = intersect(ifn.sig, non.ifn.sig)
non.ifn.sig = non.ifn.sig[!non.ifn.sig %in% both]
sig.test = list('ifn' = ifn.sig, 'non.ifn' = non.ifn.sig, 'sig' = unique(unlist(li0)))

# save combined signature genes (e2k)
sig.genes = sig.test$sig
data.table::fwrite(list(sig.genes),file = paste0(datapath,'sig.txt'), sep = '\t')


# load monocyte gated CITE-seq data from pfizer data 
s.mono = readRDS('mid_res/mrna/generated_data/s.mono.rds')
s.mono = NormalizeData(s.mono,assay = 'RNA',normalization.method = 'LogNormalize')
# define umi matrix and metadata 
umi = s.mono@assays$RNA@data
md = s.mono@meta.data
# format metadata for lme4 
md$time = factor(md$day,levels = c('0', '1', '21', '22'))
md$pt_id = factor(as.character(md$pt_id))

# module score simple average for the 3 signatures defined above. 
mscore = WeightedCellModuleScore(gene_matrix = umi, 
                                 module_list = sig.test, 
                                 threshold = 0, 
                                 cellwise_scaling = FALSE, 
                                 return_weighted = FALSE)

# combine signature scores with meta.data 
dat.fit = cbind(mscore, md)
saveRDS(dat.fit, file = paste0(datapath, 'dat.fit.mono.rds'))
dat.fit = readRDS(file = here('mid_res/mrna/generated_data/dat.fit.mono.rds'))

# note strucure with one major outlier in cell number in this dataset: 
table(dat.fit$time, dat.fit$sample_id) %>% t()

# lmer formula for the 3 signatures 
f1 = 'sig ~ 0 + time + (1|pt_id)'
f2 = 'ifn ~ 0 + time + (1|pt_id)'
f3 = 'non.ifn ~ 0 + time + (1|pt_id)'


library(emmeans)
# fit model for each signature
m1 = lme4::lmer(formula = f1,data = dat.fit)
emm1 = emmeans::emmeans(m1,specs = ~time, lmer.df = 'asymptotic')
plot(emm1)
m2 = lme4::lmer(formula = f2,data = dat.fit)
emm2 = emmeans::emmeans(m2,specs = ~time, lmer.df = 'asymptotic')

m3 = lme4::lmer(formula = f3,data = dat.fit)
emm3 = emmeans::emmeans(m3,specs = ~time, lmer.df = 'asymptotic')

# contrast time differences 
clevels = levels(dat.fit$time)

#make custom contrasts 
c0 = c(1, 0, 0 ,0)
c1 = c(0, 1, 0 ,0)
c3 = c(0, 0, 1 ,0)
c4 = c(0, 0, 0 ,1)
contrast_list = list( "time1vs0" = c1 - c0, 'time22vs21' = c4 - c3 )
clist = list('sig' = emm1, 'ifn' = emm2, 'non.ifn' = emm3)

c.res = 
lapply(clist, function(x) { 
  emmeans::contrast(object = x, method = contrast_list) %>% 
    broom::tidy()
  } ) %>% 
  bind_rows(.id = 'signature')
data.table::fwrite(x = c.res, file = paste0(datapath,'c.res.mono.txt'), sep = '\t')  

# delta contrast 
cmat = emmeans::contrast(object = emm1, method = contrast_list)
pairs(cmat, reverse = TRUE)
# contrast              estimate      SE  df z.ratio p.value
# time22vs21 - time1vs0    0.118 0.00434 Inf 27.137  <.0001 
# 
# Degrees-of-freedom method: asymptotic 


# plotsingle cell distributionn and emmeans contrasts 
em_aes = list(theme_bw(), 
              coord_flip(), 
              theme(axis.title.y = element_text(size = 7),
                    axis.text.y = element_text(size = 7)), 
              scale_color_manual('grey')
              )

plot.aes = list(theme_bw(), ylab(label ='Baseline high responder\nCD14 Mono signature'))


cu = sapply(c('grey', '#e2a359', 'grey', '#e2a359'), col.alpha, 0.8) %>% unname()
# combined signature change emm in p1 and change y value in p0
p0 = ggplot(dat.fit, aes(x = time, y = sig, fill = time)) + 
  geom_violin(show.legend = F) + 
  plot.aes +
  xlab('time') + 
  scale_fill_manual(values = cu)+
  theme(axis.title.x = element_text(size = 12))
ggsave(p0, filename = paste0(figpath, 'sig.pdf'), width = 2.5, height = 3)
p1 = plot(emm1) + em_aes
ggsave(p1, filename = paste0(figpath, 'sig.emm.pdf'), width = 1, height = 3)

# ifn -- change emm in p1 and change y value in p0
p0 = ggplot(dat.fit, aes(x = time, y = ifn, fill = time)) + 
  geom_violin(show.legend = F) + 
  plot.aes +
  xlab('time') + 
  scale_fill_manual(values = cu)+
  theme(axis.title.x = element_text(size = 12))
ggsave(p0, filename = paste0(figpath, 'ifn.pdf'),  width = 2.5, height = 3)
p1 = plot(emm2) + em_aes
ggsave(p1, filename = paste0(figpath, 'ifn.emm.pdf'), width = 1, height = 3)


# non-ifn -- change emm in p1 and change y value in p0 
p0 = ggplot(dat.fit, aes(x = time, y = non.ifn, fill = time)) + 
  geom_violin(show.legend = F) + 
  plot.aes +
  xlab('time') + 
  scale_fill_manual(values = cu)+
  theme(axis.title.x = element_text(size = 12))
ggsave(p0, filename = paste0(figpath, 'non.ifn.pdf'),  width = 2.5, height = 3)
p1 = plot(emm3) + em_aes
ggsave(p1, filename = paste0(figpath, 'non.ifn.emm.pdf'), width = 1, height = 3)
```

Run same test in DCs  
mid_res/mrna/mrna_3\_signature.test.mdc.r

``` r
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
suppressMessages(library(emmeans))

# set save paths 
datapath = file.path(here('mid_res/mrna/generated_data/'))
figpath = file.path(here('mid_res/mrna/figures/'))

# load baseline monocyte leadingedge index unique genes 
gs0 = readRDS(file = here('mid_res/baseline_response/dataV3/g0.sub.rds'))
li0 = LeadingEdgeIndexed(gsea.result.list = gs0, padj.threshold = 0.05)
li0 = li0$mDC
# define sigs 
li.0.m11  = grepl('LI.M11.0',x = names(li0))
li.0_m11 = li0[li.0.m11]
li0.other = li0[!li.0.m11]
m11.sig = li.0_m11 %>%  unlist() %>%  unique()
non.m11.sig = li0.other %>% unlist() %>% unique()
# further prune non ifn to not include overlappint genes with ifn sigs.
both = intersect(m11.sig, non.m11.sig)
non.m11.sig = non.m11.sig[!non.m11.sig %in% both]
sig.test = list('msig' = m11.sig, 'non.msig.sig' = non.m11.sig, 'sig' = unique(unlist(li0)))


# load monocyte gated CITE-seq data from pfizer data 
s.mdc = readRDS('mid_res/mrna/generated_data/s.mdc.rds')
s.mdc = NormalizeData(s.mdc,assay = 'RNA', normalization.method = 'LogNormalize')
# define umi matrix and metadata 
umi = s.mdc@assays$RNA@data
md = s.mdc@meta.data
# format metadata for lme4 
md$time = factor(md$day,levels = c('0', '1', '21', '22'))
md$pt_id = factor(as.character(md$pt_id))

# this caused a singular fit 
# add log10 n cell as covariate per sample 
# samplen = md %>% select(sx = sample_id) %>% group_by(sx) %>%  tally() %>% mutate(logncell = log10(n))
# md$log10ncell = plyr::mapvalues(x = md$sample_id, from = samplen$sx,to = samplen$logncell)

# module score for the 3 signatures defined above. 
mscore = WeightedCellModuleScore(gene_matrix = umi, 
                                 module_list = sig.test, 
                                 threshold = 0, 
                                 cellwise_scaling = FALSE, 
                                 return_weighted = FALSE)

# combine signature scores with meta.data 
dat.fit = cbind(mscore, md)
saveRDS(dat.fit, file = paste0(datapath, 'dat.fit.mdc.rds'))
dat.fit = readRDS(file = here('mid_res/mrna/generated_data/dat.fit.mdc.rds'))

# note strucure not as bad on outlier sample as mono
table(dat.fit$time, dat.fit$sample_id) %>% t()

# lmer formula for the 3 signatures 
f1 = 'msig ~ 0 + time + (1|pt_id)'
f2 = 'non.msig.sig ~ 0 + time + (1|pt_id)'
f3 = 'sig ~ 0 + time + (1|pt_id)'

library(emmeans)
# fit model for each signature
m1 = lme4::lmer(formula = f1, data = dat.fit)
emm1 = emmeans::emmeans(m1, specs = ~time, lmer.df = 'asymptotic')

m2 = lme4::lmer(formula = f2, data = dat.fit)
emm2 = emmeans::emmeans(m2,specs = ~time, lmer.df = 'asymptotic')

m3 = lme4::lmer(formula = f3, data = dat.fit)
emm3 = emmeans::emmeans(m3, specs = ~time, lmer.df = 'asymptotic')

# contrast time differences 
clevels = levels(dat.fit$time)

#make custom contrasts 
c0 = c(1, 0, 0 ,0)
c1 = c(0, 1, 0 ,0)
c3 = c(0, 0, 1 ,0)
c4 = c(0, 0, 0 ,1)
contrast_list = list( "time1vs0" = c1 - c0,
                      'time22vs21' = c4 - c3)
clist = list('msig' = emm1, 'non.msig.sig' = emm2, 'sig' = emm3)

c.res = 
  lapply(clist, function(x) { 
    emmeans::contrast(object = x, method = contrast_list) %>% 
      broom::tidy()
  } ) %>% 
  bind_rows(.id = 'signature')
data.table::fwrite(x = c.res, file = paste0(datapath,'c.res.mDC.txt'), sep = '\t')  


# delta contrast 
cmat = emmeans::contrast(object = emm1, method = contrast_list)
pairs(cmat, reverse = TRUE)
# contrast              estimate     SE  df z.ratio p.value
# time22vs21 - time1vs0  0.00908 0.0172 Inf 0.528   0.5978 

# plotsingle cell distributionn and emmeans contrasts 
em_aes = list(theme_bw(), 
              coord_flip(), 
              theme(axis.title.y = element_text(size = 7),
                    axis.text.y = element_text(size = 7)), 
              scale_color_manual('grey')
)

plot.aes = list(theme_bw(), ylab(label ='Baseline high responder\nmDC signature'))
cu = sapply(c('grey', '#e2a359', 'grey', '#e2a359'), col.alpha, 0.8) %>% unname()
# combined signature change emm in p1 and change y value in p0
p0 = ggplot(dat.fit, aes(x = time, y = msig, fill = time)) + 
  geom_violin(show.legend = F) + 
  plot.aes + 
  xlab('time') + 
  scale_fill_manual(values = cu)+
  theme(axis.title.x = element_text(size = 12))
p0
ggsave(p0, filename = paste0(figpath, 'msig_mDC.pdf'), width = 2.5, height = 3)
p1 = plot(emm1) + em_aes
ggsave(p1, filename = paste0(figpath, 'msig_mDC.emm.pdf'), width = 1, height = 3)

# ifn -- change emm in p1 and change y value in p0
p0 = ggplot(dat.fit, aes(x = time, y = non.msig.sig, fill = time)) + 
  geom_violin(show.legend = F) + 
  plot.aes +
  xlab('time') + 
  scale_fill_manual(values = cu)+
  theme(axis.title.x = element_text(size = 12))
ggsave(p0, filename = paste0(figpath, 'non.msig.sig_mDC.pdf'), width = 2.5, height = 3)
p1 = plot(emm2) + em_aes
ggsave(p1, filename = paste0(figpath, 'non.msig.sig_mDC.emm.pdf'), width = 1, height = 3)


# non-ifn -- change emm in p1 and change y value in p0 
p0 = ggplot(dat.fit, aes(x = time, y = sig, fill = time)) + 
  geom_violin(show.legend = F) + 
  plot.aes +
  xlab('time') + 
  scale_fill_manual(values = cu)+
  theme(axis.title.x = element_text(size = 12))
ggsave(p0, filename = paste0(figpath, 'sig_mDC.pdf'), width = 2.5, height = 3)
p1 = plot(emm3) + em_aes
ggsave(p1, filename = paste0(figpath, 'sig_mDC.emm.pdf'), width = 1, height = 3)
```

### Fig.5. Define and test AS03 specific cell phenotypes in high responders at baseline <a name="fig5.1"></a>

Aggregate log cpm in high and low responders at baseline  
mid_res/nat_adj/1.calc.lcpm.baselineH1_as03_modelgenes.r

``` r
# R version 4.0.5 
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(scglmmr))
source(here('functions/MattPMutils.r'))

# set save paths 
figpath = here("mid_res/nat_adj/figures/V4/"); dir.create(figpath)
datapath = here("mid_res/nat_adj/generated_data//V4/"); dir.create(datapath)

# define high responders 
high.responders = c("205","207","209","212","215","234","237","245","250","256")

# read pb data, subset to day 0 non adj, subset out day 0 metadata. 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))
pb = lapply(pb, function(x) x = x[ ,1:40])
cnames = gsub("~.*","",colnames(pb[[1]]))
pb = lapply(pb, function(x){
  x %>% as.data.frame() %>% setNames(nm = cnames) %>% as.matrix() 
})
d0 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV3/d0.rds'))
d0d = lapply(pb, function(x){ x = x[ , rownames(d0)]})

# make a list of genes indxed by celltype for genes to fit from H5 model 
av_tidy = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/git_ignore/av_tidy.rds'))
genes.test = lapply(av_tidy , function(x) unique(x$gene))

# check names are the same of indexed genes and average data
stopifnot(isTRUE(all.equal(names(d0d), names(genes.test))))

# calcualte log CPM of the baseline pseudobulk data
av = list()
for (i in 1:length(d0d)) {
  dge =  edgeR::DGEList( d0d[[i]] ) 
  dge = dge[genes.test[[i]], ]
  av[[i]]  = edgeR::cpm(dge, log = TRUE)
}
names(av) = names(d0d)
saveRDS(av,file = paste0(datapath,'av.rds'))

# tidy aggregated data 
av0 = list()
for (i in 1:length(pb)) {
  ct = names(av)[i]
  gs = rownames(av[[i]])
  av0[[i]] = GetTidySummary(
    av.exprs.list = av,
    celltype.index = i,
    genes.use = gs) %>%
    mutate(response = if_else(str_sub(sample, 1, 3) %in%  high.responders, 'High', "Low")) %>%
    mutate(response = factor(response, levels = c('Low', 'High')))
}
names(av0) = names(pb)
saveRDS(av0, file = paste0(datapath,'av0.rds'))
```

Validate AS03 specificity of combined signatures in external AS03
cohort.  
mid_res/nat_adj/2.AS03.innatesigs.vand.validationcohort.r

``` r
# R version 4.0.5 
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(scglmmr))
suppressMessages(library(magrittr))
source(here('functions/MattPMutils.r'))

# set save paths 
figpath = here("mid_res/nat_adj/figures/V4/")
datapath = here("mid_res/nat_adj/generated_data//V4/")

# define adjuvant signatures 
# leading edge combined signature 
# upregulated genes only 
# mono and mDC specific and combined genes
li.up = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li.up.rds'))
li2.up = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li2.up.rds'))
li.full = unlist(li.up,use.names = FALSE, recursive = TRUE) %>%  unique 
mono.genes  = li.up$CD14_Mono %>%  unlist() %>% unique()
mdc.genes = li2.up$mDC %>%  unlist() %>% unique()
as03.sig = list('AS03_signature' = li.full, "AS03_Monocyte" = mono.genes, 'AS03_mDC' = mdc.genes)

# vand validation for the highlighted signatures
vand.fit = readRDS(file = here('mid_res/vand/generated_data/fit1.rds'))
vand.rank = ExtractResult(model.fit.list = vand.fit,what = 'lmer.z.ranks', coefficient.number = 1, coef.name = 'delta')
gvand = FgseaList(rank.list.celltype = vand.rank, pathways = as03.sig)

# leading edge from each cell tpe 
dc.va.le = gvand$DNC %>%  filter(pathway == 'AS03_mDC') %$% leadingEdge %>%  unlist()
saveRDS(dc.va.le,file = paste0(datapath, 'dc.va.le.rds'))

mono.va.le = gvand$MNC %>%  filter(pathway == 'AS03_Monocyte') %$% leadingEdge %>%  unlist() 
saveRDS(mono.va.le, file = paste0(datapath, 'mono.va.le.rds'))


#mono
enrline = list(geom_line(color = "deepskyblue3", size = 2 ))
p = fgsea::plotEnrichment(pathway = as03.sig$AS03_Monocyte, stats = vand.rank$MNC) + enrline
# mDC
p = fgsea::plotEnrichment(pathway = as03.sig$AS03_mDC, stats = vand.rank$DNC) + enrline
ggsave(p, filename = paste0(figpath, 'dc.vand.enr.pdf'), width = 4, height = 2.5)
```

Test AS03 specific mDC and monocyte phenotypes in high vs low responders
at baseline.  
mid_res/nat_adj/3.natural.adjuvant.signatures.r

``` r
# R version 4.0.5 
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(scglmmr))
source(here('functions/MattPMutils.r'))

# set save paths 
figpath = here("mid_res/nat_adj/figures/V4/")
datapath = here("mid_res/nat_adj/generated_data/V4/")

# set plot themes to distinguisn between groups / signatures being tested 
# AS03 theme 
cu = c("grey48", "grey", "grey48",  "deepskyblue3")
cu.alpha = sapply(cu, col.alpha, alpha = 0.4) %>% unname()
mtheme = list(
  geom_boxplot(show.legend = FALSE, outlier.shape = NA), 
  theme_bw(base_size = 10.5), 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)), 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold",family = "Helvetica"), 
        axis.text.y =  element_text(size = 6), 
        axis.title.y = element_text(size = 10))
  )

cua = sapply(c('dodgerblue', 'red'), col.alpha, 0.2) %>% unname()
cub = c('dodgerblue', 'red')
baselinetheme = list(
  mtheme, 
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, color = 'black', size = 10), 
        axis.title.y = element_text(size = 14)), 
  geom_boxplot(show.legend = FALSE, outlier.shape = 21),
  scale_y_continuous( breaks= scales::pretty_breaks(), expand = c(0.15,0)), 
  scale_x_discrete(labels = c("low", "high")), 
  xlab("Antibody Response"), 
  scale_fill_manual(values = cua), 
  scale_color_manual(values = cub) 
)

# baseline bulk theme 
cu1 = sapply(c('dodgerblue','grey', 'red'), col.alpha, 0.2) %>% unname()
cu2 = c('dodgerblue','grey', 'red')
bulktheme = list(
  mtheme, 
  theme(axis.text.x=element_text(angle = 0, hjust = 0.5, color = 'black', size = 10), 
        axis.title.y = element_text(size = 14)), 
  geom_boxplot(show.legend = FALSE, outlier.shape = 21),
  scale_y_continuous( breaks= scales::pretty_breaks(), expand = c(0.15,0)), 
  scale_x_discrete(labels = c("low", "mid", "high")), 
  xlab("Antibody Response"), 
  scale_fill_manual(values = cu1), 
  scale_color_manual(values = cu2) 
)
my_compare = list(l1 = c("1", "0"), l2 = c("2", "0"))

# define adjuvant signatures 
# leading edge combined signature 
# upregulated gnees only 
# mono and mDC specific and combined genes
li.up = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li.up.rds'))
li2.up = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li2.up.rds'))
li.full = unlist(li.up,use.names = FALSE, recursive = TRUE) %>%  unique 
mono.genes  = li.up$CD14_Mono %>%  unlist() %>% unique()
mdc.genes = li2.up$mDC %>%  unlist() %>% unique()
as03.sig = list('AS03_signature' = li.full, "AS03_Monocyte" = mono.genes, 'AS03_mDC' = mdc.genes)


# leading edge signatures from validation cohort
mono.va.le = readRDS(file = here('mid_res/nat_adj/generated_data/V4/mono.va.le.rds'))
mono.va.le = list('AS03_Monocyte_LE' = mono.va.le)
dc.va.le = readRDS(file = here('mid_res/nat_adj/generated_data/V4/dc.va.le.rds'))
dc.va.le = list('AS03_mDC_LE' = dc.va.le)
as03.sig = c(as03.sig, mono.va.le, dc.va.le)


# load average day 1 comparison cohort data 
av_tidy = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gene_dist/av_tidy.rds'))

# load average baseline high  and low responders from unadjuvanted cohort 
av0 = readRDS(file = here('mid_res/nat_adj/generated_data/V4/av0.rds'))

# load bulk microarray data and process to baseline samples 
array = data.table::fread("data/CHI_H1N1_data/microarray/CHI_GE_matrix_gene.txt", data.table = F) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("gene") %>% 
  select(., matches("day0")) %>% 
  select(-matches("day70")) %>% 
  select(-matches("pre")) %>%
  select(-matches("day1"))

#############################
# CD14 Monocytes 
#############################
# Monocyte combined AS03 signature average across time between groups 
#natural adjuvant test
mono.sig.av2 = 
  av0$CD14_Mono %>% 
  filter(gene %in% as03.sig$AS03_Monocyte_LE) %>% 
  group_by(sample, response) %>% 
  summarize(meansig = mean(count))
#plot
p = ggplot(mono.sig.av2, aes(x = response, y = meansig, fill = response , color = response)) +
  baselinetheme  +
  ggpubr::stat_compare_means(method = 'wilcox',paired = FALSE, show.legend = FALSE) + 
  ylab('Monocyte AS03 Adjuvant Signature') +
  ggtitle('Baseline: CD14 Monocytes\nunadjuvanted subjects')
p
ggsave(p, filename = paste0(figpath, 'baseline_monosig_monocytes_LeadingEdgeVand.pdf'), width = 2.7, height = 4)

#############################
# mDC
#############################

# mdc Combined AS03 signature average across time between groups 
mdc.sig.av = 
  av_tidy$mDC %>% 
  filter(gene %in% mdc.genes) %>% 
  group_by(sample, group) %>% 
  summarize(meansig = mean(count))
#plot
p = ggplot(mdc.sig.av, aes(x = group, y = meansig, fill = group , color = group)) +
  mtheme + 
  theme(axis.title.x = element_blank()) +
  ylab('mDC AS03 Adjuvant Signature') +
  scale_color_manual(values = cu) +
  scale_fill_manual(values = cu.alpha) + 
  ggtitle('mDC')
ggsave(p,filename = paste0(figpath, 'as03_mDC_sig.pdf'), width = 1.9, height = 3)


# natural adjuvant test - mDC
mdc.sig.av2 = 
  av0$mDC %>% 
  filter(gene %in% as03.sig$AS03_mDC) %>% 
  group_by(sample, response) %>% 
  summarize(meansig = mean(count))
#plot
p = ggplot(mdc.sig.av2, aes(x = response, y = meansig, fill = response , color = response)) +
  baselinetheme  +
  ggpubr::stat_compare_means(method = 'wilcox',paired = FALSE, show.legend = FALSE) + 
  ylab('mDC AS03 Adjuvant Signature') +
  ggtitle('Baseline: mDCs\nunadjuvanted subjects')
p
ggsave(p, filename = paste0(figpath, 'baseline_mDC_mDC_LeadingEdgeVand.pdf'), width = 2.7, height = 4)


for (i in 1:length(as03.sig)) {
  data.table::fwrite(as03.sig[i],file = paste0(datapath, names(as03.sig)[i],'.txt'))
}
```

### Fig.5. Analysis of cell frequency of activated monocyte phenotypes in flow cytometry data <a name="fig5.2"></a>

Analyze flow cytometry data for differences at baseline in high and low
responders. Only test innate cell subsets to focus on hypothesis
generated from analysis of CITE-seq data. Further test the identified
activated monocyte phenotype for its longitudinal kinetics and day 1 vs
baseline fold change difference in high and low responders.  
mid_res/flow_kinetic/flow_cellfreq_kinetics.r

``` r
# flow activated monocyte kinetic 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(emmeans))
source(here('functions/MattPMutils.r'))

# save path 
figpath = here('mid_res/flow_kinetic/figures/'); 
dir.create(figpath,recursive = TRUE)

# test only innate subsets based on hypothesis generated from CITE-seq data
# load flow data 
fp = data.table::fread('data/CHI_H1N1_data/flow/flow_annotation.txt', header = TRUE)
id.select = paste('ID', 64:78,sep = '') %>%  as.character()
fp = fp %>% filter(`Population ID` %in% id.select)

# load flow data day 1 fold changes 
fd = 
  data.table::fread(file = here('data/CHI_H1N1_data/flow/day1-day0.log10.txt'),header = TRUE) %>% 
  filter(ID %in% fp$`Population ID`) %>%  
  mutate(subset_name = plyr::mapvalues(ID, from = fp$`Population ID`,to = fp$`Subset name`)) %>%
  select(-ID) %>% 
  column_to_rownames('subset_name') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Subject')


# map titers to data 
titer = data.table::fread(file = here('data/CHI_H1N1_data/titer/titer_processed.txt'))
fd$adjMFC_class = plyr::mapvalues(x = fd$Subject, from = titer$Subject,to = titer$adjMFC_class )
fd = fd %>% select(adjMFC_class, everything())
fd$adjMFC_class = factor(fd$adjMFC_class, levels = c('0','1','2'))
fd = fd[!is.na(fd$adjMFC_class), ]

# baseline 
# these are raw percentages so use non parametric rank stats 
fd3 = 
  data.table::fread(file = here('data/CHI_H1N1_data/flow/day0.raw.txt'),header = TRUE) %>% 
  filter(ID %in% fp$`Population ID`) %>%  
  mutate(subset_name = plyr::mapvalues(ID, from = fp$`Population ID`,to = fp$`Subset name`)) %>%
  select(-ID) %>% 
  column_to_rownames('subset_name') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Subject')
# map adj mfc class 
fd3$adjMFC_class = plyr::mapvalues(x = fd3$Subject, from = titer$Subject,to = titer$adjMFC_class )
fd3 = fd3 %>% select(adjMFC_class, everything())
fd3$adjMFC_class = factor(fd3$adjMFC_class, levels = c('0','1','2'))
fd3 = fd3[!is.na(fd3$adjMFC_class), ]
fd3 = fd3[!fd3$adjMFC_class == '1', ]

wilcox.res3 = apply(
  X =  fd3[, 3:ncol(fd3)],
  MARGIN = 2,
  FUN = function(x) {
    wilcox.test(x ~ fd3$adjMFC_class) %>%  broom::tidy()
  }) %>% 
  bind_rows(.id = 'subset')
wilcox.res3 %>%  filter(p.value < 0.1)
saveRDS(wilcox.res3,file = paste0(datapath,'wilcox.res3.rds'))


# comparison 
flow_compare = list(c('2','0'))

# color specification
cu1 = sapply(c('dodgerblue', 'red'), col.alpha, 0.2) %>% unname()
cu2 = c('dodgerblue', 'red')

# theme 
mtheme = list(
  theme_bw(), 
  geom_boxplot(show.legend = FALSE, outlier.shape = 21), 
  ggpubr::stat_compare_means(comparisons = flow_compare,method = 'wilcox', paired = FALSE),
  scale_y_continuous( breaks= scales::pretty_breaks(), expand = c(0.15,0)), 
  theme(axis.text.y = element_text(size = 8), 
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = 'black')),
  scale_x_discrete(labels = c("low", "high")),
  xlab("Antibody \n Response"),
  scale_fill_manual(values = cu1),
  scale_color_manual(values = cu2) 
)

# plot 
p3 = ggplot(fd3, aes(x = adjMFC_class, y = `activated monocyte HLA-DR+`,
                     color = adjMFC_class,
                     fill = adjMFC_class)) +  
  mtheme
ggsave(p3,filename = paste0(figpath, 'mono_HLADR.pdf'), width = 2, height = 3)



# longitudinal analysis of cell population frequency
# connection between baseline (-7 and 0) and day 1 kinetics 
d01 = data.table::fread(file = here('data/CHI_H1N1_data/flow/longitudinal/pre7.raw.txt'), header = TRUE)
d02 = data.table::fread(file = here('data/CHI_H1N1_data/flow/longitudinal/day0.raw.txt'), header = TRUE)
d1 = data.table::fread(file = here('data/CHI_H1N1_data/flow/longitudinal/day1.raw.txt'), header = TRUE)
d7 = data.table::fread(file = here('data/CHI_H1N1_data/flow/longitudinal/day7.raw.txt'), header = TRUE)
d70 = data.table::fread(file = here('data/CHI_H1N1_data/flow/longitudinal/day70.raw.txt'), header = TRUE)

d.list = list('t01' = d01, 't02' = d02, 't1' = d1, 't7' = d7, 't70' = d70 )

# format, combine and label by time 
d = lapply(d.list, function(x) 
  x %>%  
    filter(ID %in% fp$`Population ID`) %>%  
    mutate(subset_name = plyr::mapvalues(ID, from = fp$`Population ID`,to = fp$`Subset name`)) %>%
    select(-ID) %>% 
    column_to_rownames('subset_name') %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column('Subject') ) %>% 
  bind_rows(.id = 'timepoint')
# map adj mfc class 
d$adjMFC_class = plyr::mapvalues(x = d$Subject, from = titer$Subject,to = titer$adjMFC_class )
d = d %>%  filter(adjMFC_class %in% c('0','2'))
d$adjMFC_class = factor(d$adjMFC_class, levels = c('0','2'))
d$timepoint = factor(d$timepoint, levels = c("t01", "t02", "t1",  "t7",  "t70"))

# plot time course 
mtheme2 = list(
  theme_bw(), 
  theme(axis.text.y = element_text(size = 8), 
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = 'black')),
  scale_fill_manual(values = cu1),
  scale_color_manual(values = cu2) 
)

# monocyte
d2 = d[!is.na(d$`activated monocyte HLA-DR+`), ]
p3=
  ggplot(d2, aes(x = timepoint, y = `activated monocyte HLA-DR+`, color = adjMFC_class, group = Subject)) + 
  geom_line(size = 0.5, alpha = 0.2, show.legend = FALSE) + 
  geom_point(size = 0.5, alpha = 0.2, shape = 21, show.legend = FALSE) + 
  geom_smooth(data = d2, size  = 2.5, 
              method = 'loess',
              aes(x = timepoint,  y = `activated monocyte HLA-DR+`, 
                  color = adjMFC_class,
                  fill = adjMFC_class,
                  group = adjMFC_class), alpha = 0.2,
              se = TRUE,  show.legend = FALSE) +
  scale_x_discrete(expand = c(0,0.1)) +
  mtheme2 
p3
ggsave(p3,filename = paste0(figpath, 'monoDRfreqtime.pdf'), width = 3, height = 3)

### mixed model 
d2 = d %>% 
  select(timepoint, subjectid = Subject,
         drmono = `activated monocyte HLA-DR+`,
         adjmfc = adjMFC_class) %>% 
  mutate(timepoint= as.character(timepoint)) %>% 
  mutate(timepoint = ifelse(timepoint %in% c('t01','t02'),yes =  't0',no =  timepoint)) %>% 
  mutate(timepoint = factor(timepoint , levels = c( "t0",  "t1" , "t7"  ,"t70")))

# fit model 
m = lme4::lmer(drmono ~ timepoint + (1|subjectid),data = d2)
m2 = lme4::lmer(drmono ~ timepoint*adjmfc + (1|subjectid),data = d2)
anova(m,m2)
# Data: d2
# Models:
#   m: drmono ~ timepoint + (1 | subjectid)
# m2: drmono ~ timepoint * adjmfc + (1 | subjectid)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# m     6 615.95 632.58 -301.98   603.95                       
# m2   10 615.02 642.73 -297.51   595.02 8.9311  4    0.06284 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
emm1 = emmeans(m2, revpairwise ~ timepoint|adjmfc)
p = emm1$contrasts[c(1,7), ] %>%  plot() + theme_bw() + 
  xlab('effect size') + 
  ylab('group') + 
  theme(text = element_text(size = 5))
p
ggsave(p,filename = paste0(figpath, 'monoDREMM.pdf'), width = 3, height = 1)


sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] emmeans_1.5.4   scglmmr_0.1.0   forcats_0.5.1   stringr_1.4.0   dplyr_1.0.4     purrr_0.3.4     readr_1.4.0     tidyr_1.1.2    
# [9] tibble_3.1.8    ggplot2_3.3.3   tidyverse_1.3.0 here_1.0.1     
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.2                  tidyselect_1.2.0            lme4_1.1-26                 htmlwidgets_1.5.3          
# [5] RSQLite_2.2.7               AnnotationDbi_1.52.0        grid_4.0.5                  BiocParallel_1.24.1        
# [9] scatterpie_0.1.7            munsell_0.5.0               codetools_0.2-18            statmod_1.4.35             
# [13] withr_2.4.3                 colorspace_2.0-0            GOSemSim_2.16.1             Biobase_2.50.0             
# [17] knitr_1.39                  rstudioapi_0.13             stats4_4.0.5                ggsignif_0.6.0             
# [21] DOSE_3.16.0                 labeling_0.4.2              MatrixGenerics_1.2.1        Rdpack_2.1.1               
# [25] GenomeInfoDbData_1.2.4      polyclip_1.10-0             bit64_4.0.5                 farver_2.0.3               
# [29] pheatmap_1.0.12             rprojroot_2.0.2             downloader_0.4              coda_0.19-4                
# [33] vctrs_0.5.1                 generics_0.1.2              TH.data_1.0-10              xfun_0.30                  
# [37] R6_2.5.0                    doParallel_1.0.16           GenomeInfoDb_1.26.7         graphlayouts_0.7.2         
# [41] locfit_1.5-9.4              bitops_1.0-6                cachem_1.0.4                fgsea_1.16.0               
# [45] DelayedArray_0.16.3         assertthat_0.2.1            scales_1.1.1                nnet_7.3-15                
# [49] multcomp_1.4-16             ggraph_2.0.5                enrichplot_1.10.2           gtable_0.3.0               
# [53] egg_0.4.5                   tidygraph_1.2.0             sandwich_3.0-0              rlang_1.0.6                
# [57] slanter_0.2-0               splines_4.0.5               rstatix_0.7.0               checkmate_2.0.0            
# [61] broom_0.7.5                 BiocManager_1.30.10         reshape2_1.4.4              abind_1.4-5                
# [65] modelr_0.1.8                backports_1.2.1             Hmisc_4.5-0                 qvalue_2.22.0              
# [69] clusterProfiler_3.18.1      tools_4.0.5                 ellipsis_0.3.2              gplots_3.1.1               
# [73] RColorBrewer_1.1-2          BiocGenerics_0.36.1         Rcpp_1.0.9                  plyr_1.8.6                 
# [77] base64enc_0.1-3             progress_1.2.2              zlibbioc_1.36.0             RCurl_1.98-1.3             
# [81] prettyunits_1.1.1           rpart_4.1-15                ggpubr_0.4.0                viridis_0.5.1              
# [85] cowplot_1.1.1               S4Vectors_0.28.1            zoo_1.8-8                   cluster_2.1.2              
# [89] SummarizedExperiment_1.20.0 haven_2.4.3                 ggrepel_0.9.1               fs_1.5.0                   
# [93] variancePartition_1.25.6    magrittr_2.0.3              data.table_1.14.0           DO.db_2.9                  
# [97] openxlsx_4.2.3              reprex_1.0.0                mvtnorm_1.1-1               packrat_0.7.0              
# [101] matrixStats_0.58.0          hms_1.0.0                   GSVA_1.38.2                 xtable_1.8-4               
# [105] pbkrtest_0.5-0.1            RhpcBLASctl_0.21-247.1      XML_3.99-0.6                jpeg_0.1-8.1               
# [109] rio_0.5.16                  readxl_1.3.1                IRanges_2.24.1              gridExtra_2.3              
# [113] compiler_4.0.5              KernSmooth_2.23-18          crayon_1.4.1                shadowtext_0.0.9           
# [117] htmltools_0.5.2             minqa_1.2.4                 mgcv_1.8-34                 ggfun_0.0.4                
# [121] Formula_1.2-4               lubridate_1.8.0             DBI_1.1.1                   corrplot_0.84              
# [125] tweenr_1.0.2                dbplyr_2.1.0                MASS_7.3-53.1               boot_1.3-27                
# [129] Matrix_1.4-1                car_3.0-10                  cli_3.4.1                   rbibutils_2.0              
# [133] parallel_4.0.5              igraph_1.2.6                GenomicRanges_1.42.0        pkgconfig_2.0.3            
# [137] rvcheck_0.1.8               foreign_0.8-81              xml2_1.3.2                  foreach_1.5.1              
# [141] annotate_1.68.0             XVector_0.30.0              GeneOverlap_1.26.0          estimability_1.3           
# [145] rvest_0.3.6                 digest_0.6.27               graph_1.68.0                cellranger_1.1.0           
# [149] fastmatch_1.1-0             htmlTable_2.1.0             edgeR_3.32.1                GSEABase_1.52.1            
# [153] curl_4.3                    gtools_3.8.2                nloptr_1.2.2.2              lifecycle_1.0.3            
# [157] nlme_3.1-152                jsonlite_1.7.2              aod_1.3.1                   carData_3.0-4              
# [161] viridisLite_0.3.0           limma_3.46.0                fansi_0.4.2                 pillar_1.8.1               
# [165] lattice_0.20-41             fastmap_1.1.0               httr_1.4.2                  survival_3.2-10            
# [169] GO.db_3.12.1                glue_1.6.2                  UpSetR_1.4.0                zip_2.1.1                  
# [173] png_0.1-7                   iterators_1.0.13            bit_4.0.4                   ggforce_0.3.3              
# [177] stringi_1.5.3               blob_1.2.1                  org.Hs.eg.db_3.12.0         latticeExtra_0.6-29        
# [181] caTools_1.18.1              memoise_2.0.0   
```

### Fig.5. Analysis of CyTOF stimulation phenotypes <a name="fig5.3"></a>

Visualize stimulated and unstimulated cells identified by HDStim.  
mid_res/stim/visualize_stim_cells.r

``` r
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(HDStIM))
library(Rcpp)
library(emmeans)
library(lme4)
source('functions/MattPMutils.r')
# save paths 
figpath = here('mid_res/stim/figures/')
datapath = here('mid_res/stim/generated_data/')
dir.create(figpath); dir.create(datapath)

# read data from stim cell selector 
d = readRDS(file = here('data/stim/mapped_data.rds'))

d$umap_plot_data  

up = HDStIM::plot_umap(mapped_data = d)
pd = up[[4]]
pd2 = pd$data
pd2$stim = ifelse(str_sub(pd2$response_status, -5,-1) == 'Stim.', yes = 'LPS stimulated', no = 'unstimulated')
cf = ggsci::pal_d3( palette = 'category20', 0.5)(2) %>% rev
p = 
ggplot(pd2, aes(x = UMAP1, y = UMAP2, fill = stim )) + 
  geom_point(shape = 21, stroke = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = cf)  +
  ggtitle(pd$labels$title)
ggsave(p,filename = paste0(figpath, 'umap_lps_mono.png'), width = 5, height = 4)
ggsave(p,filename = paste0(figpath, 'umap_lps_mono_outline.pdf'), width = 5, height = 4)
```

Fit mixed effects model of median phospho protein marker expression in
classical monocytes pre vs post stimulation and compare effects in high
vs low responders.  
mid_res/cytof_stim/stim_test_ag.r

``` r
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(HDStIM))
suppressMessages(library(Rcpp))
suppressMessages(library(emmeans))
suppressMessages(library(lme4))
source('functions/MattPMutils.r')
# save paths 
figpath = here('mid_res/cytof_stim/figures/')
datapath = here('mid_res/cytof_stim/generated_data/')
dir.create(figpath); dir.create(datapath)

# baseline bulk theme 
cu1 = sapply(c('dodgerblue','red'), col.alpha, 0.2) %>% unname()
cu2 = c('dodgerblue', 'red')

# read data from stim cell selector 
d = readRDS(file = here('data/stim/mapped_data.rds'))

# celltype markers 
cmarkers = c('CD45',  'CD7',  'CD19',  'CD4', 'IgD', 'CD20', 
             "CD11c", "CD127", 'CD25', 'CD123', 'CD27', 'CD24', 
             'CD14', 'CD56', 'CD16', 'CD38', 'CD8', 'CD45RA',
             'CD3', 'HLA_DR')

# phenotyping markers 
pmarkers = c("pPLCg2", "pSTAT5", "AKT", "pSTAT1", "pP38",
             "pSTAT3", "IkBa", "pCREB", "pERK1_2", "pS6")

# innate cells 
innate.cells = c('CD14Mono', 'DC1', 'DC2')

########################
# make a data matrix 
dr = d$response_mapping_main
dr$cell_id = paste(dr$sample_id, rownames(dr),sep = '_')

dr = as.data.frame(dr) %>%  
  column_to_rownames('cell_id') %>%  
  mutate(sx = paste(patient_id, condition, stim_type, sep = '_')) %>% 
  mutate(response.stim = paste(condition, stim_type, sep = '_')) %>% 
  mutate(adjMFC = factor(condition,levels = c('high', 'low'))) 

# log transform the state marker matrix
mat = dr %>% select(all_of(pmarkers))
#mat = log(mat + 1 )

#make a separate dataframe of metadata 
prots = c(cmarkers, pmarkers)
met = dr %>% select(!all_of(prots))

# combine log transformed data back with md 
stopifnot(isTRUE(all.equal(rownames(met), rownames(mat))))
dr = cbind(met, mat )

# aggregate the protein markers to test the median log transformed marker intensity 
da = dr %>% 
  group_by(sx, response.stim, patient_id, batch, cell_population, stim_type) %>%  
  summarize_at(.vars = pmarkers, .funs = median) %>% 
  filter(cell_population %in% c('CD14Mono', 'DC1', 'DC2', 'CD16Mono'))


# model 
f1 = as.formula(
  prot ~ 0 + batch + response.stim + ( 1 | patient_id)
)

# separate by stim 
lps = da %>% 
  filter(stim_type %in% c('U', 'L')) %>% 
  filter(cell_population == 'CD14Mono')

# modify group factor for contrast estimand 
lps$response.stim = factor(
  lps$response.stim, 
  levels = c("high_U", "high_L", "low_U", "low_L")
)

dfit = lps 
contrast_sub = c ( '(high_L - high_U) - (low_L - low_U)' )
emmeans =  reslist =  list()
for (i in 1:length(pmarkers)) {
  
  # test 
  prot.name = pmarkers[i]
  print(prot.name)
  # # metadata 
  met = dfit %>% select(!all_of(pmarkers))
  
  # extract single protein 
  dat_vec = data.frame(dfit[ ,prot.name])
  colnames(dat_vec) = 'prot'
  # model data to fit 
  dat_fit = base::data.frame(cbind(dat_vec, met))
  
  # save quick plot 
  dplot = dat_fit %>% filter(response.stim %in% c('low_L', 'high_L'))
  dplot$response.stim = factor(dplot$response.stim, levels = c('low_L', 'high_L'))
  p = ggplot(dplot, aes(x = response.stim,y = prot , color = response.stim, fill = response.stim)) +
    theme_bw() + 
    geom_boxplot(show.legend = FALSE, outlier.shape = NA) + 
    geom_jitter(show.legend = FALSE, width = 0.2, shape = 21, color = 'black', stroke = 0.3) + 
    ylab(prot.name) + 
    scale_fill_manual(values = cu1) + 
    scale_color_manual(values = cu2) 
  p
  ggsave(p, filename = paste0(figpath,prot.name,'lps.pdf'), width = 1.75, height = 2)
  
  # fit models 
  m1 = tryCatch(
    lme4::lmer(formula = f1, data = dat_fit),
    error = function(e)
      return(NA)
  )
  
  # emmeans to apply contrast of fold changes 
  emm1 = tryCatch(
    emmeans::emmeans(
      object = m1, 
      specs = ~ 'response.stim', 
      data = dat_fit, 
      lmer.df = "asymptotic"),
    error = function(e) return(NA)
  )

  # apply contrasts 
  if (!is.na(emm1)) {
    m1.cont = contrast(emm1, method = 'revpairwise',adjust = NULL)
    m1cont = pairs(m1.cont,adjust = NULL) %>% 
      as.data.frame() %>% 
      filter( contrast == '(high_L - high_U) - (low_L - low_U)' )
    
    
    # store results 
    emmeans[[i]] = emm1 
    names(emmeans)[i] = prot.name
    
    # store contrast test 
    tidy.fit = m1cont %>%  
      mutate(prot = prot.name) %>%
      select(prot, everything())
    reslist[[i]] = tidy.fit
  } else{ 
    # store results 
    emmeans[[i]] = emm1 
    names(emmeans)[i] = prot.name
    reslist[[i]] = data.frame(
      prot = prot.name,
      contrast = NA,
      estimate = NA,
      SE = NA,
      df = NA,
      z.ratio = NA,
      p.value = NA
    )
    }
}
rd = do.call(rbind,reslist)
rd$stim = 'LPS'
data.table::fwrite(rd,file = paste0(datapath,'mono14_lps.txt'), sep = '\t')


######################
# PMA
 pma = da %>%
  filter(stim_type %in% c('U', 'P')) %>%
  filter(cell_population == 'CD14Mono')

# modify group factor for contrast estimand 
pma$response.stim = factor(
  pma$response.stim, 
  levels = c("high_U", "high_P", "low_U", "low_P")
)
dfit = pma 
contrast_sub = c ( '(high_P - high_U) - (low_P - low_U)' )
emmeans =  reslist =  list()
for (i in 1:length(pmarkers)) {
  #i = 5
  # test 
  prot.name = pmarkers[i]
  print(prot.name)
  # # metadata 
  met = dfit %>% select(!all_of(pmarkers))
  
  # extract single protein 
  dat_vec = data.frame(dfit[ ,prot.name])
  colnames(dat_vec) = 'prot'
  # model data to fit 
  dat_fit = base::data.frame(cbind(dat_vec, met))
  
  # save quick plot 
  dplot = dat_fit %>% filter(response.stim %in% c('low_P', 'high_P'))
  dplot$response.stim = factor(dplot$response.stim, levels = c('low_P', 'high_P'))
  p = ggplot(dplot, aes(x = response.stim,y = prot , color = response.stim, fill = response.stim)) +
    theme_bw() + 
    geom_boxplot(show.legend = FALSE, outlier.shape = NA) + 
    geom_jitter(show.legend = FALSE, width = 0.2, shape = 21, color = 'black', stroke = 0.3) + 
    ylab(prot.name) + 
    scale_fill_manual(values = cu1) + 
    scale_color_manual(values = cu2) 
  p
  ggsave(p, filename = paste0(figpath,prot.name,'pma.pdf'), width = 1.75, height = 2)
  
  # fit models 
  m1 = tryCatch(
    lme4::lmer(formula = f1, data = dat_fit),
    error = function(e)
      return(NA)
  )
  
  # emmeans to apply contrast of fold changes 
  emm1 = tryCatch(
    emmeans::emmeans(
      object = m1,
      specs = ~ 'response.stim',
      data = dat_fit,
      lmer.df = "asymptotic"
    ), 
    error = function(e) return(NA)
  )
  
  # apply contrasts 
  if (!is.na(emm1)) {
    m1.cont = contrast(emm1, method = 'revpairwise',adjust = NULL)
    m1cont = pairs(m1.cont,adjust = NULL) %>% 
      as.data.frame() %>% 
      filter( contrast == '(high_P - high_U) - (low_P - low_U)' )
    
    
    # store results 
    emmeans[[i]] = emm1 
    names(emmeans)[i] = prot.name
    
    # store contrast test 
    tidy.fit = m1cont %>%  
      mutate(prot = prot.name) %>%
      select(prot, everything())
    reslist[[i]] = tidy.fit
  } else{ 
    # store results 
    emmeans[[i]] = emm1 
    names(emmeans)[i] = prot.name
    reslist[[i]] = data.frame(
      prot = prot.name,
      contrast = NA,
      estimate = NA,
      SE = NA,
      df = NA,
      z.ratio = NA,
      p.value = NA
    )
  }
}
rd = do.call(rbind,reslist)
rd$stim = 'PMA'
data.table::fwrite(rd,file = paste0(datapath,'mono14_PMA.txt'), sep = '\t')

# IFN
# separate by stim 
ifn = da %>% 
  filter(stim_type %in% c('U', 'A')) %>% 
  filter(cell_population == 'CD14Mono')


# modify group factor for contrast estimand 
ifn$response.stim = factor(
  ifn$response.stim, 
  levels = c("high_U", "high_A", "low_U", "low_A")
)
dfit = ifn 
contrast_sub = c ( '(high_A - high_U) - (low_A - low_U)' )
emmeans =  reslist =  list()
for (i in 1:length(pmarkers)) {
  
  # test 
  prot.name = pmarkers[i]
  print(prot.name)
  # # metadata 
  met = dfit %>% select(!all_of(pmarkers))
  
  # extract single protein 
  dat_vec = data.frame(dfit[ ,prot.name])
  colnames(dat_vec) = 'prot'
  # model data to fit 
  dat_fit = base::data.frame(cbind(dat_vec, met))
  
  # save quick plot 
  dplot = dat_fit %>% filter(response.stim %in% c('low_A', 'high_A'))
  dplot$response.stim = factor(dplot$response.stim, levels = c('low_A', 'high_A'))
  p = ggplot(dplot, aes(x = response.stim,y = prot , color = response.stim, fill = response.stim)) +
    theme_bw() + 
    geom_boxplot(show.legend = FALSE, outlier.shape = NA) + 
    geom_jitter(show.legend = FALSE, width = 0.2, shape = 21, color = 'black', stroke = 0.3) + 
    ylab(prot.name) + 
    scale_fill_manual(values = cu1) + 
    scale_color_manual(values = cu2) 
  ggsave(p, filename = paste0(figpath,prot.name,'IFN.pdf'), width = 1.75, height = 2)
  
  # fit models
  m1 = tryCatch(
    lme4::lmer(formula = f1, data = dat_fit),
    error = function(e)
      return(NA)
  )
  
  # emmeans to apply contrast of fold changes
  emm1 = tryCatch(
    emmeans::emmeans(
      object = m1,
      specs = ~ 'response.stim',
      data = dat_fit,
      lmer.df = "asymptotic"
    ),
    error = function(e)
      return(NA)
  )
  
  # apply contrasts 
  if (!is.na(emm1)) {
    m1.cont = contrast(emm1, method = 'revpairwise',adjust = NULL)
    m1cont = pairs(m1.cont,adjust = NULL) %>% 
      as.data.frame() %>% 
      filter( contrast == '(high_A - high_U) - (low_A - low_U)' )
    
    
    # store results 
    emmeans[[i]] = emm1 
    names(emmeans)[i] = prot.name
    
    # store contrast test 
    tidy.fit = m1cont %>%  
      mutate(prot = prot.name) %>%
      select(prot, everything())
    reslist[[i]] = tidy.fit
  } else{ 
    # store results 
    emmeans[[i]] = emm1 
    names(emmeans)[i] = prot.name
    reslist[[i]] = data.frame(
      prot = prot.name,
      contrast = NA,
      estimate = NA,
      SE = NA,
      df = NA,
      z.ratio = NA,
      p.value = NA
    )
  }
}
rd = do.call(rbind,reslist)
rd$stim = 'IFN'
data.table::fwrite(rd,file = paste0(datapath,'mono14_IFN.txt'), sep = '\t')

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] lme4_1.1-26     Matrix_1.4-1    emmeans_1.5.4   Rcpp_1.0.9      HDStIM_0.1.0   
# [6] here_1.0.1      forcats_0.5.1   stringr_1.4.0   dplyr_1.0.4     purrr_0.3.4    
# [11] readr_1.4.0     tidyr_1.1.2     tibble_3.1.8    ggplot2_3.3.3   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
# [1] httr_1.4.2        jsonlite_1.7.2    splines_4.0.5     modelr_0.1.8      assertthat_0.2.1 
# [6] statmod_1.4.35    cellranger_1.1.0  yaml_2.2.1        pillar_1.8.1      backports_1.2.1  
# [11] lattice_0.20-41   glue_1.6.2        digest_0.6.27     rvest_0.3.6       minqa_1.2.4      
# [16] colorspace_2.0-0  sandwich_3.0-0    htmltools_0.5.2   plyr_1.8.6        pkgconfig_2.0.3  
# [21] broom_0.7.5       haven_2.4.3       xtable_1.8-4      mvtnorm_1.1-1     scales_1.1.1     
# [26] farver_2.0.3      generics_0.1.2    ellipsis_0.3.2    TH.data_1.0-10    withr_2.4.3      
# [31] Boruta_7.0.0      cli_3.4.1         survival_3.2-10   magrittr_2.0.3    crayon_1.4.1     
# [36] readxl_1.3.1      evaluate_0.15     estimability_1.3  fs_1.5.0          fansi_0.4.2      
# [41] nlme_3.1-152      MASS_7.3-53.1     xml2_1.3.2        rsconnect_0.8.25  tools_4.0.5      
# [46] data.table_1.14.0 hms_1.0.0         lifecycle_1.0.3   multcomp_1.4-16   munsell_0.5.0    
# [51] reprex_1.0.0      packrat_0.7.0     compiler_4.0.5    rlang_1.0.6       grid_4.0.5       
# [56] nloptr_1.2.2.2    ggridges_0.5.3    rstudioapi_0.13   rmarkdown_2.9     labeling_0.4.2   
# [61] boot_1.3-27       gtable_0.3.0      codetools_0.2-18  DBI_1.1.1         R6_2.5.0         
# [66] zoo_1.8-8         lubridate_1.8.0   knitr_1.39        fastmap_1.1.0     uwot_0.1.10      
# [71] utf8_1.2.2        rprojroot_2.0.2   stringi_1.5.3     parallel_4.0.5    vctrs_0.5.1      
# [76] xfun_0.30         dbplyr_2.1.0      tidyselect_1.2.0  coda_0.19-4  
```

### Write output <a name="output"></a>

Write results for supplementary tables.  
mid_res/data_write/Final_script_wite_table_fsc.r

``` r
# make tables 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))

datapath = file.path(here('mid_res/data_write/generated_data/')); 
dir.create(datapath)

# specify order of variables in the output table for readability 
var.order = c('contrast', 'celltype', 'pathway', 'NES', 'padj', 'leadingEdge')

# result format for gsea results
filter.gsea = function(list){
  lapply(list, function(x) x %>%  filter(padj < 0.05))
}

format.result = function(x) { 
  x %>% 
    select(all_of(var.order), everything()) %>%
    arrange(celltype, NES) %>% 
    tibble::remove_rownames()
}



# Baseline
# res text gsea curate
# gsea res raw 
g0.sub = readRDS(file = here("mid_res/baseline_response/dataV3/g0.sub.rds"))
g0.sub = do.call(rbind, g0.sub)
g0.sub$contrast = 'baseline high vs low responders'
d0.res = format.result(g0.sub) %>%
  select(-c(signal)) %>% 
  mutate(model = 'gene ~ 0 + response.group + batch + sex + age')
  

# day 1 non-adjuvanted vaccine 
g1c = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/g1c.rds'))
g1c = do.call(rbind, g1c)
g1c$contrast = '24h vs baseline unadjuvanted vaccine'
g1c$model = 'gene ~ 0 + timepoint + batch + sex + age + (1|subjectid) '

# day 7 non-adjuvanted vaccine
g7f = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/g7f.rds'))
g7f = lapply(g7f, function(x) 
  x %>%  
    filter(!str_sub(pathway, 1,5) == 'REACT' ) %>% 
    filter(NES > 0) %>%  
    filter(pval <0.1)
)
g7f = do.call(rbind, g7f)
g7f$contrast = 'day 7 vs baseline unadjuvanted vaccine'
g7f$model = 'gene ~ 0 + timepoint + batch + sex + age + (1|subjectid) '

# as03 model 
gc = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/gc.rds'))
gc = lapply(gc, function(x) x %>%  filter(padj < 0.1))
gc = do.call(rbind, gc)
gc = format.result(gc)
gc$contrast = 'AS03 vs unadjuvanted vaccine day 1 vs baseline fold change difference'
gc = format.result(gc) %>% 
  mutate(model = '0 + timepoint_vaccinegroup + sex + age')

# combine 
d = rbind(g1c, g7f, gc, d0.res)

# write results 
data.table::fwrite(d,
                   file = paste0(datapath,'combined.results.fsc.txt'), 
                   sep = '\t')


# gene signatures 
core_sigs = readRDS(file = here('signature_curation/sig_test_sub.rds'))
mann = data.table::fread(file = here('signature_curation/sig_test_sub_annotation.txt'))
dmod = mann
dmod$signature_genes = core_sigs

# add natural adjuvant signautres   
mv = readRDS(file = here('mid_res/vand/generated_data/mv.rds'))
dcv = readRDS(file = here('mid_res/vand/generated_data/dcv.rds'))
mono.sig = mv$leadingEdge  %>% unlist() %>% unique() 
dc.sig = dcv$leadingEdge  %>% unlist() %>% unique()
as03.sig = list(mono.sig,dc.sig)

dcite= data.frame(
  pathway = c("CD14_Mono_AS03", "mDC_AS03"), 
  annotation = c(rep('CITE-seq contrast model and sorted cell validated', 2))
)
dcite$signature_genes = as03.sig
# combine with signatures 
dmod = rbind(dmod, dcite)

data.table::fwrite(dmod,
                   file = paste0(datapath,'combined.modules.fsc.txt'), 
                   sep = '\t')


### Variance fractions 
vp = readRDS(file = here('mid_res/variance_partition/generated_data/vp.rds'))
vp = as.data.frame(vp) %>% 
  rownames_to_column('gene') %>% 
  select(gene, everything()) %>% 
  mutate(model = '~ age + (1|sex) + (1|subjectid) + (1|celltype) + (1|timepoint) + (1|adjmfc.group) + (1|celltype:timepoint)')
# write
data.table::fwrite(vp,
                   file = paste0(datapath,'variance.partition.across.celltypes.txt'), 
                   sep = '\t')

# per cell type results 
#pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))
dl = list.files(
  path = here('mid_res/variance_partition/generated_data/'),
  pattern = '.rds',
  recursive = TRUE,
  full.names = TRUE
)
dl = dl[-c(15,17)] # remove total bulk 

# get cell type names (file names)
cts = list.files(
  path = here('mid_res/variance_partition/generated_data/'),
  pattern = '.rds',
  recursive = TRUE,
  full.names = FALSE
)
cts = cts[-c(15,17)]
cts = str_replace_all(string = cts,pattern = 'vp.rds', replacement = '')
# read and format variance partition results 
vl = lapply(dl, readRDS)
for (i in 1:length(vl)) {
  vl[[i]] = vl[[i]] %>% 
    as.data.frame() %>% 
    rownames_to_column('gene') %>% 
    mutate(model_fit_within_this_celltype = names(vl)[i]) %>% 
    mutate(model = '~ age + (1|sex)  + (1|subjectid) + (1|timepoint) + (1|adjmfc.group) + (1|timepoint:adjmfc.group)')
}
vp_within = do.call(rbind, vl)
data.table::fwrite(vp_within,
                   file = paste0(datapath,'variance.partition.within.celltypes.txt'), 
                   sep = '\t')
```

### Low-level bioinformatic processing to generate starting data <a name="preprocessing"></a>

**PREPROCESSING WORKFLOW: Multiscale integration of human and single-cell immune
variations reveals unadjuvanted vaccine high responders are naturally
adjuvanted**

The directory `Flu_CITEseq_preprocess` downloaded from Zenodo contains
raw data needed to demultiplex cells and generate the starting data.

`Flu_CITEseq_preprocess` is another self-contained R project with its
own documentation.

**Please see the full documentation and code for these steps in the
separate readme.**

These scripts are used to document the creation of starting single cell
data object `h1h5_annotated_with_meta.rds` from raw HTO ADT and mRNA umi
reads. `h1h5_annotated_with_meta.rds` is the starting data used in the
main analysis directory `mid_res`. A visual guide below shows how the
object is built.

Scripts:

    1_write_bc_for_citeseqcount.r
    1a_RunCITEseqCount.txt
    2_merge10xlanes_hiseq_hto_adt_rna_V8.r
    2a_RUNDEMUXLET.txt
    3_hto_demux_by_batch.r
    3_Multiseq_Demultiplex_by_batch.r
    4_validate_htodmx_w_demuxlet_singlets.r
    4_validate_multiseq_w_demuxlet_singlets.r
    5_save_merged_Seurat_singlet_object.r

Visual guide:

    -Experiment                                                                           
                                                                                           
                 ┌───────────────────────┐   ┌─────────────────────┐                       
                 │ Concentrated antibody │   │    Pooled cells     │                       
                 │  based on titration   │   │      12 donors      │                       
                 │ 83 surface antibodies │   │     24 samples      │                       
                 │  + 4 isotype control  │   │      per batch      │                       
                 └──────────────┬────────┘   └────┬────────────────┘                       
                                │                 │                                        
                                │    ┌─────────┐  │                                        
                                │    │  stain  │  │                                        
                                └───▶│  cells  │──┘                                        
                                     └────┬────┘                                           
                         ┌──────────────┐ │                                                
                         │ load on 10X  │ │                                                
                        ┌┴────┬─────────┴─┘                                                
                        │     │                                                            
                        │     │                                                            
                        ▼     ▼                                                            
                       .─.   .─.   ┌──────────────────────┐                                
                      (   ) (   )  │... x 6 lanes / batch │                                
                       `─'   `─'   └──────────────────────┘                                
                                                                                           
                                 │ ┌──────────────────────┐                                
                                 │ │   ... x 3 batches    │                                
    -Processing                  │ └──────────────────────┘                                
                                 ▼                                                         
                      ┌───────────────────────────────┐                                    
                      │  HIseq 2500 sequencing GEMs   │                                    
                      └───────────────┬───────────────┘                                    
                                      │                1_write_bc_for_citeseqcount.R       
                               ┌──────┘     ┌────────▶   Save barcode "whitelist"          
                               │            │            for ADT alignment / dsb.          
                               │            │                        │                     
                               ▼            │                        ▼                     
                     ┌──────────────────┐   │              ┌──────────────────┐            
                     │Cell Ranger Count │   │              │  Cite-Seq-Count  │            
                     └──────────────────┘   │              └─────────┬────────┘            
                               │            │                   X 18 lanes                 
                          X 18 lanes        │                        │                     
                               │            │                        ▼                     
                               ▼            │                                              
                          cell              │        cell                 cell             
                       ┌────────┐           │     ┌────────┐           ┌────────┐          
                       │ ┌──────┴─┐         │     │ ┌──────┴─┐         │ ┌──────┴─┐        
                   gene│ │        │         │ HTO │ │        │     ADT │ │        │        
                       │ │mRNA UMI│         │     │ │HTO UMI │         │ │ADT UMI │        
                       │ │  .mtx  │─────────┘     │ │  .mtx  │         │ │  .mtx  │        
                       └─┤        │               └─┤        │         └─┤        │        
                         └────────┘                 └────────┘           └────────┘        
                         ... x 18 lanes            ... x 18 lanes       ... x 18 lanes     
                                │                         │                                
                                │           3_hto_demux_by_batch.R                         
      - demultiplexing          │              ┌──────────┴─────┐                          
      - doublet exclusion       │              │    HTODemux    │                          
      - cell qc                 │              └──────────┬─────┘             Saved:       
                                │    3_Multiseq_Demultiplex_by_batch.R   background drops  
                                ▼              ┌──────────┴─────┐                          
                       ┌────────────────┐      │    multiseq    │         empty droplets   
                       │    Demuxlet    │      │  deMULTIplex   │          ┌────────┐      
                       └────────────────┘      └────────────────┘          │        │      
                                │                       │              ADT │        │      
                                └────────────┬──────────┘                  │        │      
                                             ▼                             └────────┘      
                       Merge quality-controlled singlets                        │          
                       4_validate_multiseq_w_demuxlet_singlets.R                │          
                       4_validate_htodmx_w_demuxlet_singlets.R                  │          
                       5_save_merged_Seurat_singlet_object.R                    │          
                                                                                │          
                  ┌───────────────────────────────────────────────────┐         │          
                  │ h1_h5_merged_seurat_object_demultiplexed_sng.rds  │         │          
                  └───────────────────────────────────────────────────┘         │          
                                                                                │          
                                 cells             cells            ┌───────────┘          
                               ┌────────┐        ┌────────┐         │                      
                           gene│        │        │        │         │                      
                               │        │    ADT │        │         │                      
                               │        │        │        │         │                      
                               └────────┘        └────────┘         │                      
                                    │                 │             │                      
                                    │                 └─────────────┤                      
                                    │                               │                      
                                    ▼                               ▼                      
      - normalization      ┌────────────────┐      ┌─────────────────────────────────┐     
      - clustering         │ scran RNA norm │      │ dsb normalization and denoising │     
                           └────────────────┘      └─────────────────────────────────┘     
                                                                    │                      
                                1_dsbnorm_prot_scrannorm_rna.R      │                      
                                                                    ▼                      
                                               ┌────────────────────────────────────────┐  
                                               │        protein based clustering        │  
                                               │          cell type annotation          │  
                                               └────────────────────────────────────────┘  
                                                                    │                      
                                               1_merged_h1h5_adt_clustering.R              
                                               2_joint_cluster_annotation.R                
                                                                    │                      
                                       ┌────────────────────────────┘                      
                                       │                                                   
                                       │                                                   
                                       │                                                   
                                       ▼                                                   
                   ┌───────────────────────────────────────┐                               
                   │h1h5_annotated_with_meta.rds           │                               
                   │    mRNA                               │                               
                   │        - raw counts                   │                               
                   │        - scran normalized values      │                               
                   │    ADT                                │                               
                   │        - raw counts                   │                               
                   │        - dsb normalized values        │                               
                   │    cell meta data                     │                               
                   │        - donor information            │                               
                   └───────────────────────────────────────┘                               

**Software package versions**

*R version 4.0.5*

    dsb_1.0.2                  
    ggsci_2.9                
    sp_1.4-5              
    SeuratObject_4.1.0       
    Seurat_4.0.1            
    viridis_0.5.1            
    viridisLite_0.3.0        
    scglmmr_0.1.0            
    variancePartition_1.25.6 
    BiocParallel_1.24.1     
    limma_3.46.0             
    magrittr_2.0.3           
    lme4_1.1-26              
    Matrix_1.4-1             
    emmeans_1.5.4           
    Rcpp_1.0.9               
    HDStIM_0.1.0             
    here_1.0.1               
    forcats_0.5.1            
    stringr_1.4.0           
    dplyr_1.0.4              
    purrr_0.3.4              
    readr_1.4.0              
    tidyr_1.1.2              
    tibble_3.1.8            
    ggplot2_3.3.3            
    tidyverse_1.3.0         

*R version 3.5.3*

    scglmmr_0.1.0       
    ggraph_1.0.2        
    igraph_1.2.4.1      
    ggsci_2.9           
    ggridges_0.5.1      
    monocle_2.10.1     
    DDRTree_0.1.5       
    irlba_2.3.3         
    VGAM_1.1-1          
    Biobase_2.42.0      
    BiocGenerics_0.28.0 
    viridis_0.5.1      
    viridisLite_0.3.0   
    here_0.1            
    Seurat_2.3.4        
    Matrix_1.2-15       
    cowplot_0.9.4       
    magrittr_2.0.1     
    forcats_0.4.0       
    stringr_1.4.0       
    dplyr_0.8.5         
    purrr_0.3.3         
    readr_1.3.1         
    tidyr_1.0.2        
    tibble_2.1.1        
    ggplot2_3.1.1       
    tidyverse_1.2.1 
