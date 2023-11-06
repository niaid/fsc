## rev
# R version 4.0.5 -> now switching to R --
library(Seurat)
suppressMessages(library(here))
suppressMessages(library(tidyverse))
source(here('functions/MattPMutils.r'))
source(here('functions/scglmmr_functions/pseudobulk_helpers.r'))
source(here('functions/scglmmr_functions/enrichment_analysis.r'))
source(here('functions/scglmmr_functions/model_result_interaction.r'))
source(here('functions/scglmmr_functions/single_cell_gene_module_scores.r'))

# save paths
figpath = here('mid_res/rev/rev.figs/');dir.create(figpath)

##################################
# R 1 
##################################


# Sequencing saturation 
d = read.csv(file = 'data/full_metadata/full_sample_metadata.txt',sep = '\t')
# saturation 
sat = read.csv(file = here('mid_res/rev/rev.data/sequencing.saturation.txt'), header = TRUE, sep = '\t')

p = ggplot(sat, aes(x = lane, y = sequencing.saturation.cDNA)) + 
  geom_bar(fill = 'deepskyblue3', stat = 'identity') + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 11), limits = c(0,100)) + 
  ylab('CITE-seq mRNA sequencing saturation %') + 
  xlab('10x Chromium lane') 
ggsave(p, filename = paste0(figpath, 'read depth.pdf'), width = 4, height = 3)


# UMI per cell variaiton explained 
# read single cell data object 
s = readRDS(file = here('data/h1h5_annotated_with_meta.rds'))

# use models from fig 3 as example 
# scored 
md = s@meta.data %>% 
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
norm.rna = s@data[ ,rownames(md)]

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
names(mod_scores) = names(ct.md)


dd2 = data.frame(module.score = c(
  mod_scores$BC_Naive$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod_scores$CD14_Mono$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod_scores$CD16_Mono$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod_scores$CD4_CD25_Tcell$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod_scores$CD4_Efct_Mem_Tcell$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod_scores$CD8_CD161_Tcell$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod_scores$CD8_Mem_Tcell$`LI.M7.2 enriched in NK cells (I)`, 
  mod_scores$MAIT_Like$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod_scores$mDC$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod_scores$NK$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING),
  rbind(ct.md$BC_Naive,ct.md$CD14_Mono, ct.md$CD16_Mono, ct.md$CD4_CD25_Tcell, ct.md$CD4_Efct_Mem_Tcell, 
        ct.md$CD8_CD161_Tcell, ct.md$CD8_Mem_Tcell, ct.md$MAIT_Like, ct.md$mDC, ct.md$NK)
)
                 
                 
p = ggplot(dd2, aes(x = nUMI, y = module.score )) + 
  theme_bw() + 
  geom_smooth(method = 'lm', color = 'purple') +
  facet_wrap(~celltype_joint,nrow = 2)+
  ggpubr::stat_cor(method = "pearson", aes(label = ..r.label..)) + 
  geom_point(data = dd2, mapping = aes(x = nUMI, y = module.score), size = 0.4, alpha = 0.1) + 
  ylab('single cell module score') +
  ggtitle('NORMALIZED data (as used in manuscript) - variance explained by nUMI - single cells') 
p
ggsave(p, filename = paste0(figpath, 'norm.vpartumi.singlecell.png'),width = 9 , height = 4)

# non-normalized 
# subset normalized RNA data 
raw.rna = s@raw.data[ ,rownames(md)]

# module score for each cell type of specific baseline enriched leading edge genes.
mod.raw = list()
for (i in 1:length(ct.md)) {
  
  # init data for subset i 
  rna = raw.rna[ ,rownames(ct.md[[i]])]
  mod.list = li.g0[[i]]
  
  # calculate single cell score for baseline-enriched module 
  mod.raw[[i]] = WeightedCellModuleScore(
    gene_matrix = rna, 
    module_list = mod.list, 
    threshold = 0, 
    cellwise_scaling = TRUE, 
    return_weighted = FALSE 
  )
  
  # add a "null" score of Gaussian noise as a reference 
  mod.raw[[i]]$null = rnorm(n = nrow(mod_scores[[i]]), mean = 0, sd = 1)
}
names(mod.raw) = names(ct.md)

dd3 = data.frame(module.score = c(
  mod.raw$BC_Naive$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod.raw$CD14_Mono$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod.raw$CD16_Mono$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod.raw$CD4_CD25_Tcell$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod.raw$CD4_Efct_Mem_Tcell$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod.raw$CD8_CD161_Tcell$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod.raw$CD8_Mem_Tcell$`LI.M7.2 enriched in NK cells (I)`, 
  mod.raw$MAIT_Like$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod.raw$mDC$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
  mod.raw$NK$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING),
  rbind(ct.md$BC_Naive,ct.md$CD14_Mono, ct.md$CD16_Mono, ct.md$CD4_CD25_Tcell, ct.md$CD4_Efct_Mem_Tcell, 
        ct.md$CD8_CD161_Tcell, ct.md$CD8_Mem_Tcell, ct.md$MAIT_Like, ct.md$mDC, ct.md$NK)
)
  

p = ggplot(dd3, aes(x = nUMI, y = module.score )) + 
  theme_bw() + 
  geom_smooth(method = 'lm', color = 'red') +
  facet_wrap(~celltype_joint, nrow = 2)+
  ggpubr::stat_cor(method = "pearson", aes(label = ..r.label..)) + 
  geom_point(data = dd3, mapping = aes(x = nUMI, y = module.score), size = 0.4, alpha = 0.1) + 
  ylab('single cell module score') +
  ggtitle('NON NORMALIZED data - variance explained by nUMI - single cells') 
p
ggsave(p, filename = paste0(figpath, 'non.norm.vpartumi.singlecell.png'), width = 9 , height = 4)

# Mixed effects model effect size variance with and without covariates.
# data from original workflow 
# mm res 
mm1 = readRDS(file = here('../Flu_CITEseq_analysis/mid_res/baseline_response/dataV3/mm1.rds'))
mm2 = readRDS(file = here('../Flu_CITEseq_analysis/mid_res/baseline_response/dataV3/mm2.rds'))

dd = data.frame(celltype = mm1$celltype, module = mm1$module, effect1 = mm1$estimatetime0_group2vs1, effect2 = mm2$estimatetime0_group2vs1)
p = ggplot(dd %>%  filter(abs(effect1) > 0.1), aes(x = effect2, y = effect1) ) +
  theme_bw() + 
  facet_wrap(~celltype, nrow = 2) + 
  geom_point(color = 'deepskyblue3')  + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0)  +
  geom_abline(linetype = 'dashed') + 
  xlab('model 2 (multivariate) effect size: expression ~ response + cells per donor + age + sex') + 
  ylab('model 1 (unadjusted) effect size: expression ~ response')
ggsave(p, filename = paste0(figpath, 'effect.sizes.singlecell.module.score.pdf'), width = 9 , height = 4)


# calc n cells per individual per cell type 
s@meta.data %>% 
  group_by(sample, celltype_joint) %>% 
  tally() %>%  
  ungroup() %>% 
  group_by(celltype_joint) %>% 
  summarize_at(.vars = 'n', .funs = c('median' = median, 'sd' = sd)) %>% 
  print(n=50)
# output
# celltype_joint      median     sd
# <chr>                <dbl>  <dbl>
#   1 BC_Mem                72.5  54.4 
# 2 BC_Naive             206.  193.  
# 3 CD103_Tcell           12     9.83
# 4 CD14_Mono            448.  331.  
# 5 CD16_Mono             71.5  36.5 
# 6 CD38_Bcell             3     6.74
# 7 CD4Naive_Tcell       539   256.  
# 8 CD4_CD161_Mem_Tcell  124.   57.7 
# 9 CD4_CD25_Tcell        65.5  33.4 
# 10 CD4_CD56_Tcell         1    12.1 
# 11 CD4_CD57_Tcell        18    96.3 
# 12 CD4_Efct_Mem_Tcell   249   112.  
# 13 CD8_CD161_Tcell       77    59.0 
# 14 CD8_Mem_Tcell         84   102.  
# 15 CD8_NKT                3.5  60.0 
# 16 CD8_Naive_Tcell      294.  158.  
# 17 DOUBLET               22    12.8 
# 18 HSC                    5.5   3.35
# 19 IgA_CD14_Mono          8    10.7 
# 20 MAIT_Like             39.5  77.3 
# 21 NK                   160.  125.  
# 22 mDC                   29.5  20.7 
# 23 pDC                   14    10.0 

# make fig 
cdat=  s@meta.data %>% 
  group_by(sample, celltype_joint) %>% 
  tally() %>%  
  ungroup() %>% 
  group_by(celltype_joint) 
csub = cdat %>%  filter(celltype_joint %in% c('mDC', 'MAIT_Like', 'CD4_CD25_Tcell'))

p = ggplot(csub, aes(x = celltype_joint, y = n )) + 
  theme_bw() +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape = 21, fill = 'deepskyblue3', alpha = 0.7) + 
  ggtitle('lower frequency cell types') + 
  ylab('number of cells (per sample - full cohort)') +
  xlab('')
ggsave(p, filename = paste0(figpath, 'ncell.low.freq.pdf'), width = 4, height = 3)

# n cells per mDC and plasmablast 
s@meta.data %>% 
  group_by(sample) %>% 
  filter(celltype_joint == 'mDC') %>% 
  tally() %$% n %>% 
  quantile()
# 0%  25%  50%  75% 100% 
# 2.0 19.5 29.5 43.5 80.0 

s@meta.data %>% 
  group_by(sample, celltype_joint) %>% 
  summarize_at(.funs = 'tally', .vars = 'celltype_joint') %>% 
  tally %$% n %>% 
  quantile()
quantile
filter(celltype_joint == 'mDC') %>% 
  tally() %$% n %>% 
  quantile()

s@meta.data %>% 
  group_by(sample) %>% 
  filter(celltype_joint == 'CD38_Bcell') %>% 
  tally() %$% n %>% 
  quantile()
# 0%  25%  50%  75% 100% 
# 1    2    3    5   34 

 











 