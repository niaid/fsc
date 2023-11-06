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
datapath = here('mid_res/rev/rev.data/')


# read single cell data object 
s = readRDS(file = here('data/h1h5_annotated_with_meta.rds'))

##################################
# R 2  
##################################

# n cells cd4 mem tfh calc 
med.cd4mem = s@meta.data %>% 
  group_by(sample, celltype_joint) %>% 
  tally() %>% 
  filter(celltype_joint == 'CD4_Efct_Mem_Tcell') %$% 
  n %>% 
  median() 
est.ctfh = 0.2*0.1*med.cd4mem
est.ctfh
# [1] 4.98




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
my_compare_2 = list( l2 = c("Low.F", "Low.M"))

# define adjuvant signatures 
# leading edge combined signature 
# upregulated gnees only 
# mono and mDC specific and combined genes
mono.validated = readRDS(file = here('mid_res/nat_adj/generated_data/V4/mono.as03.sig.validated.rds'))
dc.validated = readRDS(file = here('mid_res/nat_adj/generated_data/V4/dc.as03.sig.validated.rds'))
adjuvant.signatures = c(list('Mono_AS03_validated' = mono.validated), list('DC_AS03_validated' = dc.validated))

# load average baseline high  and low responders from unadjuvanted cohort 
av0 = readRDS(file = here('mid_res/nat_adj/generated_data/V4/av0.rds'))


#############################
# CD14 Monocytes 
#############################
# Monocyte combined AS03 signature average across time between groups 
#natural adjuvant test
mono.sig.av2 = 
  av0$CD14_Mono %>% 
  filter(gene %in% adjuvant.signatures$Mono_AS03_validated) %>% 
  group_by(sample, response) %>% 
  summarize(meansig = mean(count))

# append sex and age 
meta=read.csv(file = here('data/full_metadata/full_sample_metadata.txt'),sep = '\t')

# append subject id 
mono.sig.av2 = mono.sig.av2 %>% mutate(subjectid = str_sub(sample, 1,3))

# mapo sex and age to average data based on metadata match 
mono.sig.av2$sex = plyr::mapvalues(x = mono.sig.av2$subjectid,from = meta$subjectid,to = meta$gender)
mono.sig.av2$age = plyr::mapvalues(x = mono.sig.av2$subjectid,from = meta$subjectid,to = meta$age)
mono.sig.av2$age = as.numeric(mono.sig.av2$age)
#mono.sig.av2$response.sex = paste(mono.sig.av2$response, mono.sig.av2$sex,sep = '.')

p = ggplot(mono.sig.av2, aes(x = sex, y = meansig, fill = sex , color = sex)) +
  theme_bw() +
  geom_boxplot() +
  ggsci::scale_color_d3() +
  ggsci::scale_fill_d3(alpha = 0.5) +
  ggpubr::stat_compare_means(paired = FALSE, show.legend = FALSE, method = 'wilcox') + 
  geom_jitter() + 
  ylab('Monocyte AS03 Adjuvant Signature') +
  ggtitle('Baseline') + 
  theme(legend.position = 'none') + 
  xlab('sex')
p
ggsave(p, filename = paste0(figpath, 'MONO.baseline.naturaladjuvant.sex.pdf'), width = 2.1, height = 4)


## Repeat for mDC 
mdc.sig.av = 
  av0$mDC %>% 
  filter(gene %in% adjuvant.signatures$DC_AS03_validated) %>% 
  group_by(sample, response) %>% 
  summarize(meansig = mean(count))

# append sex and age 
# meta=read.csv(file = here('data/full_metadata/full_sample_metadata.txt'),sep = '\t') # done above

# append subject id 
mdc.sig.av = mdc.sig.av %>% mutate(subjectid = str_sub(sample, 1,3))

# mapo sex and age to average data based on metadata match 
mdc.sig.av$sex = plyr::mapvalues(x = mdc.sig.av$subjectid,from = meta$subjectid,to = meta$gender)
mdc.sig.av$age = plyr::mapvalues(x = mdc.sig.av$subjectid,from = meta$subjectid,to = meta$age)
mdc.sig.av$age = as.numeric(mdc.sig.av$age)
#mdc.sig.av$response.sex = paste(mdc.sig.av$response, mdc.sig.av$sex,sep = '.')


#plot
p = ggplot(mdc.sig.av, aes(x = sex, y = meansig, fill = sex , color = sex)) +
  theme_bw() +
  geom_boxplot() +
  ggsci::scale_color_d3() +
  ggsci::scale_fill_d3(alpha = 0.5) +
  ggpubr::stat_compare_means(paired = FALSE, show.legend = FALSE, method = 'wilcox') + 
  geom_jitter() + 
  ylab('mDC AS03 Adjuvant Signature') +
  ggtitle('Baseline') + 
  theme(legend.position = 'none') + 
  xlab('sex')
p
ggsave(p, filename = paste0(figpath, 'mDC.baseline.naturaladjuvant.sex.pdf'), width = 2.1, height = 4)


##################
# Age effects 
##################

# age vs nat adj signature 
p = ggplot(mono.sig.av2, aes(x = age, y = meansig )) + 
  theme_bw() +
  geom_point() + 
  geom_smooth(method = 'lm') +
  facet_wrap(~response) + 
  ggpubr::stat_cor(method = 'spearman') + 
  ylab('CD14 Mono Natural adjuvant signature') + 
  ggtitle('CD14 monocytes baseline unadjuvanted subjects')
p
ggsave(p, filename = paste0(figpath, 'mono.natadj.vs.age.pdf'), width = 6, height = 3)

p2 = ggplot(mono.sig.av2, aes(x = age, y = meansig)) + 
  theme_bw() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = 'spearman') + 
  geom_point(aes(x = age, y = meansig, color = response)) +
  scale_color_manual(values = c('dodgerblue', 'red'))+ 
  ylab('CD14 Mono Natural adjuvant signature') +  
  ggtitle('CD14 monocytes ')
p2
p3 = cowplot::plot_grid(p, p2,rel_widths = c(1.3,1))
ggsave(p3, filename = paste0(figpath, 'mono.natadj.vs.age.combined.pdf'), width = 9, height = 3)

## 
p = ggplot(mdc.sig.av, aes(x = age, y = meansig )) + 
  theme_bw() +
  geom_point() + 
  geom_smooth(method = 'lm') +
  facet_wrap(~response) + 
  ggpubr::stat_cor(method = 'spearman') + 
  ylab('mDC Natural adjuvant signature') + 
  ggtitle('mDC monocytes baseline unadjuvanted subjects')
ggsave(p, filename = paste0(figpath, 'mDC.natadj.vs.age.pdf'), width = 6, height = 3)

p2 = ggplot(mdc.sig.av, aes(x = age, y = meansig)) + 
  theme_bw() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = 'spearman') + 
  geom_point(aes(x = age, y = meansig, color = response)) +
  scale_color_manual(values = c('dodgerblue', 'red'))+ 
  ylab('mdc Natural adjuvant signature') +  
  ggtitle('CD14 monocytes ')
p2
p3 = cowplot::plot_grid(p, p2,rel_widths = c(1.3,1))
ggsave(p3, filename = paste0(figpath, 'mdc.natadj.vs.age.combined.pdf'), width = 9, height = 3)



## Variance partition analysis 
# age vs age-related signatures From Figure 1f within CD8 naive CD161 T cells and CD8 naive T cells 
age.pos.cd161 = readRDS(file = here('mid_res/variance_partition/generated_data2/age.pos.cd161.rds'))
age.pos.cd8n = readRDS(file = here('mid_res/variance_partition/generated_data2/age.pos.cd8n.rds'))

# CD161 age associated signature stratified by response 
cd161av = 
  av0$CD8_CD161_Tcell %>% 
  filter(gene %in% age.pos.cd161) %>% 
  group_by(sample, response) %>% 
  summarize(meansig = mean(count))
cd161av = cd161av %>% mutate(subjectid = str_sub(sample, 1,3))
cd161av$age = plyr::mapvalues(x = cd161av$subjectid,from = meta$subjectid,to = meta$age)
cd161av$age = as.numeric(cd161av$age)

# analyze correlation stratified by response
p = ggplot(cd161av, aes(x = age, y = meansig )) + 
  theme_bw() +
  geom_point() + 
  geom_smooth(method = 'lm') +
  facet_wrap(~response) + 
  ggpubr::stat_cor(method = 'spearman') + 
  ylab('CD8+CD161+ T Age-assoc. gene average') + 
  ggtitle('CD8 CD161 T cell baseline unadjuvanted subjects')
ggsave(p, filename = paste0(figpath, 'CD161agesig.vs.response.pdf'), width = 6, height = 3)

cd8nav = 
  av0$CD8_Naive_Tcell %>% 
  filter(gene %in% age.pos.cd8n) %>% 
  group_by(sample, response) %>% 
  summarize(meansig = mean(count))
cd8nav = cd8nav %>% mutate(subjectid = str_sub(sample, 1,3))
cd8nav$age = plyr::mapvalues(x = cd8nav$subjectid,from = meta$subjectid,to = meta$age)
cd8nav$age = as.numeric(cd8nav$age)

# analyze correlation stratified by response
p = ggplot(cd8nav, aes(x = age, y = meansig )) + 
  theme_bw() +
  geom_point() + 
  geom_smooth(method = 'lm') +
  facet_wrap(~response) + 
  ggpubr::stat_cor(method = 'spearman') + 
  ylab('CD8 Naive T Age-assoc. gene average') + 
  ggtitle('CD8 Naive T cell baseline unadjuvanted subjects')
p
ggsave(p, filename = paste0(figpath, 'CD8naive.agesig.vs.response.pdf'), width = 6, height = 3)
# 
