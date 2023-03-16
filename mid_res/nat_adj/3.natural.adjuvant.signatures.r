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
