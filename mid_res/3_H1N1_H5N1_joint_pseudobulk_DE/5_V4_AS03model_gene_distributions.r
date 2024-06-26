suppressMessages(library(here))
suppressMessages(library(tidyverse))
#suppressMessages(library(scglmmr))
source('functions/scglmmr.functions.R')
source(here('functions/MattPMutils.r'))

# save path 
datapath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/as03fig/")
figpath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/figuresV4/as03fig/")

# load aggregated tidy data 
av_tidy = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gene_dist/av_tidy.rds'))

# load mixed model fit res v4
fit12e = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/fit12e.rds'))
fitres = ExtractResult(model.fit.list = fit12e, coefficient.number = 1,coef.name = 'delta')

# load gsea results 
gc = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/gc.rds'))

# load leading edge indexed 
li = readRDS(file = 'mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/li.rds')

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


for (i in 1:length(li$CD14_Mono)) {
  topr = fitres$CD14_Mono %>%  filter(P.Value < 0.05)
  names(li$CD14_Mono)[i] %>%  print()
  dplyr::intersect(topr$gene, unlist(li$CD14_Mono[i],use.names = FALSE) ) %>% print()
}

# Format monocyte subset plot 
#gene_highlight = c('MB21D1', 'FPR2', 'P2RY13', 'TLR4')
gene_highlight = c('MB21D1', 'FPR2', 'P2RY13', 'IFIT3')
fitres$CD14_Mono %>% 
  filter(gene %in% gene_highlight) %>% 
  select(gene, P.Value, adj.P.Val, z.std, contrast,  celltype)
# gene      P.Value    adj.P.Val    z.std contrast  celltype
# 1   FPR2 4.950912e-08 0.0004864845 5.453064    delta CD14_Mono
# 2 P2RY13 8.134318e-07 0.0012724398 4.932117    delta CD14_Mono
# 3 MB21D1 4.226564e-04 0.0351672126 3.525517    delta CD14_Mono
# 4  IFIT3 8.919633e-04 0.0520128594 3.322557    delta CD14_Mono
mplt = av_tidy$CD14_Mono %>% filter(gene %in% gene_highlight)
mplt$gene[mplt$gene == "MB21D1"] = "CGAS"
mplt$gene = factor(mplt$gene, levels = c("CGAS",'FPR2', 'P2RY13', 'IFIT3'))

# monocyte subset 
p = ggplot(mplt, aes(x =group, y = count, fill = group , color = group)) +
  mtheme +
  box_gg + 
  facet_wrap(~gene, scales = "free", nrow = 1) + 
  scale_color_manual(values = cu) +
  scale_fill_manual(values = cu.alpha) 
p
ggsave(p, filename = paste0(figpath, 'mono.subset.pdf'), width = 3.5, height = 2.3)



# mDC 
for (i in 1:length(li$mDC)) {
  topr = fitres$mDC %>%  filter(P.Value < 0.05)
  names(li$mDC)[i] %>%  print()
  dplyr::intersect(topr$gene, unlist(li$mDC[i],use.names = FALSE) ) %>% print()
  
}

#gene_highlight = c('FPR1', 'CCR1', 'P2RY13', 'TLR4')
gene_highlight = c('FPR1', 'CCR1', 'P2RY13', 'SCIMP')
fitres$mDC %>% 
  filter(gene %in% gene_highlight)%>% 
  select(gene, P.Value, adj.P.Val, z.std, contrast,  celltype)
# gene      P.Value adj.P.Val    z.std contrast celltype
# 1  SCIMP 0.0001868746 0.1750548 3.736129    delta      mDC
# 2 P2RY13 0.0011814319 0.3280209 3.243325    delta      mDC
# 3   FPR1 0.0133877926 0.5692877 2.473283    delta      mDC
# 4   CCR1 0.0160817638 0.5787543 2.407055    delta      mDC

mplt = av_tidy$mDC %>% filter(gene %in% gene_highlight)
mplt$gene = factor(mplt$gene , levels =c('FPR1', 'CCR1', 'P2RY13', 'SCIMP'))
# mDC subset 
p = ggplot(mplt, aes(x =group, y = count, fill = group , color = group)) +
  mtheme +
  box_gg + 
  facet_wrap(~gene, scales = "free", nrow = 1) + 
  scale_color_manual(values = cu) +
  scale_fill_manual(values = cu.alpha) 
p
ggsave(p, filename = paste0(figpath, 'mdc.subset.pdf'), width = 3., height = 2.3)


# B cell Naive
for (i in 1:length(li$BC_Naive)) {
  names(li$BC_Naive)[i] %>%  print()
  dplyr::intersect(topgene$BC_Naive, unlist(li$BC_Naive[i],use.names = FALSE) ) %>%
    print()
}
#bplt = av_tidy$BC_Naive %>% filter(gene == 'PMAIP1')
gene_highlight = c('PMAIP1', 'BTG1')
fitres$BC_Naive %>% 
  filter(gene %in% gene_highlight) %>% 
  select(gene, P.Value, adj.P.Val, z.std, contrast,  celltype)
bplt = av_tidy$BC_Naive %>% filter(gene %in%  c('PMAIP1', 'BTG1')) 
bplt$gene[bplt$gene == "PMAIP1"] = "NOXA (PMAIP1)"

# B cell subset 
p = ggplot(bplt, aes(x =group, y = count, fill = group , color = group)) +
  mtheme +
  box_gg + 
  facet_wrap(~gene, scales = "free", nrow = 1) + 
  scale_color_manual(values = cu) +
  scale_fill_manual(values = cu.alpha) 
p
ggsave(p, filename = paste0(figpath, 'Bcell.PMAIP1.pdf'), width = 2.5, height = 2.5)
