suppressMessages(library(here))
suppressMessages(library(tidyverse))
source(here('functions/scglmmr.functions.R'))
source(here('functions/MattPMutils.r'))
set.seed(1990)
# set save paths 
figpath = here("mid_res/nat_adj/figures/V4/")
datapath = here("mid_res/nat_adj/generated_data/V4/")


# set theme 
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


# load average day 1 comparison cohort data 
av_tidy = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gene_dist/av_tidy.rds'))

# AS03 adjuvant signatures (no ifn)
as03.sig.list = readRDS(file = here('mid_res/nat_adj/generated_data/V4/as03.sig.list.rds'))

# mdc Combined AS03 signature average across time between groups 
mdc.sig.av = 
  av_tidy$mDC %>% 
  filter(gene %in% as03.sig.list$AS03_mDC) %>% 
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
ggsave(p,filename = paste0(figpath, 'as03_mDC_sig.2.pdf'), width = 1.9, height = 3)

# monocyte Combined AS03 signature average across time between groups 
mono.sig.av = 
  av_tidy$CD14_Mono %>% 
  filter(gene %in% as03.sig.list$AS03_Monocyte) %>% 
  group_by(sample, group) %>% 
  summarize(meansig = mean(count))
#plot
p = ggplot(mono.sig.av, aes(x = group, y = meansig, fill = group , color = group)) +
  mtheme + 
  theme(axis.title.x = element_blank()) +
  ylab('CD14 Mono AS03 Adjuvant Signature') +
  scale_color_manual(values = cu) +
  scale_fill_manual(values = cu.alpha) + 
  ggtitle('CD14 Monocytes')
ggsave(p,filename = paste0(figpath, 'as03_mono_sig.2.pdf'), width = 1.9, height = 3)


##########################
## GSEA of NA signatures 
##########################
# enrichment of validated signatures (Vand cohort) in baseline high vs low responders  
mono.as03.sig.validated = readRDS(file = here('mid_res/nat_adj/generated_data/V4/mono.as03.sig.validated.rds'))
dc.as03.sig.validated = readRDS(file = here('mid_res/nat_adj/generated_data/V4/dc.as03.sig.validated.rds'))

#  high vs low model gene ranks within mono and mDC age and sex adjusted 
cont0 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/cont0.rds'))
r0 = ExtractResult(model.fit.list = cont0, what = 'gene.t.ranks',coefficient.number = 1, coef.name = 'adjmfc')

# enrichment 
mono.na.gsea = fgsea::fgsea(pathways =  list('AS03.mono' = mono.as03.sig.validated), stats = r0$CD14_Mono)
# pathway         pval         padj  log2err        ES     NES size                                    leadingEdge
# 1: AS03.mono 7.515986e-13 7.515986e-13 0.921426 0.6263153 2.74048   78 S100A11,S100A12,APOBEC3A,DDX60,DDX58,CREB5,...

p = fgsea::plotEnrichment(pathway = mono.as03.sig.validated, stats = r0$CD14_Mono) +
  geom_line(size = 1.5, color = 'red') + 
  theme(plot.title = element_text(size = 10))
p
ggsave(p, filename =paste0(figpath, 'mono.natadj.gsea.2.pdf'), width = 5, height = 3)

mono.na.gsea$leadingEdge
# [1] "S100A11"  "S100A12"  "APOBEC3A" "DDX60"    "DDX58"    "CREB5"    "TFEC"     "S100A9"   "GCA"      "FGL2"     "NACC2"   
# [12] "KYNU"     "SLC16A3"  "MS4A4A"   "S100A8"   "FCGR3A"   "SAMHD1"   "TNFSF13B" "C19orf59" "CCR2"     "MARCO"    "P2RY13"  
# [23] "RSAD2"    "SERPINA1" "FPR1"     "FGR"      "PLSCR1"   "SIGLEC9"  "IFIH1"    "ACSL1"    "LMNB1"    "LRRK2"    "MNDA"    
# [34] "PLBD1"    "KIAA0513" "AQP9"     "SLC31A2"  "LILRB1"   "VCAN"    

data.table::fwrite(mono.na.gsea$leadingEdge,file = paste0(datapath,'mono.na.gsea.leadingEdge.txt'))


mdc.na.gsea = fgsea::fgsea(pathways = list('AS03.mdc' = dc.as03.sig.validated), stats = r0$mDC)
# pathway       pval       padj   log2err        ES      NES size                                  leadingEdge
# 1: AS03.mdc 0.02383525 0.02383525 0.3524879 0.3713362 1.521304   58 S100A8,S100A9,PSMB6,SERPINA1,RB1,SLC31A2,...
p = fgsea::plotEnrichment(pathway = dc.as03.sig.validated,stats = r0$mDC)  +
  geom_line(size = 1.5, color = 'red') + 
  theme(plot.title = element_text(size = 10))
p
ggsave(p, filename =paste0(figpath, 'mdc.natadj.gsea.pdf'), width = 5, height = 3)

mdc.na.gsea$leadingEdge
# [1] "S100A8"   "S100A9"   "PSMB6"    "SERPINA1" "RB1"      "SLC31A2"  "TYMP"     "S100A11"  "MS4A4A"   "FCN1"     "LILRB2"  
# [12] "PSMA7"    "CDC26"    "RBX1"     "PSMB3"    "PLBD1"    "LMNB1"    "PPP2R5E"  "KYNU"    

data.table::fwrite(mdc.na.gsea$leadingEdge,file = paste0(datapath,'mdc.na.gsea.leadingEdge.txt'))

