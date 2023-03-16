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



