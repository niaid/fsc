suppressMessages(library(tidyverse))
suppressMessages(library(here))
source(here('functions/scglmmr.functions.R'))
source(here('functions/MattPMutils.r'))
# output directories 
figpath = here("mid_res/combined_contrast/figures/")
datapath = here('mid_res/combined_contrast/generated_data/')
dir.create(figpath); dir.create(datapath)

# load CITE results 
gc = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/gc.rds'))

#######
# CD14 monocytes cite-seq 
######
mo = gc$CD14_Mono %>%
  as.data.frame() %>% 
  filter(padj < 0.05 & NES > 0)

# shorten names 
mo$pathway[mo$pathway == "LI.M4.3 myeloid cell enriched receptors and transporters" ] <- 'LI.M4.3 myeloid receptors and transporters'
mo$pathway[mo$pathway == "REACTOME_INTERFERON_SIGNALING" ] <- 'reactome interferon signaling'
mo$pathway[mo$pathway == "LI.M37.0 immune activation - generic cluster" ] <- 'LI.M37.0 immune activation'
mo$pathway[mo$pathway == "SLE_SIG" ] <- 'IFN Sig (SLE)'
mo$pathway[mo$pathway == "IFN1_DCACT" ] <- 'IFN I DCACT'

# save adjuvant signatures without IFN sigs 
mo.noifn = mo %>% filter(! pathway %in% c('IFN I DCACT', 'IFN Sig (SLE)', 'reactome interferon signaling', "LI.M127 type I interferon response"))
saveRDS(mo.noifn,file = paste0(datapath, 'mo.noifn.rds'))
mo.noifn$pathway
# [1] "LI.M11.0 enriched in monocytes (II)"        "LI.M16 TLR and inflammatory signaling"     
# [3] "LI.S4 Monocyte surface signature"           "LI.M4.0 cell cycle and transcription"      
# [5] "LI.M37.0 immune activation"                 "LI.M4.3 myeloid receptors and transporters"
# [7] "LI.M118.1 enriched in monocytes (surface)"  "LI.M169 mitosis (TF motif CCAATNNSNNNGCG)" 
# [9] "LI.M67 activated dendritic cells"           "LI.M194 TBA"                               
# [11] "LI.M118.0 enriched in monocytes (IV)"      

###### 
# mDCs cite-seq
######
dc = gc$mDC %>%
  as.data.frame() %>% 
  filter(padj<0.05 & NES > 0)

# shorten names 
dc$pathway[dc$pathway == "REACTOME_MITOTIC_M_M_G1_PHASES" ] <- 'reactome mitosis M G1 phase'
dc$pathway[dc$pathway == "REACTOME_CELL_CYCLE_CHECKPOINTS"  ] <- 'reactome cell cycle checkpoints'
dc$pathway[dc$pathway == "REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_CDC20_AND_OTHER_APC_C_CDH1_TARGETED_PROTEINS_IN_LATE_MITOSIS_EARLY_G1" ] <- 'reactome CDC20 APC C degradation late mitosis early G1'
dc$pathway[dc$pathway == "LI.M4.3 myeloid cell enriched receptors and transporters" ] <- 'LI.M4.3 myeloid receptors and transporters'
dc$pathway[dc$pathway == "REACTOME_VIF_MEDIATED_DEGRADATION_OF_APOBEC3G" ] <- 'reactome VIF APOBEC3G'
dc$pathway[dc$pathway == "REACTOME_DNA_REPLICATION" ] <- 'reactome DNA replication'
dc$pathway[dc$pathway == "REACTOME_CROSS_PRESENTATION_OF_SOLUBLE_EXOGENOUS_ANTIGENS_ENDOSOMES" ] <- 'reactome cross presentation of exogenous antigen'
dc$pathway[dc$pathway == "REACTOME_AUTODEGRADATION_OF_CDH1_BY_CDH1_APC_C" ] <- 'reactome autodegradation of CDH1 by APC C'
dc$pathway[dc$pathway == "REACTOME_MITOTIC_M_M_G1_PHASES" ] <- 'reactome mitosis M G1 phase'
dc$pathway[dc$pathway == "REACTOME_APC_C_CDC20_MEDIATED_DEGRADATION_OF_MITOTIC_PROTEINS" ] <- 'reactome APC C CDC20 degradation of mitotic proteins'
dc$pathway[dc$pathway == "REACTOME_M_G1_TRANSITION" ] <- 'reactome M G1 transition'
dc$pathway[dc$pathway == "REACTOME_EXTRINSIC_PATHWAY_FOR_APOPTOSIS" ] <- 'reactome apoptotic extrinsic pathway'



# IFN these pathways are not present in the dc enrichments. 
dc.noifn = dc 
saveRDS(dc.noifn,file = paste0(datapath, 'dc.noifn.rds'))
dc.noifn$pathway
# [1] "LI.M4.0 cell cycle and transcription"                 "LI.M11.0 enriched in monocytes (II)"                 
# [3] "reactome mitosis M G1 phase"                          "reactome cell cycle checkpoints"                     
# [5] "reactome mitosis late mitosis early G1"               "LI.M4.3 myeloid receptors and transporters"          
# [7] "reactome VIF APOBEC3G"                                "reactome DNA replication"                            
# [9] "LI.M194 TBA"                                          "reactome cross presentation of exogenous antigen"    
# [11] "reactome autodegradation of CDH1 by APC C"            "reactome APC C CDC20 degradation of mitotic proteins"
# [13] "LI.M177.0 TBA"                                        "reactome M G1 transition"                            
# [15] "reactome apoptotic extrinsic pathway"      

#####
# B cells cite-seq 
#####
bn = gc$BC_Naive %>%
  as.data.frame() %>% 
  filter(padj < 0.05) # do not apply NES filter

# add string for cohort 
mo$cohort = 'CITE-seq'
dc$cohort = 'CITE-seq'
bn$cohort = 'CITE-seq'

##################################
# load validation cohort data 
##################################
mv = readRDS(file = here('mid_res/vand/generated_data/mv.rds')) 
dcv = readRDS(file = here('mid_res/vand/generated_data/dcv.rds')) 
bcv = readRDS(file = here('mid_res/vand/generated_data/bcv.rds')) 
mv = as.data.frame(mv$MNC)
dcv = as.data.frame(dcv$DNC)
bcv = as.data.frame(bcv$BCL)

# shorten names 
mv$pathway[mv$pathway == "LI.M4.3 myeloid cell enriched receptors and transporters" ] <- 'LI.M4.3 myeloid receptors and transporters'
mv$pathway[mv$pathway == "LI.M37.0 immune activation - generic cluster" ] <- 'LI.M37.0 immune activation'
mv$pathway[mv$pathway == "SLE_SIG" ] <- 'IFN Sig (SLE)'
mv$pathway[mv$pathway == "IFN1_DCACT" ] <- 'IFN I DCACT'
dcv$pathway[dcv$pathway == "REACTOME_MITOTIC_M_M_G1_PHASES" ] <- 'reactome mitosis M G1 phase'
dcv$pathway[dcv$pathway == "REACTOME_CELL_CYCLE_CHECKPOINTS"  ] <- 'reactome cell cycle checkpoints'
dcv$pathway[dcv$pathway == "REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_CDC20_AND_OTHER_APC_C_CDH1_TARGETED_PROTEINS_IN_LATE_MITOSIS_EARLY_G1" ] <- 'reactome CDC20 APC C degradation late mitosis early G1'
dcv$pathway[dcv$pathway == "LI.M4.3 myeloid cell enriched receptors and transporters" ] <- 'LI.M4.3 myeloid receptors and transporters'
dcv$pathway[dcv$pathway == "REACTOME_VIF_MEDIATED_DEGRADATION_OF_APOBEC3G" ] <- 'reactome VIF APOBEC3G'
dcv$pathway[dcv$pathway == "REACTOME_DNA_REPLICATION" ] <- 'reactome DNA replication'
dcv$pathway[dcv$pathway == "REACTOME_CROSS_PRESENTATION_OF_SOLUBLE_EXOGENOUS_ANTIGENS_ENDOSOMES" ] <- 'reactome cross presentation of exogenous antigen'
dcv$pathway[dcv$pathway == "REACTOME_AUTODEGRADATION_OF_CDH1_BY_CDH1_APC_C" ] <- 'reactome autodegradation of CDH1 by APC C'
dcv$pathway[dcv$pathway == "REACTOME_MITOTIC_M_M_G1_PHASES" ] <- 'reactome mitosis M G1 phase'
dcv$pathway[dcv$pathway == "REACTOME_APC_C_CDC20_MEDIATED_DEGRADATION_OF_MITOTIC_PROTEINS" ] <- 'reactome APC C CDC20 degradation of mitotic proteins'
dcv$pathway[dcv$pathway == "REACTOME_M_G1_TRANSITION" ] <- 'reactome M G1 transition'
dcv$pathway[dcv$pathway == "REACTOME_EXTRINSIC_PATHWAY_FOR_APOPTOSIS" ] <- 'reactome apoptotic extrinsic pathway'

# append with cohort 
dcv$cohort = 'validation'
dcv$celltype = 'sorted DC'

bcv$cohort = 'validation'
bcv$celltype = 'sorted B cells'

mv$cohort = 'validation'
mv$celltype = 'sorted monocytes'

# pathways in CITE hyp set. 
dcv = dcv %>% filter(pathway %in% dc$pathway)
mv = mv %>% filter(pathway %in% mo$pathway)


# combine
col.keep = c('pathway', 'pval', 'padj', 'NES', 'celltype', 'cohort') 
r.list = list(dcv, bcv, mv, mo, dc, bn)
r.list = lapply(r.list, function(x) x %>% select(all_of(col.keep)))
d = bind_rows(r.list)

# group
d$main = ifelse(d$celltype %in% c('mDC', 'sorted DC'), yes = 'DC', no = d$celltype)
d$main = ifelse(d$celltype %in% c('CD14_Mono', 'sorted monocytes'), yes = 'Mono', no = d$main)
d$main = ifelse(d$celltype %in% c('BC_Naive', 'sorted B cells'), yes = 'BC', no = d$main)

# unnate subset 
d2 = d %>% filter(!celltype %in% c('BC_Naive', 'sorted B cells'))
d2$main = factor(d2$main, levels = c('Mono', 'DC'))


# add asterisk for significant validation 
d2 = d2 %>%  filter(!pathway == 'combined.signature')
d3 = 
  d2 %>% 
  mutate(padj.validation = ifelse(cohort == 'validation', padj, no = Inf)) %>% 
  mutate(padj.citeseq = ifelse(cohort == 'CITE-seq', padj,no = Inf)) %>% 
  mutate(pathway.new = ifelse( padj.validation < 0.01, yes = paste0(' * ', pathway), no = pathway))
d3 %>% filter(cohort == 'validation')
d2$pathway = plyr::mapvalues(d2$pathway,from = d3$pathway,to = d3$pathway.new)

d2$cohort = factor(d2$cohort, levels = c("validation", "CITE-seq"))
p = 
  ggplot(d2, aes(x = NES, y = reorder(pathway, NES), size = -log10(padj), fill=cohort)) +
  theme_bw() +
  geom_jitter(stroke = 0.3, shape = 21, height = 0.1, width = 0) + # visualize directly overlapping points 
  scale_fill_manual(values = c(col.alpha('#90C983',0.7), col.alpha('deepskyblue3', 0.5))) + 
  ylab("") +
  xlab('Normalized Enrichment Score') + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  facet_grid(vars(main), scales = 'free', space = 'free') +
  theme_bw(base_size = 9) + 
  theme(axis.text = element_text(color = 'black')) +
  guides(fill = guide_legend(override.aes = list(size = 5))) 
p
ggsave(p, filename = paste0(figpath,'combined_as03_model_withifn.pdf'), width = 5.5, height = 3.8)
  

# B cells 
d3 = d %>% 
  filter(celltype %in% c('BC_Naive', 'sorted B cells')) %>% 
  filter(pathway %in% c(bn$pathway, 'apoptosis.signature', 'CD40_ACT') )


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



