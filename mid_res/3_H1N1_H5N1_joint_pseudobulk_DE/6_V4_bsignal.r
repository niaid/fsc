# b cell figures 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
source('functions/scglmmr.functions.R')
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
  filter(padj < 0.05) 

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
ggsave(p,filename = paste0(figpath, 'BCNaive.as03.enrichment.pdf'), width = 6, height = 3)


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
mods = c("CD40_ACT", 
         "REACTOME_ACTIVATION_OF_BH3_ONLY_PROTEINS", 
        "LI.M160 leukocyte differentiation", 
        "REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS") 
cd40 = readRDS('signature_curation/combined_sig_sub.rds')['CD40_ACT']

# Define apoptosis signature
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
sig.test=readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/bsig_data/sig.test.rds'))

data.table::fwrite(list(sig.test$apoptosis.signature),file = paste0(datapath, 'bsig.apoptosis.txt'),sep = '\t')
data.table::fwrite(list(sig.test$CD40_ACT),file = paste0(datapath, 'bsig.cd40.txt'),sep = '\t')

# gsea 
# load contrast fit results 
fit12e = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/fit12e.rds'))
toprank = ExtractResult(
  model.fit.list = fit12e,
  what = 'lmer.z.ranks',
  coefficient.number = 1,
  coef.name = 'delta'
)

gs.bsig = fgsea::fgsea(pathways = list('apoptosis.signature' = sig.test$apoptosis.signature), stats = toprank$BC_Naive)
gs.bsig$leadingEdge

p = fgsea::plotEnrichment(pathway = sig.test$apoptosis.signature,stats = toprank$BC_Naive) + 
  geom_line(size = 2, color = 'deepskyblue3')
ggsave(p, filename = paste0(figpath, 'apoptosis.sig.gsea.citeseq.pdf'), width = 4, height = 3)


##################
# fit single cell model
# score modules 
ms = WeightedCellModuleScore(gene_matrix = s@assays$RNA@data, 
                             module_list = sig.test, 
                             cellwise_scaling = FALSE,
                             return_weighted = FALSE)
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
p
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

