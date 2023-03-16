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

  