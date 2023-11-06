# Save results table 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
#suppressMessages(library(scglmmr))
source('functions/scglmmr.functions.R')
source(here('functions/MattPMutils.r'))
datapath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/bsig_data/")
figpath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/figuresV4/bsig/")

# read gsea and mixed model results 
fit12e = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/fit12e.rds'))
res = ExtractResult(model.fit.list = fit12e, coefficient.number = 1, coef.name = 'delta')
gc = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/gc.rds'))
apop = c('PMAIP1 (NOXA)', 'BCL2', 'BTG2', 'BTG1')
res$BC_Naive$gene[res$BC_Naive$gene == 'PMAIP1'] <- 'PMAIP1 (NOXA)'
sig.test = readRDS(here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/bsig_data/sig.test.rds'))
apop.genes = sig.test$apoptosis.signature
p = 
  ggplot(res$BC_Naive, aes(x = logFC, y = z.std)) + 
  theme_bw() + 
  xlim(-3,2.1) + 
  geom_bin2d(bins = 400, show.legend = FALSE, fill = 'black') +
  ylab('Mixed model contrast \n standardized z statistic') + 
  xlab('Difference in day 1 log fold change\nAS03 vs unadjuvanted') +
  geom_point(data = res$BC_Naive %>%  filter(gene %in% c(apop , apop.genes)), aes(x = logFC, y = z.std),
             color = 'deepskyblue3', size = 2) +
  ggrepel::geom_text_repel(data = res$BC_Naive %>%  filter(gene %in% apop), 
                           aes(x = logFC, y = z.std, label = gene),
                           color = 'red',size = 5,
                           nudge_x = -0.6, seed = 2,
                           nudge_y = 0, box.padding = 1,
                           segment.size = 0.1) + 
  theme(title = element_text(size = 18)) + 
  geom_vline(xintercept = 0, linetype ='dashed') + 
  ggtitle('Naive B cells')  
p
ggsave(p,filename = paste0(figpath,'bcn_contrast_genes.pdf'), width = 5, height = 5)
