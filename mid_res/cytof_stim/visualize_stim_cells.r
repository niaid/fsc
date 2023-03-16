suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(HDStIM))
library(Rcpp)
library(emmeans)
library(lme4)
source('functions/MattPMutils.r')
# save paths 
figpath = here('mid_res/stim/figures/')
datapath = here('mid_res/stim/generated_data/')
dir.create(figpath); dir.create(datapath)

# read data from stim cell selector 
d = readRDS(file = here('data/stim/mapped_data.rds'))

d$umap_plot_data  

up = HDStIM::plot_umap(mapped_data = d)
pd = up[[4]]
pd2 = pd$data
pd2$stim = ifelse(str_sub(pd2$response_status, -5,-1) == 'Stim.', yes = 'LPS stimulated', no = 'unstimulated')
cf = ggsci::pal_d3( palette = 'category20', 0.5)(2) %>% rev
p = 
ggplot(pd2, aes(x = UMAP1, y = UMAP2, fill = stim )) + 
  geom_point(shape = 21, stroke = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = cf)  +
  ggtitle(pd$labels$title)
ggsave(p,filename = paste0(figpath, 'umap_lps_mono.png'), width = 5, height = 4)
ggsave(p,filename = paste0(figpath, 'umap_lps_mono_outline.pdf'), width = 5, height = 4)
