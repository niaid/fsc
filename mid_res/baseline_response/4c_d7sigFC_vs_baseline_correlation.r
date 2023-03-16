# correlate module expression avz with d7 FC in antibody predictive signature. 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))
source(here('functions/MattPMutils.r'))
# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/d7cor/");
dir.create(figpath, recursive = TRUE)

# day 7 response signature fc 
d7res = readRDS(file = here('mid_res/baseline_response/dataV3/d7res.rds'))

# baseline expression correlation 
ds = readRDS(file = here('mid_res/baseline_response/dataV3/ds.rds'))

# created shorter names 
new.names = data.table::fread(
  file = here('mid_res/baseline_response/dataV3/baseline.module.name.shortened.txt'),
  sep = '\t'
)
new.names$cname2 = paste(new.names$celltype, new.names$shortname, sep = ' :: ')
colnames(ds) = new.names$cname2

# format 
d7form = ds %>%
  as.data.frame() %>% 
  rownames_to_column('subject') %>% 
  mutate(subject = str_sub(subject, 1,3)) %>% 
  filter(subject %in% colnames(d7res)) %>% 
  column_to_rownames('subject')
saveRDS(d7form, file = paste0(datapath,'d7form.rds'))
d7form = readRDS(file = here('mid_res/baseline_response/dataV3/d7form.rds'))

# pairwise correlation with d7 response 
dd = cbind(t(d7res), as.data.frame(d7form))
saveRDS(dd, file = paste0(datapath,'dd.rds'))
d7.cor = Hmisc::rcorr(as.matrix(dd),type = 'spearman')
saveRDS(d7.cor, file = paste0(datapath,'d7.cor.rds'))
d7.cor = readRDS(file = here('mid_res/baseline_response/dataV3/d7.cor.rds'))

# aes set 
plotattr = list(
  theme_bw(),
  geom_point(shape = 21, size = 3.5, stroke = 0.8), 
  #geom_text(nudge_y = 0.05, size = 3), 
  ylab('Day 7 FC Antibody response signature'),
  theme(axis.title.y = element_text(size = 11)), 
  theme(axis.title.x = element_text(size = 8)), 
  scale_fill_manual(values = c(col.alpha("red",0.7), col.alpha("dodgerblue", 0.7))),
  theme(legend.position = 'none', legend.key.size = unit(0.29, units = 'cm')), 
  theme(aspect.ratio = 1) 
)

# specify response in vis.  
high.responders = c("205","207","209","212","215","234","237","245","250","256")

# calculate and vis correlation between day7 response signature fold change (bulk)
# versus log cpm of the day 7 signatures associated with adjMFC group 
for (i in 1:length(colnames(d7form))) {
  #i = 1 
  mod.names = c(colnames(d7form)[i], colnames(t(d7res)))
  cplot = cbind(as.data.frame(d7form[, i]), t(d7res))
  cplot$subject = rownames(d7form)
  cplot$response = ifelse(cplot$subject %in% high.responders, 'high', 'low')
  colnames(cplot)[1:2] = mod.names
  
  # -- for correlation -- 
  dsub = as.data.frame(cbind(v1 = cplot[ ,1], v2 = cplot[ ,2]))
  
  # calculate and vis. correlation across all subjects; color by response. 
  p = ggpubr::ggscatter(cplot, x = mod.names[1], y = mod.names[2],
                        color = col.alpha('white','0.01'),
                        add.params = list(color = "black", fill = "grey")
                        )  +
    plotattr +
    aes(fill = response) + 
    ggpubr::stat_cor(data = dsub, aes(x = v1, y = v2), method = 'spearman',
                     inherit.aes = FALSE, label.x.npc = "left",label.y.npc = "top", cor.coef.name =  "rho",) 
  # save name by cell type first 
  module = sub("^[^::]*::", "", mod.names[1])
  celltype = gsub("::.*", "", mod.names[1])
  ggsave(p, filename = paste0(figpath, celltype, module, 'd7cor.pdf'), width = 3.1, height = 3.1)
}
  


