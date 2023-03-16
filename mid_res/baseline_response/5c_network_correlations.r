# baseline nettwork correlation 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(igraph))
source('functions/MattPMutils.r')

# set paths 
datapath = here("mid_res/baseline_response/dataV3")
figpath = here("mid_res/baseline_response/figuresV3/network_correlations/");
dir.create(figpath,recursive = TRUE)


# load of baseline module expression across donors only of 
# leading edge genes from baseline enrichments 
ds = readRDS(here('mid_res/baseline_response/dataV3/ds.rds'))
#data.table::fwrite(ds,file = paste0(here('git_ignore/ds.csv')),sep = ',')
# created shorter names -- read these in 
new.names = data.table::fread(
  file = here('mid_res/baseline_response/dataV3/baseline.module.name.shortened.txt'),
  sep = '\t'
)
new.names$cname2 = paste(new.names$celltype, new.names$shortname, sep = ' :: ')
colnames(ds) = new.names$cname2

# fix subject names 
dp = ds %>%
  as.data.frame() %>% 
  rownames_to_column('subject') %>% 
  mutate(subject = str_sub(subject, 1,3)) %>% 
  column_to_rownames('subject')

# read matrix and network 
net = readRDS(file = here('mid_res/baseline_response/dataV3/net.rds'))
edf = as_long_data_frame(net)

# readadj p vas for comparison 
padj = readRDS(file = here('mid_res/baseline_response/dataV3/padj.rds'))

# aes set 
plotattr = list(
  theme_bw(),
  geom_point(shape = 21, size = 3.5, stroke = 0.8), 
  theme(axis.title.y = element_text(size = 8)), 
  theme(axis.title.x = element_text(size = 8)), 
  scale_fill_manual(values = c(col.alpha("red",0.7), col.alpha("dodgerblue", 0.7))),
  theme(legend.position = 'none', legend.key.size = unit(0.29, units = 'cm')), 
  theme(aspect.ratio = 1) 
)

# specify response in vis.  
high.responders = c("205","207","209","212", "215",
                    "234","237","245","250","256")

for (i in 1:nrow(edf)) {
 
  # get edge to plot from the data framed network  
  mod.names = c(edf[i, ]$from_name, edf[i, ]$to_name)
  
  cplot = dp %>% 
    select(all_of(mod.names)) %>% 
    rownames_to_column('subject') %>% 
    mutate(response = ifelse(subject %in% high.responders, 'high', 'low'))
  
  
  ctp = cor.test(cplot[ , 2], cplot[ ,3], method = 'spearman', exact = FALSE)$p.value
  adjusted.p = padj[mod.names[1], mod.names[2]]
  print(ctp < adjusted.p)

  # calculate and vis. correlation across all subjects; color by response. 
  p = ggpubr::ggscatter(cplot, 
                        x = mod.names[1], 
                        y = mod.names[2],
                        conf.int = FALSE, 
                        color = col.alpha('white','0.01'),
                        add.params = list(color = "black", fill = "grey"))  +
    plotattr +
    aes(fill = response) 
  
  # save name by cell type first 
  modsave = str_replace_all(mod.names,pattern = ' :: ', replacement = '..')
  ggsave(p, filename = paste0(figpath, modsave[1], '___', modsave[2], 'cor.pdf'), width = 3.1, height = 3.1)
  
}

# draw a legend 
fp2 = here('mid_res/baseline_response/figuresV3/')
p2 = p + theme(legend.position = 'top')
legend <- cowplot::get_legend(p2)
pdf(file = paste0(fp2, 'LEGEND.pdf'),width = 2, height = 1)
grid::grid.draw(legend)
dev.off()

