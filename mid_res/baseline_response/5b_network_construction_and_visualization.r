# baseline nettwork correlation 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
source('functions/MattPMutils.r')

# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")

# load SLI corrected matrix 
mat = readRDS(file = here('mid_res/baseline_response/dataV3/mat.rds'))

# Create mDC and innate sub-network 
dn = data.frame(mods = rownames(mat)) %>%
  mutate(name = mods) %>% 
  separate(name, into = c('celltype', 'module'),sep = ' :: ')
innate.sub = dn %>% 
  filter(celltype %in% c('CD14_Mono', 'CD16_Mono', 'MAIT_Like', 'mDC', 'BC_Naive'))
ms = innate.sub$mods

# rm QCd modules 
#m.rm = readRDS(here('mid_res/baseline_response/dataV3/m.rm.rds'))
m.rm = c(
  "CD14_Mono :: M111.1 viral sensing IRF2",
  "MAIT_Like :: Kegg Ag Presentation",           
  "MAIT_Like :: reactome interferon alpha beta"
)
ms = ms[!ms %in% m.rm]

# calculate adjusted p values across the spearman correlation matrix 
# load uncorrected matrix Hmisc obejct containing correlation p values 
spearmanmat = readRDS(file = here('mid_res/baseline_response/dataV3/spearmanmat.rds'))
padj = p.adjust.cormat(hmisc.cor =  spearmanmat, method = 'fdr')
saveRDS(padj, file = paste0(datapath,'padj.rds'))

# filter to the innate subnetwork 
# filter edges in SLI adjusted network
# only include correlations with adjusted p < 0.05
mat2 = mat
mat2 = mat2[ms, ms]
padj = padj[ms, ms]
stopifnot(isTRUE(all.equal(colnames(padj), colnames(mat2))))
# filter based on adjusted p 
mat2[padj > 0.05] <- 0

# save the pruned mat2 
saveRDS(mat2,file = paste0(datapath,'mat2.rds'))
mat2 = readRDS(file = here('mid_res/baseline_response/dataV3/mat2.rds'))


# make graph of the strongly linked edges pruned above
net <- graph_from_adjacency_matrix(
  mat2, weighted = TRUE,
  mode = 'undirected',
  diag = FALSE
  )


# prune the graph further to retain links above the median weight
med.weight <- median(E(net)$weight)
mat3 = mat2
mat3[mat3 < med.weight] <- 0
saveRDS(mat3,file = paste0(datapath,'mat3.rds'))
mat3 = readRDS(file = here('mid_res/baseline_response/dataV3/mat3.rds'))

# make a subhraph with stonger connections above prev median weight.
net <- graph_from_adjacency_matrix(
  mat3,
  weighted = TRUE,
  mode = 'undirected',
  diag = FALSE
  )

# create network annotations frame for vertices 
d = data.frame(signal = colnames(mat3)) %>%
  mutate(s = signal) %>%
  separate(s,into = c('celltype', 'module'),sep = ' :: ')

# specify vertex attributes 
V(net)$celltype = d$celltype
V(net)$module = d$module

# calculate network degree and hubs / authority scores (same for undirected)
V(net)$degree <- degree(net)                        
V(net)$hubs <- hub.score(net)$vector                
V(net)$authorities <- authority.score(net)$vector   

# add vertex property information to `d``
d$degree = V(net)$degree
d$hubscore = V(net)$hubs

# add d7 correlation to d 
d7.cor = readRDS(file = here('mid_res/baseline_response/dataV3/d7.cor.rds'))
d7.cor.p = d7.cor$P[1, -1][ms]
d7.cor.rho = d7.cor$r[1, -1][ms]
# check orders correct 
stopifnot(isTRUE(all.equal(names(d7.cor.p), d$signal)))
d$d7cor.p = d7.cor.p
d$d7cor.rho = d7.cor.rho

# add this information to the network 
V(net)$d7cor.p = d7.cor.p
V(net)$d7cor.rho = d7.cor.rho

# specify edge width as the weight 
E(net)$width = E(net)$weight 

# save network 
saveRDS(net,file = paste0(datapath,'net.rds'))
net = readRDS(file = here('mid_res/baseline_response/dataV3/net.rds'))

############################
# plot hubs 
signal.highlight = d %>% filter(celltype == 'CD14_Mono') %$% signal
signal.highlight2 = d %>% filter(celltype == 'CD16_Mono') %$% signal

cu3 = c('#FFD38F', '#F4A69B', '#A7DDEA', '#8ACFC3', '#9FABC4')

p = 
  ggplot(d, aes(y = reorder(module, hubscore), x = hubscore ,  fill = celltype, label = signal )) + 
  theme_bw() + 
  geom_point(shape =21, size = 3.5) + 
  #ggsci::scale_fill_npg(alpha = 0.8) + 
  scale_fill_manual(values = cu3) + 
  ylab('') + 
  theme(legend.position = c(0.8,0.15), legend.key.size = unit(0.2,units = 'cm')) +
  theme(axis.text = element_text(color = 'black')) +
  xlab('Hub Score') + 
  theme(aspect.ratio = 1.1)  +
  ggrepel::geom_text_repel(data = d %>% filter(signal %in% signal.highlight & hubscore > 0.75), 
                           size = 2.5, nudge_y = 0, nudge_x = -0.3, seed = 2, segment.size = 0.1,
                           force = 40,
                           max.overlaps = 10) + 
  ggrepel::geom_text_repel(data = d %>% filter(signal %in% signal.highlight2 & hubscore > 0.85), 
                           size = 2.5, nudge_y = 0, nudge_x = -0.3, box.padding = 0.4, seed = 1, 
                           max.overlaps =10,force = 40,
                           segment.size = 0.1) 
p
ggsave(p, filename = paste0(figpath, 'innate.subnetwork.hub.pdf'), width = 6.5, height = 6.5)

# specify colors for nodes 
cu = c( col.alpha('orange', alpha = 0.5), ggsci::pal_npg(alpha = 0.5)(4))
c.celltype = cu[factor(V(net)$celltype)]
# layout network.
lay <- layout_in_circle(net)

# specify celltypes for leend 
cts = str_replace_all(string = levels(factor(V(net)$celltype)),pattern = '_',replacement = ' ')


# version with vertiices highlighted 
# specify sve path for subgraph plots 
figpath3 = here('mid_res/baseline_response/figuresV3/subgraphsLABELED/'); dir.create(figpath3, recursive = TRUE)
# plot the subgraphs 
for (i in 1:length(unique(V(net)))) {
  
  # highlight edges 
  # specify subset highlighted 
  edge.highlight.t = incident(net, v = V(net)[i], mode="all")
  # for savig 
  signal = str_replace_all( names(V(net)[i]), pattern = ' :: ', replacement = '  ')
  
  # make a new network 
  net.sp = net 
  # set size for highlighted edge 
  E(net.sp)$width = 1.1
  # remove the other edges for visualization 
  ot = E(net.sp)[!E(net.sp) %in% edge.highlight.t]
  net.sp <- delete_edges(net.sp, edges = ot)
  E(net.sp)$color = col.alpha('black',0.7)
  # plot network 
  pdf(file = paste0(figpath3,signal,'.subnetwork.pdf'),width = 10, height = 10)
  plot(net.sp, 
       layout = lay, 
       vertex.label = names(V(net.sp)),
       vertex.label.cex=0.3, 
       edge.size = E(net.sp)$weight,
       vertex.size = log(V(net.sp)$degree+ 1)*4,
       edge.curved = 0.1,
       vertex.color = c.celltype,
       vertex.size=degree(net.sp)) 
  legend(x=0.75, y=1.2, 
         cts,
         pch=21, 
         pt.bg=cu,
         pt.cex=1,  cex=.5, bty="n",ncol=1)
  dev.off()
  
}

#### This not run in published workflow 
### Commented out for published workflow bc redundant w code below -- this used to make figs. 
# specify sve path for subgraph plots 
# figpath2 = here('mid_res/baseline_response/figuresV3/subgraphs/'); dir.create(figpath2, recursive = TRUE)
# # plot the subgraphs 
# for (i in 1:length(unique(V(net)))) {
#   
#   # highlight edges 
#   # specify subset highlighted 
#   edge.highlight.t = incident(net, v = V(net)[i], mode="all")
#   # for savig 
#   signal = str_replace_all( names(V(net)[i]), pattern = ' :: ', replacement = '  ')
# 
#   # make a new network 
#   net.sp = net 
#   # set size for highlighted edge 
#   E(net.sp)$width = 1.1
#   # remove the other edges for visualization 
#   ot = E(net.sp)[!E(net.sp) %in% edge.highlight.t]
#   net.sp <- delete_edges(net.sp, edges = ot)
#   E(net.sp)$color = col.alpha('black',0.7)
#   # plot network 
#   pdf(file = paste0(figpath2,signal,'.subnetwork.pdf'),width = 4.5, height = 4.5)
#   plot(net.sp, 
#        layout = lay, 
#        vertex.label = NA,
#        vertex.label.cex=0.3, 
#        edge.size = E(net.sp)$weight,
#        vertex.size = log(V(net.sp)$degree+ 1)*4,
#        edge.curved = 0.1,
#        vertex.color = c.celltype,
#        vertex.size=degree(net.sp)) 
#   legend(x=0.75, y=1.2, 
#          cts,
#          pch=21, 
#          pt.bg=cu,
#          pt.cex=1,  cex=.5, bty="n",ncol=1)
#   dev.off()
# }
# 