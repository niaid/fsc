suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(variancePartition))
source('functions/MattPMutils.r')
suppressMessages(library(magrittr))
# figpath
figpath = here('mid_res/variance_partition/figures_vars/'); dir.create(figpath, recursive = TRUE)
datapath = here('mid_res/variance_partition/generated_data2/'); dir.create(datapath, recursive = TRUE) 

# parallel opts
# register(SnowParam(4))
# pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)
col = list(celltype = c(
  "BC_Mem" = "lightslateblue",
  "BC_Naive" = "#2B3D26",       
  "CD103_Tcell" = "#E25822",       
  "CD14_Mono"= "red",       
  "CD16_Mono"  = "firebrick4",       
  "CD38_Bcell" = "#882D17",       
  "CD4_CD161_Mem_Tcell" = "navy",       
  "CD4_CD25_Tcell"= "#B3446C",       
  "CD4_CD56_Tcell" = "maroon1",       
  "CD4_CD57_Tcell" = "#604E97",       
  "CD4_Efct_Mem_Tcell" ="#F99379",       
  "CD4Naive_Tcell" = "#0067A5",       
  "CD8_CD161_Tcell" = "olivedrab", 
  "CD8_Mem_Tcell" = "#008856",       
  "CD8_Naive_Tcell" = "#848482",       
  "CD8_NKT" = "#C2B280",       
  "HSC" = "#BE0032",       
  "IgA_CD14_Mono" = "#A1CAF1",       
  "MAIT_Like" = "#F38400",       
  "mDC" = "#875692",       
  "NK" = "#F3C300",     
  "pDC" = "#222222"))
ccu = structure(col[[1]]) 
names(ccu) = str_replace_all(string = names(ccu), pattern = '_', replacement = ' ')
ccu2 = sapply(ccu, col.alpha, 0.5)


# load bulk data 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))
dl = list.files(path = here('mid_res/variance_partition/generated_data/'),
                pattern = '.rds', recursive = TRUE,full.names = TRUE)
dl = dl[-c(15,17)] # remove total bulk 

# get cell type names (file names)
cts = list.files(path = here('mid_res/variance_partition/generated_data/'),
                 pattern = '.rds', recursive = TRUE,full.names = FALSE)
cts = cts[-c(15,17)]
cts = str_replace_all(string = cts,pattern = 'vp.rds', replacement = '')


# read and format variance partition results 
vl = lapply(dl, readRDS)
names(vl) = cts
dl = list()
for (i in 1:length(vl)) {
  p = plotVarPart(vl[[i]])
  p$data$celltype = names(vl[i])
  d = p$data 
  dl[[i]] = d
}

# combine results across cell types 
dl[[i]] %>% head 
test = reduce(dl, .f = rbind)
test2 = test %>% select(-c(gene))
test2$celltype = factor(test2$celltype, levels = cts)

# rename metadata vars 
levels(test2$variable) =
  list(
    `response group` = 'adjmfc.group',
    sex = 'gender',
    subjectID = 'sampleid',
    timepoint = 'timepoint',
    `timepoint:response` = "timepoint:adjmfc.group",
    age = 'age',
    residuals = 'Residuals'
  )


# visualize full results 
test2$celltype = str_replace_all(test2$celltype, pattern = '_',replacement = ' ')
d = test2 %>% filter(variable %in% c('age', 'sex', 'subjectID', 'timepoint'))
d$variable = factor(d$variable,levels = c('subjectID', 'timepoint', 'age', 'sex'))

# add outlier designation 
d = d %>% 
  group_by(celltype, variable) %>%
  mutate(outlier = value > quantile(value, 0.75) + IQR(value) * 1.5) %>%
  ungroup


p = ggplot(d, aes(x = reorder(celltype,value), y = value , color = celltype, fill= celltype)) +
  facet_wrap(~variable, nrow = 1) +
  theme_bw() + 
  theme(axis.text = element_text(color = 'black', size = 10)) +
  geom_boxplot(varwidth = TRUE, 
               outlier.color = 'red', 
               outlier.alpha = 0.2, 
               outlier.shape = NA, 
               show.legend = FALSE, size = 0.3) + 
  geom_point(data = function(x) dplyr::filter_(x, ~ outlier),
             position = 'jitter', 
             shape = 21, size = 1.1, stroke = 0, alpha = 1/3,
             show.legend = FALSE) + 
  scale_color_manual(values = ccu) +
  scale_fill_manual(values = ccu2) +
  ylab('% variance explained') + xlab('') +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 18, color = 'black'),
        axis.text.y = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 16, color = 'black')) + 
  coord_flip() 
p
ggsave(p, filename = paste0(figpath, 'vpartfull_.pdf'), width = 10, height = 5.1) 
ggsave(p, filename = paste0(figpath, 'vpartfull_.png'), width = 10, height = 5.1) 



######## 
# monnocyte 
m = vl$CD14_Mono %>%
  as.data.frame() %>%
  rownames_to_column('gene')

colnames(m) = c(
  'gene',
  'response.group',
  'sex',
  'subjectID',
  'timepoint',
  'timepoint:response',
  'age',
  'residuals'
)

# rank genes by vars 
ds = m[order(desc(m$sampleid)), ]

# get sampele meta data to add to gene data 
samplemd = readRDS(file = here('data/samplemd.rds'))
mdat = pb$CD14_Mono
mdat = edgeR::cpm(mdat, log = TRUE)

# plot genes 
p = plotPercentBars(vl$CD14_Mono[c('DDX3Y', 'TMEM176B',  'STAT1','PPARGC1',    'TP53RK'), ] ) 
levels(p$data$variable) = c('response.group', 'sex', 'SubjectID', 'timepoint', 'timepoint:response', 'age', 'residuals')
levels(p$data$variable) = c('response.group', 'sex', 'subjectID', 'timepoint', 'timepoint:response', 'age', 'residuals')
p = p + 
  ggsci::scale_fill_jama(alpha = 0.9) +
  theme_bw() + 
  theme(axis.text = element_text(size = 15, color = 'black')) + 
  theme(axis.title = element_text(size = 20, color = 'black'))+ 
  theme(legend.position = 'top') +
  guides(fill=guide_legend(nrow=4, byrow=TRUE)) + 
  theme(legend.text = element_text(size = 18), legend.key.size = unit(0.8,units = 'cm')) + 
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(hjust=1)) 
ggsave(p, filename = paste0(figpath,'genesubpct.pdf'), width = 5.6, height = 5)


# genes 
mgene.highlight = c(
  'PPARGC1B',
  'TMEM176B',
  'LILRA3',
  'TP53RK',
  'PRPF19',
  'HLA-DRB5',
  'GBP2',
  'PSME2',
  'VAMP5',
  'STAT1',
  'CD69',
  'MAP3K8',
  'DDX3Y'
)
# make matrix
d2 = as.data.frame(as.matrix(t(mdat[ mgene.highlight, ])))
d2 = d2 %>%
  rownames_to_column('sid') %>% 
  separate(sid, into = c('sample', 'celltype'), sep = '~') 
d2 = full_join(d2, samplemd, by = 'sample')
d2$sampleid = str_sub(d2$sampleid, -3,-1)

# subject 
p = 
  ggplot(d2, aes(x = reorder(sampleid, TMEM176B), y = TMEM176B, fill = timepoint)) +
  theme_bw() +
  xlab('Subject ID') +
  geom_point(shape = 21, size = 2.5, fill = col.alpha('black', 0.7) , color = 'white') + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 5))  + 
  theme(axis.title = element_text(color = 'black', size = 14)) +
  theme(legend.position = c(0.74, 0.24)) +
  theme(legend.key.size = unit(0.2,units = 'cm')) +
  scale_fill_manual(values = c('black', 'black', 'black'))
p
ggsave(p, filename = paste0(figpath, 'TMEM176B.pdf'), width = 2.5, height = 2.5)


# sex 
p = 
  ggplot(d2, aes(x = gender, y = DDX3Y)) +
  theme_bw() + 
  xlab('Sex') +
  theme(axis.title = element_text(color = 'black', size = 14)) +
  geom_boxplot(show.legend = FALSE, fill = col.alpha(acol = 'black', alpha = 0.3)) 
 ggsave(p, filename = paste0(figpath, 'sex_gene.pdf'), width = 2.5, height = 2.5)

# time 
p = 
  ggplot(d2, aes(x = timepoint, y = STAT1)) +
  theme_bw() + 
  xlab('Time') +
  theme(axis.title = element_text(color = 'black', size = 14)) +
  geom_boxplot(show.legend = FALSE, fill = col.alpha(acol = 'black', alpha = 0.3)) 
ggsave(p, filename = paste0(figpath, 'timegene.pdf'), width = 2.5, height = 2.5)

# Age 
p = 
  ggplot(d2, aes(x = age, y = TP53RK)) +
  theme_bw() + 
  theme(axis.title = element_text(color = 'black', size = 14)) +
  geom_point(shape = 21, size = 2.5, fill = col.alpha('black', 0.7) , color = 'white') + 
  scale_fill_manual(values = c('black')) + 
  geom_smooth(method = 'lm', color= 'black') + 
  ggpubr::stat_cor(label.x.npc = 0.1, label.y.npc = 0.1) 
p
ggsave(p, filename = paste0(figpath, 'agegene2.pdf'), width = 2.5, height = 2.5)

# Age 2 
p = 
  ggplot(d2, aes(x = age, y = PPARGC1B)) +
  theme_bw() + 
  theme(axis.title = element_text(color = 'black', size = 14)) +
  geom_point(shape = 21, size = 2.5, fill = col.alpha('black', 0.7) , color = 'white') + 
  scale_fill_manual(values = c('black')) + 
  geom_smooth(method = 'lm', color= 'black') + 
  ggpubr::stat_cor(label.x.npc = 0.1, label.y.npc = 0.1) 
p
ggsave(p, filename = paste0(figpath, 'agegene.pdf'), width = 2.5, height = 2.5)


