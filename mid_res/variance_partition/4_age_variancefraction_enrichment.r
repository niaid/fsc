suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(variancePartition))
source('functions/MattPMutils.r')
library(magrittr)
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

# set theme 
mtheme = list(
  theme_bw(), 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
)

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


# rank genes by variance fraction assciated with age 
hlmk = readRDS(file = here('signature_curation/hallmark.rds'))
# parallel opts
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)


agelist = slist = list()
#vars = c('age', 'timepoint', 'gender', 'sampleid')
for (u in 1:length(vl)) {
  
  # get variance fractions
  m = vl[[u]] %>% 
    as.data.frame() %>%  
    rownames_to_column('gene') 
  
  # rank genes by sex 
  da = m[order(desc(m$age)), ] 
  rank.age = structure(da$age, names= da$gene)
  
  # rank genes by subject 
  ds = m[order(desc(m$sampleid)), ] 
  rank.subject = structure(ds$sampleid, names= ds$gene)
  
  # format lists 
  agelist[[u]] = rank.age
  slist[[u]] = rank.subject
}

# run gsea for age and subject 
age.gs = scglmmr::FgseaList(rank.list.celltype = age.gsea,pathways = hlmk, scoreType = "pos", BPPARAM= pparam)

age.gsea = sgsea = list()
for (u in 1:length(agelist)) {
  age.gsea[[u]] = fgsea::fgsea(hlmk, agelist[[u]], scoreType = "pos",  BPPARAM = pparam)
  sgsea[[u]] = fgsea::fgsea(hlmk, slist[[u]], scoreType = "pos",  BPPARAM = pparam)
}
for (u in 1:length(agelist)) {
  age.gsea[[u]]$celltype = names(vl)[u]
}
p = scglmmr::PlotFgsea(gsea_result_list = age.gsea, padj_filter = 0.05 )


age.gsea.sub = lapply(age.gsea, function(x)
  x %>%  filter(
    pathway %in% c(
      'HALLMARK_IL2_STAT5_SIGNALING',
      'HALLMARK_IL6_JAK_STAT3_SIGNALING',
      'HALLMARK_INFLAMMATORY_RESPONSE',
      'HALLMARK_ALLOGRAFT_REJECTION'
    )
  ))

li = scglmmr::LeadingEdgeIndexed(gsea.result.list = age.gsea.sub,padj.threshold = 0.05)
li = Filter(li, f = length)
li$CD8_CD161_Tcell


dage = bind_rows(age.gea.sub, .id = 'celltype')
dage = dage %>%  filter(padj < 0.05)



dage$pathway = str_replace_all(dage$pathway,pattern = 'HALLMARK_',replacement = '')
dage$pathway = str_replace_all(dage$pathway,pattern = '_',replacement = ' ')
cu = c("olivedrab", "#848482")
cu2 = sapply(cu, col.alpha, 0.5)
p = 
  ggplot(dage %>% 
           filter(celltype %in% c('CD8_Naive_Tcell', 'CD8_CD161_Tcell')), 
         aes(x = NES, y = pathway, label = celltype, 
             group = celltype, fill = celltype, color = celltype)) + 
  mtheme + 
  geom_linerange(aes(x = NES, color = celltype, xmin = 0, xmax = NES),
                 position = position_dodge(width = 0.35)) +
  geom_point(aes(x = NES, color = celltype),position = position_dodge(width = 0.35)) +
  geom_point(shape = 21, size = 2.5, position = position_dodge(width = 0.35)) +
  ggsci::scale_fill_npg() + 
  theme(axis.text = element_text(color = 'black')) + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  ylab('') + 
  xlab('Normalized Enrichment Score') + 
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle('Age associated variance enrichment') 

p
ggsave(p,filename = paste0(figpath, 'tcell_age.pdf'), width = 6, height =2.2)




# CD8 Naive 
# logcpm matrix
cd8.genes = unique(unlist(li$CD8_Naive_Tcell, use.names = FALSE))
mdat = pb$CD8_Naive_Tcell
mdat = edgeR::cpm(mdat, log = TRUE)
d2 = as.data.frame(as.matrix(t(mdat[cd8.genes,])))
d2 = d2 %>% 
  rownames_to_column('sid') %>% 
  separate(sid, into = c('sample', 'celltype'), sep = '~')
d2 = full_join(d2, samplemd, by = 'sample')
d2$sampleid = str_sub(d2$sampleid, -3, -1)

# scale age
scale.simple = function(x) {
  (x - mean(x)) / sd(x)
}
d2$age = scale.simple(d2$age)

# fit models
dmat = d2 %>% select(age, all_of(cd8.genes))
age.scaled = d2$age

age.coef = apply(dmat[, 2:ncol(dmat)], MARGIN = 2, function(x) {
  y = lm(x ~ 0 + age.scaled)
  return(y)
})

age.res = lapply(age.coef, broom::tidy)
age.res = bind_rows(age.res,.id = 'gene')

age.pos = age.res %>% 
  filter(estimate > 0) %$% 
  gene
plot(agelist[[11]][age.pos])

age.var = data.frame(age.var = agelist[[11]][age.pos]) %>% 
  rownames_to_column('gene')

p = 
  ggplot(age.var, aes(y = reorder(gene, age.var) , x = age.var*100, label = gene )) +
  mtheme + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab('') + 
  xlab('% variance explained by age') +
  xlim(c(-5, 40)) +
  geom_point(color = "#4DBBD5FF") +
  ggrepel::geom_text_repel(data = age.var %>% 
                             filter(gene %in% c('KLF6', 'RGS16','TNF', 'IL12A', 'HLA-DQA1', 
                                                'ABCA1', 'CD70','GZMB', 'IL10A', "FASLG", 
                                                'CCL5', 'NOD2', 'CD74', 'IFNG', 'CCL4', 'FAS')),
                           size = 3,
                           force        = 1,
                           nudge_x      = -20,
                           direction    = "y",
                           hjust        = 1,
                           segment.size = 0.001) + 
  ggtitle('Positive association with age\nCD8 Naive T cells') + 
  theme(axis.title = element_text(size = 18))

ggsave(p,filename = paste0(figpath, 'cd8n.varexp.age.positive.pdf'), width = 4, height = 4)




# CD8 CD161 
# logcpm matrix
cd161.genes = unique(unlist(li$CD8_CD161_Tcell,use.names = FALSE))
mdat = pb$CD8_CD161_Tcell
mdat = edgeR::cpm(mdat, log = TRUE)
d2 = as.data.frame(as.matrix(t(mdat[cd161.genes, ])))
d2 = d2 %>%
  rownames_to_column('sid') %>%
  separate(sid, into = c('sample', 'celltype'), sep = '~') 
d2 = full_join(d2, samplemd, by = 'sample')
d2$sampleid = str_sub(d2$sampleid, -3,-1)

# scale age 
scale.simple = function(x) { (x - mean(x)) / sd(x)}
d2$age = scale.simple(d2$age)

# fit models 
dmat = d2 %>% select(age, all_of(cd161.genes))
age.scaled = d2$age

age.coef = apply(dmat[ ,2:ncol(dmat)],MARGIN = 2, function(x) { 
  y = lm(x ~ 0 + age.scaled)
  return(y)
} )

age.res = lapply(age.coef, broom::tidy)
age.res = bind_rows(age.res,.id = 'gene')

age.pos = age.res %>%  filter(estimate > 0) %$% gene
names(vl)
age.var = data.frame(age.var = agelist[[9]][age.pos]) %>%  rownames_to_column('gene')
p = 
  ggplot(age.var, aes(y = reorder(gene, age.var) , x = age.var*100, label = gene )) +
  mtheme + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab('') + 
  xlab('% variance explained by age') +
  xlim(c(-5, 40)) +
  geom_point(color = "#E64B35FF") +
  ggrepel::geom_text_repel(data = age.var %>% 
                             filter(age.var > 0.08) %>% 
                             filter(gene %in% c('HLA-DQA1', 'KRT1', 'IFNG', 'CCL4', 'KLRD1', 
                                                "CCL5", "CFP", "CD74", "HLA-DRA",  "IL2RB",  
                                                "IL17RA",  "TNFRSF8", "MAP3K8",   "SERPINB6", "CD38",   
                                                "TYK2", "CD70",     "PROK2" ,   "RGS1"    )),
                           size = 3,
                           force        = 1,
                           nudge_x      = -20,
                           direction    = "y",
                           hjust        = 1,
                           segment.size = 0.001) + 
  ggtitle('Positive association with age\nCD8 CD161+ T cells') + 
  theme(axis.title = element_text(size = 18))
ggsave(p,filename = paste0(figpath, 'cd161T.varexp.age.positive.pdf'), width = 4, height = 4)



# write lists 
plot(age.coef)
pos.age = age.coef[age.coef > 0]
neg.age = age.coef[age.coef < 0]
data.table::fwrite(list(names(pos.age)),file = paste0(datapath,'pos.age.cd8n.txt'))
data.table::fwrite(list(names(neg.age)),file = paste0(datapath,'neg.age.cd8n.txt'))



