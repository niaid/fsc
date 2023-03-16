# make tables 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))

datapath = file.path(here('mid_res/data_write/generated_data/')); 
dir.create(datapath)

# specify order of variables in the output table for readability 
var.order = c('contrast', 'celltype', 'pathway', 'NES', 'padj', 'leadingEdge')

# result format for gsea results
filter.gsea = function(list){
  lapply(list, function(x) x %>%  filter(padj < 0.05))
}

format.result = function(x) { 
  x %>% 
    select(all_of(var.order), everything()) %>%
    arrange(celltype, NES) %>% 
    tibble::remove_rownames()
}



# Baseline
# res text gsea curate
# gsea res raw 
g0.sub = readRDS(file = here("mid_res/baseline_response/dataV3/g0.sub.rds"))
g0.sub = do.call(rbind, g0.sub)
g0.sub$contrast = 'baseline high vs low responders'
d0.res = format.result(g0.sub) %>%
  select(-c(signal)) %>% 
  mutate(model = 'gene ~ 0 + response.group + batch + sex + age')
  

# day 1 non-adjuvanted vaccine 
g1c = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/g1c.rds'))
g1c = do.call(rbind, g1c)
g1c$contrast = '24h vs baseline unadjuvanted vaccine'
g1c$model = 'gene ~ 0 + timepoint + batch + sex + age + (1|subjectid) '

# day 7 non-adjuvanted vaccine
g7f = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/g7f.rds'))
g7f = lapply(g7f, function(x) 
  x %>%  
    filter(!str_sub(pathway, 1,5) == 'REACT' ) %>% 
    filter(NES > 0) %>%  
    filter(pval <0.1)
)
g7f = do.call(rbind, g7f)
g7f$contrast = 'day 7 vs baseline unadjuvanted vaccine'
g7f$model = 'gene ~ 0 + timepoint + batch + sex + age + (1|subjectid) '

# as03 model 
gc = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/gsea/gc.rds'))
gc = lapply(gc, function(x) x %>%  filter(padj < 0.1))
gc = do.call(rbind, gc)
gc = format.result(gc)
gc$contrast = 'AS03 vs unadjuvanted vaccine day 1 vs baseline fold change difference'
gc = format.result(gc) %>% 
  mutate(model = '0 + timepoint_vaccinegroup + sex + age')

# combine 
d = rbind(g1c, g7f, gc, d0.res)

# write results 
data.table::fwrite(d,
                   file = paste0(datapath,'combined.results.fsc.txt'), 
                   sep = '\t')


# gene signatures 
core_sigs = readRDS(file = here('signature_curation/sig_test_sub.rds'))
mann = data.table::fread(file = here('signature_curation/sig_test_sub_annotation.txt'))
dmod = mann
dmod$signature_genes = core_sigs

# add natural adjuvant signautres   
mv = readRDS(file = here('mid_res/vand/generated_data/mv.rds'))
dcv = readRDS(file = here('mid_res/vand/generated_data/dcv.rds'))
mono.sig = mv$leadingEdge  %>% unlist() %>% unique() 
dc.sig = dcv$leadingEdge  %>% unlist() %>% unique()
as03.sig = list(mono.sig,dc.sig)

dcite= data.frame(
  pathway = c("CD14_Mono_AS03", "mDC_AS03"), 
  annotation = c(rep('CITE-seq contrast model and sorted cell validated', 2))
)
dcite$signature_genes = as03.sig
# combine with signatures 
dmod = rbind(dmod, dcite)

data.table::fwrite(dmod,
                   file = paste0(datapath,'combined.modules.fsc.txt'), 
                   sep = '\t')


### Variance fractions 
vp = readRDS(file = here('mid_res/variance_partition/generated_data/vp.rds'))
vp = as.data.frame(vp) %>% 
  rownames_to_column('gene') %>% 
  select(gene, everything()) %>% 
  mutate(model = '~ age + (1|sex) + (1|subjectid) + (1|celltype) + (1|timepoint) + (1|adjmfc.group) + (1|celltype:timepoint)')
# write
data.table::fwrite(vp,
                   file = paste0(datapath,'variance.partition.across.celltypes.txt'), 
                   sep = '\t')

# per cell type results 
#pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))
dl = list.files(
  path = here('mid_res/variance_partition/generated_data/'),
  pattern = '.rds',
  recursive = TRUE,
  full.names = TRUE
)
dl = dl[-c(15,17)] # remove total bulk 

# get cell type names (file names)
cts = list.files(
  path = here('mid_res/variance_partition/generated_data/'),
  pattern = '.rds',
  recursive = TRUE,
  full.names = FALSE
)
cts = cts[-c(15,17)]
cts = str_replace_all(string = cts,pattern = 'vp.rds', replacement = '')
# read and format variance partition results 
vl = lapply(dl, readRDS)
for (i in 1:length(vl)) {
  vl[[i]] = vl[[i]] %>% 
    as.data.frame() %>% 
    rownames_to_column('gene') %>% 
    mutate(model_fit_within_this_celltype = names(vl)[i]) %>% 
    mutate(model = '~ age + (1|sex)  + (1|subjectid) + (1|timepoint) + (1|adjmfc.group) + (1|timepoint:adjmfc.group)')
}
vp_within = do.call(rbind, vl)
data.table::fwrite(vp_within,
                   file = paste0(datapath,'variance.partition.within.celltypes.txt'), 
                   sep = '\t')



