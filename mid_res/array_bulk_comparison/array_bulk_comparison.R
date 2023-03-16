# R 4.0.5 
suppressMessages(library(tidyverse))
suppressMessages(library(scglmmr))
suppressMessages(library(here))
suppressMessages(library(BiocParallel))
suppressMessages(library(edgeR))
suppressMessages(library(variancePartition))

# set paths 
figpath = here('mid_res/array_bulk_comparison/figures/')
dir.create(figpath)
datapath = here('mid_res/array_bulk_comparison/generated_data/')
dir.create(datapath)

# set parallel options 
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

#pb = readRDS(file = here('mid_res/pb.ds'))
#fit7 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/fit7.rds'))

core7 = readRDS(file = here('signature_curation/core_d7.rds'))
core7$`LI.M156.0 plasma cell b cell Ig`


# read processed pseudobulk data 
# subset to unadjuvanted cohort and remove cell type string from sample names 
s = readRDS(here('data/h1h5_annotated_with_meta.rds'))


# day1.cohort = c("200", "205","207", "236", "237", "250", "273" ,"279")
day7.cohort = c("201", "209", "212", "215", "229", "233",
                "234", "245", "256",  "261",  "268", "277")

md7 = s@meta.data %>%  
  filter(cohort == 'H1N1') %>% 
  filter(sampleid %in% day7.cohort) %>% 
  arrange(sampleid, timepoint)
umi7 = s@raw.data[ ,rownames(md7)]

# pseudobulk all cells
scell = lapply(X = split(md7, f = md7$sample), FUN = rownames)
csample = lapply(scell, function(x) Matrix::rowSums(umi7[ ,x]))
pbmat = as.data.frame(t(do.call(cbind, csample))) %>% t()


#define day 7 metadata 
met = md7 %>% 
  select(sample,timepoint , subjectid = sampleid) %>% 
  group_by(sample, subjectid, timepoint) %>% 
  distinct() %>%
  ungroup() %>% 
  column_to_rownames('sample')
met$timepoint = factor(met$timepoint, levels = c('d0', 'd7'))

# check order 
stopifnot(isTRUE(all.equal(colnames(pbmat), rownames(met))))

# filter features  
gene.keep = filterByExpr(y = pbmat, design = met$timepoint,min.count = 3)
pbmat = pbmat[gene.keep, ]

# fit model 
f1 = ~ 0 + timepoint + (1|subjectid)
L1 = makeContrastsDream(formula = f1, data =  met,
                        contrasts =  "timepointd7 - timepointd0")
v7 = voomWithDreamWeights(counts = pbmat,formula = f1,
                          BPPARAM = pparam,data = met)
result7 = dream(exprObj =  v7,formula = f1,data = met,
                L = L1, BPPARAM = pparam, useWeights = TRUE)
# save 
saveRDS(result7,file = paste0(datapath,'result7.rds'))


## Part II fit same model on microarray data 

# array data coefficient 
# read array data 
array = data.table::fread("data/CHI_H1N1_data/microarray/CHI_GE_matrix_gene.txt", data.table = F) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("gene") %>% 
  select(-matches("day70")) %>% 
  select(., matches("day0|day7")) %>% 
  select(-matches("pre")) 

# day 7 samples; no day 7 data for subject 209 
array7 =
  array %>% 
  select(which(substr(names(.),1,3) %in% day7.cohort)) %>% 
  select(-matches("209"))

# Metadata 
d7md  = 
  colnames(array7) %>% 
  base::as.data.frame() %>% 
  rename(sample = ".") %>% 
  mutate(timepoint = str_sub(sample, -4,-1)) %>% 
  mutate(subjectid = str_sub(sample, 1,3)) %>% 
  column_to_rownames("sample")

d7md$timepoint = factor(d7md$timepoint, levels = c('day0', 'day7'))

# check order 
stopifnot(isTRUE(all.equal(colnames(array7), rownames(d7md))))

# test same genes
gene.sub = rownames(pbmat)
array7 = array7[gene.sub, ]
gene.keep2 = !is.na(Matrix::rowSums(array7))
array7 = array7[gene.keep2, ]


# fit model
L1.1 = makeContrastsDream(formula = f1, data =  d7md,
                          contrasts =  "timepointday7 - timepointday0")
# no weights for normalized microarray data ; same formula
result7.1 = dream(exprObj =  array7, formula = f1,data = d7md,
                L = L1.1, BPPARAM = pparam, useWeights = FALSE)
# save 
saveRDS(result7.1,file = paste0(datapath,'result7.1.rds'))


# comparison 

# extract results for array data 
ra = ExtractResult(model.fit.list = list('array' = result7.1), 
                   coefficient.number = 1,
                   coef.name = 'timepointday7 - timepointday0')
ra$array$logFC.array = ra$array$logFC
ra = ra$array %>%  
  select(gene,  logFC.array)
# extract results for CITE-seq bulk data 
rc = ExtractResult(model.fit.list = list('CITE-seq Bulk' = result7),
                   coefficient.number = 1, 
                   coef.name = 'timepointd7 - timepointd0')
rc$`CITE-seq Bulk`$logFC.CITEseq = rc$`CITE-seq Bulk`$logFC
rc = rc$`CITE-seq Bulk` %>% 
  select(gene,  logFC.CITEseq)


# visualize correlation between signals 
d = full_join(ra, rc)
dsub = d %>% filter(gene %in% core7$`LI.M156.0 plasma cell b cell Ig`)
p = 
  ggplot(dsub, aes(x = logFC.array , y = logFC.CITEseq)) + 
  theme_bw() + 
  geom_smooth(method = 'lm', color = 'black') + 
  geom_point(shape = 21, color = 'white', fill = 'black') + 
  ggrepel::geom_text_repel(data = dsub, mapping = aes(label = gene), 
                           segment.size = 0.1, box.padding = 0.1, 
                           max.overlaps = 4,  size = 2.6) + 
  ggpubr::stat_cor(method = "pearson") + 
  ylab("LI.M156 Module \n CITE-seq Day 7 log2 FC") + 
  xlab("LI.M156 Module \n Microarray Day 7 log2 FC") 
ggsave(p, filename = paste0(figpath, "/m156_gene_correlation.pdf"),width = 3, height = 3)


