# H1N1 differential expression testing and gene set enrichment analysis
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(variancePartition))
suppressMessages(library(scglmmr))
source("functions/analysis_functions.R")

# make output directories 
datapath = here("mid_res/1_H1N1_pseudobulk_DE/dataV4/")
dir.create(datapath, recursive = TRUE)
figpath = here("mid_res/1_H1N1_pseudobulk_DE/figuresV4/")
dir.create(figpath, recursive = TRUE)

# parallel options 
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# read processed pseudobulk data 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))

# subset to unadjuvanted cohort and remove cell type string from sample names 
pb = lapply(pb, function(x) x = x[ ,1:40])
cnames = gsub("~.*","",colnames(pb[[1]]))
pb = lapply(pb, function(x) {
  x %>% 
    as.data.frame() %>% 
    setNames(nm = cnames) %>% 
    as.matrix()
})


# sample metadata 
samplemd = readRDS(file = here('data/samplemd.rds')) %>% filter(! adjmfc.group %in% 'AS03')
samplemd$scaledage = as.vector(scale(samplemd$age))
names(samplemd)[names(samplemd) == 'sampleid'] <- 'subjectid'
names(samplemd)[names(samplemd) == 'adjmfc.group'] <- 'group'
samplemd$group = str_replace_all(string = samplemd$group,pattern = ' ', replacement = '')

# format 
samplemd = samplemd %>% 
  mutate(time.group = paste(timepoint, group,sep = "_")) %>% 
  remove_rownames() %>% 
  column_to_rownames('sample')

# relevel combined factor 
samplemd$time.group = factor(samplemd$time.group, 
                             levels = c('d0_high', 'd1_high', 'd7_high', 
                                        'd0_low', 'd1_low', 'd7_low'))

samplemd$timepoint = factor(samplemd$timepoint, levels = c('d0', 'd1', 'd7'))

# create separate model metadata for the separate cohorts being tested 
d1 = samplemd[samplemd$time_cohort == 'd1', ] %>% droplevels()
d7 = samplemd[samplemd$time_cohort == 'd7', ] %>% droplevels()
d0 = samplemd[samplemd$timepoint == 'd0', ] %>% droplevels()


# subset the bulk lists for each time cohort 
d1d = lapply(pb, function(x){ x = x[ , rownames(d1)]})
d7d = lapply(pb, function(x){ x = x[ , rownames(d7)]})
d0d = lapply(pb, function(x){ x = x[ , rownames(d0)]})


################################
# fit day 1 model 
f1 <- ~ 0 + timepoint + batch + gender + age + (1|subjectid) 

# set up contrast matrix (based on first element of list) 
d = edgeR::DGEList(counts = d1d[[1]], samples = d1)
cmat = getContrast(exprObj = d, formula = f1, data = d1, coefficient = c( 'timepointd1', 'timepointd0'))
plotContrasts(cmat)

# run on each subset 
fit1 = v1 = list()
for (i in 1:length(d1d)) {
  # init data 
  meta = d1 
  form = f1 
  contrast_matrix = cmat
  counts = d1d[[i]]
  
  # dge list 
  d = edgeR::DGEList(counts = counts, samples = meta)
  
  # filter cell type specific lowly expressed genes and calc norm factors 
  gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, design = as.factor(d$samples$timepoint))
  print(names(d1d)[i]);print(table(gtable))
  d = d[gtable, keep.lib.sizes=FALSE]
  d = edgeR::calcNormFactors(object = d)
  
  # get voom observation level weights 
  v = voomWithDreamWeights(counts = d, 
                           formula = form,
                           data = meta, 
                           BPPARAM = pparam, 
                           plot = TRUE, save.plot = TRUE)
  # fit contrast mixed model 
  fitmm = dream(exprObj = v, formula = form, data = meta,
                L = contrast_matrix, useWeights = TRUE,
                BPPARAM = pparam, REML = TRUE)
  
  fitmm = variancePartition::eBayes(fit = fitmm)
  # save results 
  v1[[i]] = v
  fit1[[i]] = fitmm
}
names(v1) = names(fit1) = names(d1d)


################################
# fit day 7 model (uses same formula)
f1 <- ~ 0 + timepoint + batch + gender + age + (1|subjectid) 

# set up contrast matrix (based on first element of list) 
d = edgeR::DGEList(counts = d7d[[1]], samples = d7)
cmat = getContrast(exprObj = d, formula = f1, data = d7, coefficient = c( 'timepointd7', 'timepointd0'))
plotContrasts(cmat)

# run on each subset 
fit7 = v7 = list()
for (i in 1:length(d7d)) {
  # init data 
  meta = d7 
  form = f1 
  contrast_matrix = cmat
  counts = d7d[[i]]
  
  # dge list 
  d = edgeR::DGEList(counts = counts, samples = meta)
  
  # filter cell type specific lowly expressed genes 
  gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, design = as.factor(d$samples$timepoint))
  table(gtable)
  d = d[gtable, keep.lib.sizes=FALSE]
  d = edgeR::calcNormFactors(object = d)
  
  # get voom observation level weights 
  v = voomWithDreamWeights(counts = d, 
                           formula = form, 
                           data = meta, 
                           BPPARAM = pparam, 
                           plot = TRUE, save.plot = TRUE)
  # fit contrast mixed model 
  fitmm = dream(exprObj = v, formula = form, data = meta,
                L = contrast_matrix, useWeights = TRUE, 
                BPPARAM = pparam, REML = TRUE)
  
  fitmm = variancePartition::eBayes(fit = fitmm)
  # save results 
  v7[[i]] = v
  fit7[[i]] = fitmm
}
names(v7) = names(fit7) = names(d7d)

# run baseline model using limma 
# set up fixed effects model to run with limma 
mod0 <- model.matrix(~ 0 + group + batch + gender + age, data = d0)
colnames(mod0) = c("high", "low", 'batch2', "genderM", "age")
c0 = makeContrasts(adjmfc = high - low, levels = colnames(mod0))

fit0 = v0 = cont0 = list()
for (i in 1:length(d0d)) {
  # init data 
  meta = d0
  # form = f1 
  contrast_matrix = c0
  counts = d0d[[i]]
  
  # dge list 
  d = edgeR::DGEList(counts = counts, samples = meta)
  
  # filter cell type specific lowly expressed genes ** Change grouping factor for filter by expression to group
  gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, design = as.factor(d$samples$group))
  table(gtable)
  d = d[gtable, keep.lib.sizes=FALSE]
  d = edgeR::calcNormFactors(object = d)
  
  # get voom observation level weights
  v = voom(counts = d, design = mod0, save.plot = TRUE, plot = TRUE)
  #v = voomWithDreamWeights(counts = d, formula = form, data = meta, BPPARAM = pparam, plot = TRUE, save.plot = TRUE)
  # fit contrast mixed model 
  fit = limma::lmFit(object = v,design = mod0)
  cfit = contrasts.fit(fit = fit, contrasts = c0)
  eb = limma::eBayes(fit = cfit)
  # save results 
  v0[[i]] = v
  fit0[[i]] = fit
  cont0[[i]] = eb
}
names(v0) = names(fit0) = names(cont0) = names(d0d)



# save model fitting data 
saveRDS(object = samplemd, file = paste0(datapath, 'samplemd.rds'))
saveRDS(object = d1, file = paste0(datapath, 'd1.rds'))
saveRDS(object = d7, file = paste0(datapath, 'd7.rds'))
saveRDS(object = d0, file = paste0(datapath, 'd0.rds'))
saveRDS(object = d1d, file = paste0(datapath, 'd1d.rds'))
saveRDS(object = d7d, file = paste0(datapath, 'd7d.rds'))
saveRDS(object = d0d, file = paste0(datapath, 'd0d.rds'))

# save model fits 
# d0
saveRDS(object = fit0, file = paste0(datapath, 'fit0.rds'))
saveRDS(object = cont0, file = paste0(datapath, 'cont0.rds'))
saveRDS(object = v0, file = paste0(datapath, 'v0.rds'))
# day 1 
saveRDS(object = fit1, file = paste0(datapath, 'fit1.rds'))
saveRDS(object = v1, file = paste0(datapath, 'v1.rds'))
# day 7
saveRDS(object = fit7, file = paste0(datapath, 'fit7.rds'))
saveRDS(object = v7, file = paste0(datapath, 'v7.rds'))

# sessioninfo
sessionInfo()