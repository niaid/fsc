# create baseline expression leading edge module correlation matrix 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))
# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")

# read pb data, subset to day 0 non adj, subset out day 0 metadata. 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))
pb = lapply(pb, function(x) x = x[ ,1:40])
cnames = gsub("~.*","",colnames(pb[[1]]))
pb = lapply(pb, function(x){
  x %>% as.data.frame() %>% setNames(nm = cnames) %>% as.matrix() 
  })
d0 = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV3/d0.rds'))
d0d = lapply(pb, function(x){ x = x[ , rownames(d0)]})


# convert pb data to log counts per million
d.norm = list()
for (i in 1:length(d0d)) {
  d = edgeR::DGEList(counts = d0d[[i]], samples = d0)
  gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, 
                               design = as.factor(d$samples$group))
  d = d[gtable, keep.lib.sizes=FALSE]
  d.norm[[i]] = edgeR::cpm(y = d, log = TRUE, prior.count = 1)
}
names(d.norm) = names(d0d)

# get leading edge genes from cur. baseline mods 
g0.sub = readRDS(file = here('mid_res/baseline_response/dataV3/g0.sub.rds'))
li.g0 = LeadingEdgeIndexed(gsea.result.list = g0.sub, padj.threshold = 0.05)
li.g0 = base::Filter(length, li.g0)

# subset normalized expression to subsets with baseline enrichments 
d.norm = d.norm[names(li.g0)]

res = list()
for (i in 1:length(d.norm)) {
  stopifnot(all.equal( names(d.norm[i]), names(li.g0[i]) ))
  zscore = scglmmr::calc_avg_module_zscore(
    module.list = li.g0[[i]], average.data.frame = d.norm[[i]]
  )
  rownames(zscore) = paste(rownames(zscore), names(d.norm[i]), sep = '~')
  res[[i]] = zscore
}

ds = do.call(rbind, res) %>% t()
saveRDS(ds, file = paste0(datapath, 'ds.rds'))

sessionInfo()