suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
#source(file = here('functions/scglmmr.functions.R'))
suppressMessages(library(emmeans))
source('functions/MattPMutils.r')

# set paths 
datapath = file.path(here('mid_res/mrna/generated_data/'))
figpath = file.path(here('mid_res/mrna/figures/'))

# load baseline monocyte leadingedge index unique genes 
gs0 = readRDS(file = here('mid_res/baseline_response/dataV3/g0.sub.rds'))
li0 = LeadingEdgeIndexed(gsea.result.list = gs0, padj.threshold = 0.05)
li0 = li0$CD14_Mono

# define ifn sigs 
sig.test = list('sig' = unique(unlist(li0)))

# save combined signature genes (e2k)
sig.genes = sig.test$sig
data.table::fwrite(list(sig.genes),file = paste0(datapath,'sig.txt'), sep = '\t')

# load monocyte gated CITE-seq data from pfizer data 
s.mono = readRDS('mid_res/mrna/generated_data/s.mono.rds')
s.mono = NormalizeData(s.mono,assay = 'RNA',normalization.method = 'LogNormalize')
# define umi matrix and metadata 
umi = s.mono@assays$RNA@data
md = s.mono@meta.data
# format metadata for lme4 
md$time = factor(md$day,levels = c('0', '1', '21', '22'))
md$pt_id = factor(as.character(md$pt_id))

# module score simple average for the 3 signatures defined above. 
mscore = WeightedCellModuleScore(gene_matrix = umi, 
                                 module_list = sig.test, 
                                 threshold = 0, 
                                 cellwise_scaling = FALSE, 
                                 return_weighted = FALSE)

# combine signature scores with meta.data 
dat.fit = cbind(mscore, md)
saveRDS(dat.fit, file = paste0(datapath, 'dat.fit.mono.rds'))
dat.fit = readRDS(file = here('mid_res/mrna/generated_data/dat.fit.mono.rds'))

# note strucure with one major outlier in cell number in this dataset: 
table(dat.fit$time, dat.fit$sample_id) %>% t()

# lmer formula for the 3 signatures 
f1 = 'sig ~ 0 + time + (1|pt_id)'

library(emmeans)
# fit model for each signature
m1 = lme4::lmer(formula = f1,data = dat.fit)
emm1 = emmeans::emmeans(m1,specs = ~time, lmer.df = 'asymptotic')

# contrast time differences 
clevels = levels(dat.fit$time)

#make custom contrasts 
c0 = c(1, 0, 0 ,0)
c1 = c(0, 1, 0 ,0)
c3 = c(0, 0, 1 ,0)
c4 = c(0, 0, 0 ,1)
contrast_list = list( "time1vs0" = c1 - c0, 'time22vs21' = c4 - c3 )
clist = list( 'sig' = emm1 )

c.res = 
lapply(clist, function(x) { 
  emmeans::contrast(object = x, method = contrast_list) %>% 
    broom::tidy()
  } ) %>% 
  bind_rows(.id = 'signature')
data.table::fwrite(x = c.res, file = paste0(datapath,'c.res.mono.txt'), sep = '\t')  

# delta contrast 
cmat = emmeans::contrast(object = emm1, method = contrast_list)
pairs(cmat, reverse = TRUE)
# contrast              estimate      SE  df z.ratio p.value
# time22vs21 - time1vs0   0.0967 0.00642 Inf  15.065  <.0001
cmat
# contrast   estimate      SE  df z.ratio p.value
# time1vs0     0.0846 0.00492 Inf  17.203  <.0001
# time22vs21   0.1813 0.00414 Inf  43.835  <.0001

# plotsingle cell distributionn and emmeans contrasts 
em_aes = list(theme_bw(), 
              coord_flip(), 
              theme(axis.title.y = element_text(size = 7),
                    axis.text.y = element_text(size = 7)), 
              scale_color_manual('grey')
              )

plot.aes = list(theme_bw(), ylab(label ='Baseline high responder\nCD14 Mono signature'))


cu = sapply(c('grey', '#e2a359', 'grey', '#e2a359'), col.alpha, 0.8) %>% unname()
# combined signature change emm in p1 and change y value in p0
p0 = ggplot(dat.fit, aes(x = time, y = sig, fill = time)) + 
  geom_violin(show.legend = F) + 
  plot.aes +
  xlab('time') + 
  scale_fill_manual(values = cu)+
  theme(axis.title.x = element_text(size = 12))
ggsave(p0, filename = paste0(figpath, 'sig.pdf'), width = 2.5, height = 3)
p1 = plot(emm1) + em_aes
ggsave(p1, filename = paste0(figpath, 'sig.emm.pdf'), width = 1, height = 3)



