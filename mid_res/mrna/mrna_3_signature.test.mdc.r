suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
suppressMessages(library(emmeans))

# set save paths 
datapath = file.path(here('mid_res/mrna/generated_data/'))
figpath = file.path(here('mid_res/mrna/figures/'))

# load baseline monocyte leadingedge index unique genes 
gs0 = readRDS(file = here('mid_res/baseline_response/dataV3/g0.sub.rds'))
li0 = LeadingEdgeIndexed(gsea.result.list = gs0, padj.threshold = 0.05)
li0 = li0$mDC
# define sigs 
li.0.m11  = grepl('LI.M11.0',x = names(li0))
li.0_m11 = li0[li.0.m11]
li0.other = li0[!li.0.m11]
m11.sig = li.0_m11 %>%  unlist() %>%  unique()
non.m11.sig = li0.other %>% unlist() %>% unique()
# further prune non ifn to not include overlappint genes with ifn sigs.
both = intersect(m11.sig, non.m11.sig)
non.m11.sig = non.m11.sig[!non.m11.sig %in% both]
sig.test = list('msig' = m11.sig, 'non.msig.sig' = non.m11.sig, 'sig' = unique(unlist(li0)))


# load monocyte gated CITE-seq data from pfizer data 
s.mdc = readRDS('mid_res/mrna/generated_data/s.mdc.rds')
s.mdc = NormalizeData(s.mdc,assay = 'RNA', normalization.method = 'LogNormalize')
# define umi matrix and metadata 
umi = s.mdc@assays$RNA@data
md = s.mdc@meta.data
# format metadata for lme4 
md$time = factor(md$day,levels = c('0', '1', '21', '22'))
md$pt_id = factor(as.character(md$pt_id))

# this caused a singular fit 
# add log10 n cell as covariate per sample 
# samplen = md %>% select(sx = sample_id) %>% group_by(sx) %>%  tally() %>% mutate(logncell = log10(n))
# md$log10ncell = plyr::mapvalues(x = md$sample_id, from = samplen$sx,to = samplen$logncell)

# module score for the 3 signatures defined above. 
mscore = WeightedCellModuleScore(gene_matrix = umi, 
                                 module_list = sig.test, 
                                 threshold = 0, 
                                 cellwise_scaling = FALSE, 
                                 return_weighted = FALSE)

# combine signature scores with meta.data 
dat.fit = cbind(mscore, md)
saveRDS(dat.fit, file = paste0(datapath, 'dat.fit.mdc.rds'))
dat.fit = readRDS(file = here('mid_res/mrna/generated_data/dat.fit.mdc.rds'))

# note strucure not as bad on outlier sample as mono
table(dat.fit$time, dat.fit$sample_id) %>% t()

# lmer formula for the 3 signatures 
f1 = 'msig ~ 0 + time + (1|pt_id)'
f2 = 'non.msig.sig ~ 0 + time + (1|pt_id)'
f3 = 'sig ~ 0 + time + (1|pt_id)'

library(emmeans)
# fit model for each signature
m1 = lme4::lmer(formula = f1, data = dat.fit)
emm1 = emmeans::emmeans(m1, specs = ~time, lmer.df = 'asymptotic')

m2 = lme4::lmer(formula = f2, data = dat.fit)
emm2 = emmeans::emmeans(m2,specs = ~time, lmer.df = 'asymptotic')

m3 = lme4::lmer(formula = f3, data = dat.fit)
emm3 = emmeans::emmeans(m3, specs = ~time, lmer.df = 'asymptotic')

# contrast time differences 
clevels = levels(dat.fit$time)

#make custom contrasts 
c0 = c(1, 0, 0 ,0)
c1 = c(0, 1, 0 ,0)
c3 = c(0, 0, 1 ,0)
c4 = c(0, 0, 0 ,1)
contrast_list = list( "time1vs0" = c1 - c0,
                      'time22vs21' = c4 - c3)
clist = list('msig' = emm1, 'non.msig.sig' = emm2, 'sig' = emm3)

c.res = 
  lapply(clist, function(x) { 
    emmeans::contrast(object = x, method = contrast_list) %>% 
      broom::tidy()
  } ) %>% 
  bind_rows(.id = 'signature')
data.table::fwrite(x = c.res, file = paste0(datapath,'c.res.mDC.txt'), sep = '\t')  


# delta contrast 
cmat = emmeans::contrast(object = emm1, method = contrast_list)
pairs(cmat, reverse = TRUE)
# contrast              estimate     SE  df z.ratio p.value
# time22vs21 - time1vs0  0.00908 0.0172 Inf 0.528   0.5978 

# plotsingle cell distributionn and emmeans contrasts 
em_aes = list(theme_bw(), 
              coord_flip(), 
              theme(axis.title.y = element_text(size = 7),
                    axis.text.y = element_text(size = 7)), 
              scale_color_manual('grey')
)

plot.aes = list(theme_bw(), ylab(label ='Baseline high responder\nmDC signature'))
cu = sapply(c('grey', '#e2a359', 'grey', '#e2a359'), col.alpha, 0.8) %>% unname()
# combined signature change emm in p1 and change y value in p0
p0 = ggplot(dat.fit, aes(x = time, y = msig, fill = time)) + 
  geom_violin(show.legend = F) + 
  plot.aes + 
  xlab('time') + 
  scale_fill_manual(values = cu)+
  theme(axis.title.x = element_text(size = 12))
p0
ggsave(p0, filename = paste0(figpath, 'msig_mDC.pdf'), width = 2.5, height = 3)
p1 = plot(emm1) + em_aes
ggsave(p1, filename = paste0(figpath, 'msig_mDC.emm.pdf'), width = 1, height = 3)

# ifn -- change emm in p1 and change y value in p0
p0 = ggplot(dat.fit, aes(x = time, y = non.msig.sig, fill = time)) + 
  geom_violin(show.legend = F) + 
  plot.aes +
  xlab('time') + 
  scale_fill_manual(values = cu)+
  theme(axis.title.x = element_text(size = 12))
ggsave(p0, filename = paste0(figpath, 'non.msig.sig_mDC.pdf'), width = 2.5, height = 3)
p1 = plot(emm2) + em_aes
ggsave(p1, filename = paste0(figpath, 'non.msig.sig_mDC.emm.pdf'), width = 1, height = 3)


# non-ifn -- change emm in p1 and change y value in p0 
p0 = ggplot(dat.fit, aes(x = time, y = sig, fill = time)) + 
  geom_violin(show.legend = F) + 
  plot.aes +
  xlab('time') + 
  scale_fill_manual(values = cu)+
  theme(axis.title.x = element_text(size = 12))
ggsave(p0, filename = paste0(figpath, 'sig_mDC.pdf'), width = 2.5, height = 3)
p1 = plot(emm3) + em_aes
ggsave(p1, filename = paste0(figpath, 'sig_mDC.emm.pdf'), width = 1, height = 3)



