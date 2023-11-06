# 4.2 
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(simr))
source(here('functions/MattPMutils.r'))
source(here('functions/scglmmr_functions/model_result_interaction.r'))


# save paths
figpath = here('mid_res/rev/rev.figs/');dir.create(figpath)
datapath = here('mid_res/rev/rev.data/');dir.create(datapath)
s = readRDS(file = here('data/h1h5_annotated_with_meta.rds'))
s@meta.data %>% 
  group_by(sample, celltype_joint) %>% 
  tally() %>%  
  filter(celltype_joint %in% c( "CD14_Mono", "mDC", "BC_Naive")) %>%
  group_by(celltype_joint) %>% 
  summarize_at(.vars = 'n',.funs = median)
  median(n)
  

  # A tibble: 3 × 2
  # celltype_joint     n
  # <chr>          <dbl>
  # 1 BC_Naive       206. 
  # 2 CD14_Mono      448. 
  # 3 mDC             29.5
rm(s); gc()

# load contrast samplemd 
samplemd = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/samplemd12.rds'))
pb = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/pb12.rds'))

# designmat for edgeR normalization 
met = samplemd[ ,c('gender', 'scaledage', 'time.group')]
mat = model.matrix( ~ 0 + time.group +  gender + scaledage, data = met)
betas = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/fit12e.rds'))
betas$CD14_Mono$coefficients
deltas = ExtractResult(model.fit.list = betas,what = 'lmer.z.ranks',coefficient.number = 1, 'delta')


# rank by absolute effect size 
deltas = lapply(deltas, function(x) x %>% abs() %>% sort(decreasing = TRUE))

# pb data 
power.calc = list()
names(pb)
pb = pb[c("CD14_Mono", "mDC", "BC_Naive")]
deltas = deltas[c( "CD14_Mono", "mDC", "BC_Naive")]

# specify contrast 
c00 = c(1,0,0,0)
c01 = c(0,1,0,0)
c10 = c(0,0,1,0)
c11 = c(0,0,0,1)
contrast_list = list(
  "time1vs0_group2vs1" = (c01 - c00) - (c11 - c10) 
)

# write a custom function for simr to extract FC delta p values from marginal means contrast 
test_delta <- function(m1, ...) {
  # m1 to be the lme4 fit per celltype & gene 
  contrast_fit = emmeans::emmeans(m1, specs = ~time.group) %>% 
    emmeans::contrast(method = contrast_list) %>% 
    summary()
  contrast_fit$p.value
}
attr(test_delta, "text") <- function(...) NULL

power.calc = list()
fit.list = list()
for (i in 1:length(pb)) {

  ## process data 
  d= edgeR::DGEList(counts = pb[[i]], samples = samplemd)
  gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, design = mat)
  print(names(pb)[i]);print(table(gtable))
  d = d[gtable, keep.lib.sizes=FALSE]
  d = edgeR::calcNormFactors(object = d)
  cpm = edgeR::cpm(d, log = TRUE)
  
  
  #select 200 genes 
  top.genes = names(deltas[[i]])[1:200] # can tune
  cpm.top = cpm[top.genes, ]
  

  # contrast mixed model 
  f1 <- gene ~ 0 + time.group + gender + scaledage + (1|subjectid) 
  
  fit.list = list()
  for (u in 1:length(top.genes)) {
    
    # fit a single gene 
    gene.name = top.genes[u]
    dat = data.frame(gene = cpm[top.genes[u],], samplemd)
    m1 = lme4::lmer(formula = f1, data = dat)
    
    # simulate data based on model 
    ps = powerSim(fit = m1, test = test_delta, seed = 1990, nsim = 50)
    
    # custom function output is not captured; grab it 
    ps.res = capture.output(ps,type = 'output')
    
    # reformat and reformat again outside 
    pow = gsub(" ", "", ps.res[2])
    
    # get model stuff 
    random_effect_var = insight::get_variance_random(m1, tolerance = 0)
    residual_var = insight::get_variance_residual(m1)
    
    
    # store results 
    fit.list[[u]] = data.frame(
      celltype = names(pb)[i], 
      gene = gene.name, 
      lmer.z.statistic = deltas[[i]][[u]], 
      power.sim = pow, 
      random_effect_var, 
      residual_var
      )
    
  }
  power.calc[[i]] = bind_rows(fit.list)
}
# names(power.calc) = names(pb)
power.data = bind_rows(power.calc)
saveRDS(power.data,file = paste0(datapath, 'power.data.rds'))
power.data= readRDS(file = here('mid_res/rev/rev.data/power.data.rds'))

# reformat annoying powersim captured output from the custom function. 
extract.power.ci = function(x) {
  # Remove "%" 
  x <- gsub("%", "", x)
  
  # Extract numbers between "(" and ","
  left_paren <- gregexpr("\\(", x)[[1]]
  comma <- gregexpr(",", x)[[1]]
  power.lower <- as.numeric(substr(x, left_paren + 1, comma - 1))
  
  # Extract numbers between "," and ")"
  right_paren <- gregexpr("\\)", x)[[1]]
  power.upper <- as.numeric(substr(x, comma + 1, right_paren - 1))
  
  # Extract numbers before "%"
  x <- gsub("%", "", power.data$power.sim)
  power = as.numeric(sub("\\(.*", "", x))
  
  return(data.frame(power, power.lower, power.upper))
}

# Applying the function to the dataframe
power.df <- cbind(power.data, extract.power.ci(power.data$power.sim))
saveRDS(power.df,file = paste0(datapath, 'power.df.rds'))
power.df = readRDS(file = here('mid_res/rev/rev.data/power.df.rds'))

# plot of power analysis 
p = ggplot(power.df, aes(x = reorder(gene, power), y = power, color = celltype))+ 
  theme_minimal() + 
  xlab('gene') + 
  geom_point(size = 1, show.legend = FALSE) +
  geom_errorbar(aes(ymax = power.upper, ymin = power.lower), lwd = 0.1, show.legend = FALSE) + 
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) + 
  facet_wrap(~celltype) + 
  ggsci::scale_color_d3()

ggsave(p,filename = paste0(figpath, 'power.analysis.contrastmodel.pdf'), width = 8, height = 4)
power.df %>% filter(power > 70) %>% group_by(celltype) %>% tally()
# # A tibble: 3 × 2
# 1 BC_Naive     77
# 2 CD14_Mono   155
# 3 mDC          72

p =
ggplot(power.df, aes(x = lmer.z.statistic, y = power)) + 
  theme_bw() + 
  geom_errorbar(aes(ymax = power.upper, ymin = power.lower),
                lwd = 0.2, alpha = 0.5, show.legend = FALSE) + 
  geom_point(size = 1, show.legend = FALSE, alpha = 0.4) +
  facet_wrap(~celltype, scales = 'free') + 
  ggsci::scale_color_d3() + 
  xlab('effect size') 
p
ggsave(p,filename = paste0(figpath, 'power.analysis.contrastmodel.effectsize.pdf'), width = 10, height = 4)

p = 
ggplot(power.df, aes(x = power)) + 
  theme_bw() + 
  facet_wrap(~celltype, scales = 'free') + 
  geom_histogram(fill = 'deepskyblue3', color = 'black', bins = 12) 
ggsave(p,filename = paste0(figpath, 'power.analysis.contrastmodel.distribution.pdf'), width = 10, height = 4)  


p =
  ggplot(power.df, aes(x = random_effect_var, y = power)) + 
  theme_minimal() + 
  geom_point(size = 1, show.legend = FALSE) +
  facet_wrap(~celltype) + 
  ggsci::scale_color_d3()
ggsave(p,filename = paste0(figpath, 'power.analysis.contrastmodel.ranef.pdf'), width = 8, height = 4)


p =
  ggplot(power.df, aes(x = residual_var, y = power)) + 
  theme_minimal() + 
  geom_point(size = 1, show.legend = FALSE) +
  facet_wrap(~celltype) + 
  ggsci::scale_color_d3()
ggsave(p,filename = paste0(figpath, 'power.analysis.contrastmodel.resid.pdf'), width = 8, height = 4)
