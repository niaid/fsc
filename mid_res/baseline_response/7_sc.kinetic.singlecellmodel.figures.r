# Early kinetics of baseline states 
# visualize model results
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))

# set paths 
datapath = here("mid_res/baseline_response/dataV3/")
figpath = here("mid_res/baseline_response/figuresV3/")
source('functions/MattPMutils.r')

mm2 = readRDS(file = here('mid_res/baseline_response/dataV3/mm2.rds')) %>% 
  filter(!module == 'null') %>% 
  filter(! singular_fit == 1)
mm2$cm = paste(mm2$celltype, mm2$module,sep = ' :: ')

# change to shorter names 
new.names = data.table::fread(file = here('mid_res/baseline_response/dataV3/baseline.module.name.shortened.txt'), sep = '\t')
new.names$cname2 = paste(new.names$celltype, new.names$shortname, sep = ' :: ')
new.names$cname1 = paste(new.names$celltype, new.names$module, sep = ' :: ')
mm2$cm = plyr::mapvalues(x = mm2$cm, from = new.names$cname1, to = new.names$cname2)

# assign to 'd'
d = mm2

# plot innate subset 
ds = d %>% filter(celltype  %in% c( 'CD14_Mono', 'CD16_Mono', 'mDC', "MAIT_Like" )) 

# filter the modules that did not have a effect in the single cell model 
m.rm = ds %>%
  filter(estimatetime0_group2vs1 > 0.1) %>%  
  mutate(m1s = estimatetime0_group2vs1 - std.errortime0_group2vs1) %>%
  filter(m1s > -0.01) %$% cm
saveRDS(m.rm, file = paste0(datapath, 'm.rm.rds'))
ds = ds %>% filter(cm %in% m.rm)

# remove cell label from module name to reduce clutter 
ds$cm = gsub(".*:","",ds$cm)
ds$celltype = str_replace_all(ds$celltype, pattern = '_',replacement =' ')


pl = list(
  # baseline 
  geom_point(shape = 23, size = 3, color = 'black', fill = col.alpha('red', 0.5)),
  geom_segment(aes(x = (estimatetime0_group2vs1 + -1*std.errortime0_group2vs1),
                   xend = estimatetime0_group2vs1 + 1*std.errortime0_group2vs1,
                   yend = cm),
               color = col.alpha('red', 0.5), 
               size = 2), 
  # day 1 
  geom_point(data = ds, aes(x = estimatetime1vs0, y = cm), size = 3, shape = 23, 
               color = 'black', fill = col.alpha('#e2a359', 0.5)),
  geom_segment(aes(x = (estimatetime1vs0 + -1*std.errortime1vs0),
                     xend = estimatetime1vs0 + 1*std.errortime1vs0,
                     yend = cm), 
                 color = col.alpha('#e2a359', 0.5),
                 size = 2), 
  theme(
    axis.text.y = element_text(color = 'black', size = 10),
    strip.background = element_blank(),
    axis.title.x = element_text(size = 14), 
    axis.text.x = element_text(size = 14), 
    strip.text = element_text(size = 12), 
    panel.spacing.x=unit(2, "lines")
    )
  )

# mdc, CD14 mono, CD16 mono, mait-like
p = 
  ggplot(ds, aes(x = estimatetime0_group2vs1, y = reorder(cm, estimatetime1vs0 ))) + 
  facet_grid(vars(celltype), scales = 'free', space = 'free') +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  pl + 
  scale_x_continuous( breaks= scales::pretty_breaks(), expand = c(0.15,0)) + 
  ylab('') + 
  xlab('contrast effect size ') 
p
ggsave(p, filename = paste0(figpath, 'mm2.innate.mait.pdf'), width = 5, height = 6.5)

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] magrittr_2.0.1  here_1.0.1      forcats_0.5.1   stringr_1.4.0   dplyr_1.0.4     purrr_0.3.4     readr_1.4.0    
# [8] tidyr_1.1.2     tibble_3.0.6    ggplot2_3.3.3   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.6        plyr_1.8.6        pillar_1.4.7      compiler_4.0.5    cellranger_1.1.0  dbplyr_2.1.0     
# [7] tools_4.0.5       digest_0.6.27     packrat_0.7.0     jsonlite_1.7.2    lubridate_1.7.9.2 lifecycle_1.0.0  
# [13] gtable_0.3.0      pkgconfig_2.0.3   rlang_0.4.10      reprex_1.0.0      cli_2.5.0         rstudioapi_0.13  
# [19] DBI_1.1.1         haven_2.3.1       withr_2.4.3       xml2_1.3.2        httr_1.4.2        fs_1.5.0         
# [25] generics_0.1.0    vctrs_0.3.6       hms_1.0.0         rprojroot_2.0.2   grid_4.0.5        tidyselect_1.1.0 
# [31] data.table_1.14.0 glue_1.4.2        R6_2.5.0          readxl_1.3.1      farver_2.0.3      modelr_0.1.8     
# [37] backports_1.2.1   scales_1.1.1      ellipsis_0.3.1    rvest_0.3.6       assertthat_0.2.1  colorspace_2.0-0 
# [43] labeling_0.4.2    stringi_1.5.3     munsell_0.5.0     broom_0.7.5       crayon_1.4.1     







