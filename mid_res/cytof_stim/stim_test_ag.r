suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(HDStIM))
suppressMessages(library(Rcpp))
suppressMessages(library(emmeans))
suppressMessages(library(lme4))
source('functions/MattPMutils.r')
# save paths 
figpath = here('mid_res/cytof_stim/figures/')
datapath = here('mid_res/cytof_stim/generated_data/')
dir.create(figpath); dir.create(datapath)

# baseline bulk theme 
cu1 = sapply(c('dodgerblue','red'), col.alpha, 0.2) %>% unname()
cu2 = c('dodgerblue', 'red')

# read data from stim cell selector 
d = readRDS(file = here('data/stim/mapped_data.rds'))

# celltype markers 
cmarkers = c('CD45',  'CD7',  'CD19',  'CD4', 'IgD', 'CD20', 
             "CD11c", "CD127", 'CD25', 'CD123', 'CD27', 'CD24', 
             'CD14', 'CD56', 'CD16', 'CD38', 'CD8', 'CD45RA',
             'CD3', 'HLA_DR')

# phenotyping markers 
pmarkers = c("pPLCg2", "pSTAT5", "AKT", "pSTAT1", "pP38",
             "pSTAT3", "IkBa", "pCREB", "pERK1_2", "pS6")

# innate cells 
innate.cells = c('CD14Mono', 'DC1', 'DC2')

########################
# make a data matrix 
dr = d$response_mapping_main
dr$cell_id = paste(dr$sample_id, rownames(dr),sep = '_')

dr = as.data.frame(dr) %>%  
  column_to_rownames('cell_id') %>%  
  mutate(sx = paste(patient_id, condition, stim_type, sep = '_')) %>% 
  mutate(response.stim = paste(condition, stim_type, sep = '_')) %>% 
  mutate(adjMFC = factor(condition,levels = c('high', 'low'))) 

# log transform the state marker matrix
mat = dr %>% select(all_of(pmarkers))
#mat = log(mat + 1 )

#make a separate dataframe of metadata 
prots = c(cmarkers, pmarkers)
met = dr %>% select(!all_of(prots))

# combine log transformed data back with md 
stopifnot(isTRUE(all.equal(rownames(met), rownames(mat))))
dr = cbind(met, mat )

# aggregate the protein markers to test the median log transformed marker intensity 
da = dr %>% 
  group_by(sx, response.stim, patient_id, batch, cell_population, stim_type) %>%  
  summarize_at(.vars = pmarkers, .funs = median) %>% 
  filter(cell_population %in% c('CD14Mono', 'DC1', 'DC2', 'CD16Mono'))


# model 
f1 = as.formula(
  prot ~ 0 + batch + response.stim + ( 1 | patient_id)
)

# separate by stim 
lps = da %>% 
  filter(stim_type %in% c('U', 'L')) %>% 
  filter(cell_population == 'CD14Mono')

# modify group factor for contrast estimand 
lps$response.stim = factor(
  lps$response.stim, 
  levels = c("high_U", "high_L", "low_U", "low_L")
)

dfit = lps 
contrast_sub = c ( '(high_L - high_U) - (low_L - low_U)' )
emmeans =  reslist =  list()
for (i in 1:length(pmarkers)) {
  
  # test 
  prot.name = pmarkers[i]
  print(prot.name)
  # # metadata 
  met = dfit %>% select(!all_of(pmarkers))
  
  # extract single protein 
  dat_vec = data.frame(dfit[ ,prot.name])
  colnames(dat_vec) = 'prot'
  # model data to fit 
  dat_fit = base::data.frame(cbind(dat_vec, met))
  
  # save quick plot 
  dplot = dat_fit %>% filter(response.stim %in% c('low_L', 'high_L'))
  dplot$response.stim = factor(dplot$response.stim, levels = c('low_L', 'high_L'))
  p = ggplot(dplot, aes(x = response.stim,y = prot , color = response.stim, fill = response.stim)) +
    theme_bw() + 
    geom_boxplot(show.legend = FALSE, outlier.shape = NA) + 
    geom_jitter(show.legend = FALSE, width = 0.2, shape = 21, color = 'black', stroke = 0.3) + 
    ylab(prot.name) + 
    scale_fill_manual(values = cu1) + 
    scale_color_manual(values = cu2) 
  p
  ggsave(p, filename = paste0(figpath,prot.name,'lps.pdf'), width = 1.75, height = 2)
  
  # fit models 
  m1 = tryCatch(
    lme4::lmer(formula = f1, data = dat_fit),
    error = function(e)
      return(NA)
  )
  
  # emmeans to apply contrast of fold changes 
  emm1 = tryCatch(
    emmeans::emmeans(
      object = m1, 
      specs = ~ 'response.stim', 
      data = dat_fit, 
      lmer.df = "asymptotic"),
    error = function(e) return(NA)
  )

  # apply contrasts 
  if (!is.na(emm1)) {
    m1.cont = contrast(emm1, method = 'revpairwise',adjust = NULL)
    m1cont = pairs(m1.cont,adjust = NULL) %>% 
      as.data.frame() %>% 
      filter( contrast == '(high_L - high_U) - (low_L - low_U)' )
    
    
    # store results 
    emmeans[[i]] = emm1 
    names(emmeans)[i] = prot.name
    
    # store contrast test 
    tidy.fit = m1cont %>%  
      mutate(prot = prot.name) %>%
      select(prot, everything())
    reslist[[i]] = tidy.fit
  } else{ 
    # store results 
    emmeans[[i]] = emm1 
    names(emmeans)[i] = prot.name
    reslist[[i]] = data.frame(
      prot = prot.name,
      contrast = NA,
      estimate = NA,
      SE = NA,
      df = NA,
      z.ratio = NA,
      p.value = NA
    )
    }
}
rd = do.call(rbind,reslist)
rd$stim = 'LPS'
data.table::fwrite(rd,file = paste0(datapath,'mono14_lps.txt'), sep = '\t')


######################
# PMA
 pma = da %>%
  filter(stim_type %in% c('U', 'P')) %>%
  filter(cell_population == 'CD14Mono')

# modify group factor for contrast estimand 
pma$response.stim = factor(
  pma$response.stim, 
  levels = c("high_U", "high_P", "low_U", "low_P")
)
dfit = pma 
contrast_sub = c ( '(high_P - high_U) - (low_P - low_U)' )
emmeans =  reslist =  list()
for (i in 1:length(pmarkers)) {
  #i = 5
  # test 
  prot.name = pmarkers[i]
  print(prot.name)
  # # metadata 
  met = dfit %>% select(!all_of(pmarkers))
  
  # extract single protein 
  dat_vec = data.frame(dfit[ ,prot.name])
  colnames(dat_vec) = 'prot'
  # model data to fit 
  dat_fit = base::data.frame(cbind(dat_vec, met))
  
  # save quick plot 
  dplot = dat_fit %>% filter(response.stim %in% c('low_P', 'high_P'))
  dplot$response.stim = factor(dplot$response.stim, levels = c('low_P', 'high_P'))
  p = ggplot(dplot, aes(x = response.stim,y = prot , color = response.stim, fill = response.stim)) +
    theme_bw() + 
    geom_boxplot(show.legend = FALSE, outlier.shape = NA) + 
    geom_jitter(show.legend = FALSE, width = 0.2, shape = 21, color = 'black', stroke = 0.3) + 
    ylab(prot.name) + 
    scale_fill_manual(values = cu1) + 
    scale_color_manual(values = cu2) 
  p
  ggsave(p, filename = paste0(figpath,prot.name,'pma.pdf'), width = 1.75, height = 2)
  
  # fit models 
  m1 = tryCatch(
    lme4::lmer(formula = f1, data = dat_fit),
    error = function(e)
      return(NA)
  )
  
  # emmeans to apply contrast of fold changes 
  emm1 = tryCatch(
    emmeans::emmeans(
      object = m1,
      specs = ~ 'response.stim',
      data = dat_fit,
      lmer.df = "asymptotic"
    ), 
    error = function(e) return(NA)
  )
  
  # apply contrasts 
  if (!is.na(emm1)) {
    m1.cont = contrast(emm1, method = 'revpairwise',adjust = NULL)
    m1cont = pairs(m1.cont,adjust = NULL) %>% 
      as.data.frame() %>% 
      filter( contrast == '(high_P - high_U) - (low_P - low_U)' )
    
    
    # store results 
    emmeans[[i]] = emm1 
    names(emmeans)[i] = prot.name
    
    # store contrast test 
    tidy.fit = m1cont %>%  
      mutate(prot = prot.name) %>%
      select(prot, everything())
    reslist[[i]] = tidy.fit
  } else{ 
    # store results 
    emmeans[[i]] = emm1 
    names(emmeans)[i] = prot.name
    reslist[[i]] = data.frame(
      prot = prot.name,
      contrast = NA,
      estimate = NA,
      SE = NA,
      df = NA,
      z.ratio = NA,
      p.value = NA
    )
  }
}
rd = do.call(rbind,reslist)
rd$stim = 'PMA'
data.table::fwrite(rd,file = paste0(datapath,'mono14_PMA.txt'), sep = '\t')

# IFN
# separate by stim 
ifn = da %>% 
  filter(stim_type %in% c('U', 'A')) %>% 
  filter(cell_population == 'CD14Mono')


# modify group factor for contrast estimand 
ifn$response.stim = factor(
  ifn$response.stim, 
  levels = c("high_U", "high_A", "low_U", "low_A")
)
dfit = ifn 
contrast_sub = c ( '(high_A - high_U) - (low_A - low_U)' )
emmeans =  reslist =  list()
for (i in 1:length(pmarkers)) {
  
  # test 
  prot.name = pmarkers[i]
  print(prot.name)
  # # metadata 
  met = dfit %>% select(!all_of(pmarkers))
  
  # extract single protein 
  dat_vec = data.frame(dfit[ ,prot.name])
  colnames(dat_vec) = 'prot'
  # model data to fit 
  dat_fit = base::data.frame(cbind(dat_vec, met))
  
  # save quick plot 
  dplot = dat_fit %>% filter(response.stim %in% c('low_A', 'high_A'))
  dplot$response.stim = factor(dplot$response.stim, levels = c('low_A', 'high_A'))
  p = ggplot(dplot, aes(x = response.stim,y = prot , color = response.stim, fill = response.stim)) +
    theme_bw() + 
    geom_boxplot(show.legend = FALSE, outlier.shape = NA) + 
    geom_jitter(show.legend = FALSE, width = 0.2, shape = 21, color = 'black', stroke = 0.3) + 
    ylab(prot.name) + 
    scale_fill_manual(values = cu1) + 
    scale_color_manual(values = cu2) 
  ggsave(p, filename = paste0(figpath,prot.name,'IFN.pdf'), width = 1.75, height = 2)
  
  # fit models
  m1 = tryCatch(
    lme4::lmer(formula = f1, data = dat_fit),
    error = function(e)
      return(NA)
  )
  
  # emmeans to apply contrast of fold changes
  emm1 = tryCatch(
    emmeans::emmeans(
      object = m1,
      specs = ~ 'response.stim',
      data = dat_fit,
      lmer.df = "asymptotic"
    ),
    error = function(e)
      return(NA)
  )
  
  # apply contrasts 
  if (!is.na(emm1)) {
    m1.cont = contrast(emm1, method = 'revpairwise',adjust = NULL)
    m1cont = pairs(m1.cont,adjust = NULL) %>% 
      as.data.frame() %>% 
      filter( contrast == '(high_A - high_U) - (low_A - low_U)' )
    
    
    # store results 
    emmeans[[i]] = emm1 
    names(emmeans)[i] = prot.name
    
    # store contrast test 
    tidy.fit = m1cont %>%  
      mutate(prot = prot.name) %>%
      select(prot, everything())
    reslist[[i]] = tidy.fit
  } else{ 
    # store results 
    emmeans[[i]] = emm1 
    names(emmeans)[i] = prot.name
    reslist[[i]] = data.frame(
      prot = prot.name,
      contrast = NA,
      estimate = NA,
      SE = NA,
      df = NA,
      z.ratio = NA,
      p.value = NA
    )
  }
}
rd = do.call(rbind,reslist)
rd$stim = 'IFN'
data.table::fwrite(rd,file = paste0(datapath,'mono14_IFN.txt'), sep = '\t')

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] lme4_1.1-26     Matrix_1.4-1    emmeans_1.5.4   Rcpp_1.0.9      HDStIM_0.1.0   
# [6] here_1.0.1      forcats_0.5.1   stringr_1.4.0   dplyr_1.0.4     purrr_0.3.4    
# [11] readr_1.4.0     tidyr_1.1.2     tibble_3.1.8    ggplot2_3.3.3   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
# [1] httr_1.4.2        jsonlite_1.7.2    splines_4.0.5     modelr_0.1.8      assertthat_0.2.1 
# [6] statmod_1.4.35    cellranger_1.1.0  yaml_2.2.1        pillar_1.8.1      backports_1.2.1  
# [11] lattice_0.20-41   glue_1.6.2        digest_0.6.27     rvest_0.3.6       minqa_1.2.4      
# [16] colorspace_2.0-0  sandwich_3.0-0    htmltools_0.5.2   plyr_1.8.6        pkgconfig_2.0.3  
# [21] broom_0.7.5       haven_2.4.3       xtable_1.8-4      mvtnorm_1.1-1     scales_1.1.1     
# [26] farver_2.0.3      generics_0.1.2    ellipsis_0.3.2    TH.data_1.0-10    withr_2.4.3      
# [31] Boruta_7.0.0      cli_3.4.1         survival_3.2-10   magrittr_2.0.3    crayon_1.4.1     
# [36] readxl_1.3.1      evaluate_0.15     estimability_1.3  fs_1.5.0          fansi_0.4.2      
# [41] nlme_3.1-152      MASS_7.3-53.1     xml2_1.3.2        rsconnect_0.8.25  tools_4.0.5      
# [46] data.table_1.14.0 hms_1.0.0         lifecycle_1.0.3   multcomp_1.4-16   munsell_0.5.0    
# [51] reprex_1.0.0      packrat_0.7.0     compiler_4.0.5    rlang_1.0.6       grid_4.0.5       
# [56] nloptr_1.2.2.2    ggridges_0.5.3    rstudioapi_0.13   rmarkdown_2.9     labeling_0.4.2   
# [61] boot_1.3-27       gtable_0.3.0      codetools_0.2-18  DBI_1.1.1         R6_2.5.0         
# [66] zoo_1.8-8         lubridate_1.8.0   knitr_1.39        fastmap_1.1.0     uwot_0.1.10      
# [71] utf8_1.2.2        rprojroot_2.0.2   stringi_1.5.3     parallel_4.0.5    vctrs_0.5.1      
# [76] xfun_0.30         dbplyr_2.1.0      tidyselect_1.2.0  coda_0.19-4  


