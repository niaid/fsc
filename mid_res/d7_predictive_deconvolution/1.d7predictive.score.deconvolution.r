# Day 7 signature deconvolution direct with module score.  
# R 3.5.1 
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(magrittr))
suppressMessages(library(here))
source("functions/analysis_functions.R")
source('functions/MattPMutils.r')

# save path 
figpath = here("mid_res/d7_predictive_deconvolution/figures/")
dir.create(figpath)

# day 7 core signatures. 
sig7 = readRDS("signature_curation/core_d7.rds")

# h1 data baseline and day 7 cells. 
h1 = ReadCohort(joint_object_dir = "data/h1h5_annotated_with_meta.rds", cohort = "H1N1")
h1 = SetAllIdent(h1, id = "time_cohort") %>% 
  SubsetData(ident.use = "d7")

# add module score 
h1 = AddModuleScore(h1, genes.list = sig7, seed.use = 1, enrich.name = names(sig7))

# get long for mfor visualization of module score distribtion. 
df_sig = h1@meta.data %>% select(
  celltype_joint,
  timepoint,
  LI.M156_Plasma_Cell = `LI.M156.0 plasma cell b cell Ig6`,
  CHI_4 = 'CHI 4 d710',
  CHI_5 = `CHI 5 d711`,
  CHI_d7_Response = `CHI d7 Response9`
) %>%
  mutate(celltype_b = if_else(
    celltype_joint %in% c("BC_Mem", "BC_Naive",  "CD38_Bcell", "pDC"),
    true = celltype_joint,
    false = "other"
  )) %>%
  gather(module, module_score, LI.M156_Plasma_Cell:CHI_d7_Response) %>%
  mutate(timepoint = factor(timepoint, levels = c("d0", "d7"))) %>%
  mutate(celltype_b = str_replace_all(
    string = celltype_b,
    pattern = "_",
    replacement = " "
  )) %>%
  mutate(celltype_joint = str_replace_all(
    string = celltype_joint,
    pattern = "_",
    replacement = " "
  )) %>%
  mutate(celltype_b = factor(
    celltype_b,
    levels = c("CD38 Bcell", "pDC", "BC Naive", "BC Mem", "other")
  ))

# plot CHI predictive sig from array 
subplot = df_sig %>% filter(timepoint == "d7" &
                              module %in% c("CHI_d7_Response", "CHI_4", "LI.M156_Plasma_Cell"))
subplot = subplot %>% mutate(
  module = dplyr::recode(
    module,
    "LI.M156_Plasma_Cell" = "M156",
    "CHI_d7_Response" = "Antibody sig",
    "CHI_4" = "CHI 4"
  )
)

grey1 = col.alpha(acol = 'grey',alpha = 0.2)

p = ggplot(subplot,   aes(
  x = reorder(celltype_joint, module_score),
  y = module_score,
  fill = celltype_joint
)) +
  theme_bw() +
  geom_violin(
    show.legend = FALSE,
    scale = "width", 
    size = 0.2, 
    draw_quantiles = c(0.5)
  ) +
  facet_wrap( ~ module,
              as.table = TRUE,
              scales = "free_x",
              ncol = 4) +
  geom_hline(yintercept  = 0,
             size = 0.5,
             linetype = "dashed") +
  coord_flip() +
  scale_fill_manual(values = c(rep(grey1, 5), "red", rep(grey1, 17)))  +
  theme(
    panel.spacing = unit(0.2, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 7, color = "black"),
    panel.border = element_blank()
  ) +
  xlab("") + ylab("single cell score distribution \n day 7 bulk predictive signatures") +
  ggtitle("Day 7 post-vaccination")
ggsave(p, filename = paste0(figpath,"d7core_Tirosh_mod_score_celltypes_sub.pdf"), width = 5, height = 5)

