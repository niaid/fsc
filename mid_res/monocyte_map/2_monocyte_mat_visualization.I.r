# R 3.5 
# make visualization of main monocyte pseudotime axis. 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
suppressMessages(library(here))
suppressMessages(library(monocle))
source("functions/analysis_functions.R")
source('functions/MattPMutils.r')
btm = readRDS("signature_curation/BTM_li.rds")
figpath = here("mid_res/monocyte_map/figures/"); dir.create(figpath, recursive = TRUE)
sm = readRDS(file = here("mid_res/monocyte_map/generated_data/sm_cd14_cd16_d1_monocle_object.rds"))

######
# visualization 
p = plot_cell_trajectory(sm, color_by = "Pseudotime")
df = ggplot_build(p)[["plot"]][["data"]] %>% 
  rename(component_1 = data_dim_1  ,component_2 = data_dim_2)
df$adjmfc.time = factor(df$adjmfc.time, levels = c("d0 low", "d1 low", "d0 high", "d1 high"))
library(cowplot) 

# set plot theme 
theme.set = list(theme_bw(), 
                 theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(), 
                       axis.ticks.y = element_blank(),
                       axis.ticks.x = element_blank()))
# reverse order to make the root node go left to right 
df$component_1 = -1* df$component_1

# create main plot
p1 = ggplot(df, aes(x = component_1, y = component_2, fill = Pseudotime)) +
  theme.set + 
  geom_point(shape = 21, size = 2.2, color = "grey" ,stroke = 0.1) +
  theme(legend.position = c(0.1, 0.7)) + 
  theme(strip.background = element_blank()) + 
  ylab("mRNA trajectory component 2") + 
  xlab("mRNA trajectory component 1") + 
  scale_fill_viridis_c(option = "B") + 
  scale_color_manual(values = c("grey", "black"))
ggsave(p1, filename = paste0(figpath, "mono_trajectory_only.png"), width = 6, height = 5)

# make a background plot on which to add the canvas marginal plots
pnull = ggplot(df, aes(x = component_1, y = component_2, fill = Pseudotime)) +
  theme.set + 
  theme(legend.position = c(0.1, 0.7)) + 
  theme(strip.background = element_blank()) + 
  theme(axis.title.x = element_text(size = 17)) + 
  theme(axis.title.y = element_text(size = 17)) +
  ylab("mRNA trajectory component 2") + 
  xlab("mRNA trajectory component 1") + 
  scale_fill_viridis_c(option = "plasma") + 
  scale_color_manual(values = c("grey", "black"))

# real time for top margin 
time.col = c( 
  col.alpha(acol = 'black', 0.1), 
  col.alpha(acol = ggsci::pal_jama()(2)[2], 0.4) 
  ) 
timecol2 = c(
  'black',
  ggsci::pal_jama()(2)[2] 
)
xd = axis_canvas(pnull, axis = "x") +
  geom_density(data = df, aes(x = component_1, color = timepoint), size = 1) + 
  scale_color_manual(values = time.col)

## test 
xd = axis_canvas(pnull, axis = "x") +
  geom_density(data = df, aes(x = component_1, fill = timepoint, color = timepoint), size = 1) + 
  scale_fill_manual(values = time.col) + 
  scale_color_manual(values = timecol2)
xd
## 


# CD14 vs CD16 protein for bottom margin 
x16 = axis_canvas(pnull, axis = "x") +
  geom_smooth(data = df, 
              aes(x = component_1, y = CD16_PROT), 
              method = "loess", se = TRUE,
              color = 'black') + 
  geom_smooth(data = df,
              aes(x = component_1, y = CD14_PROT),
              method = "loess", se = TRUE, 
              color = 'grey') 
# add to plot   
p2 <- insert_xaxis_grob(pnull, xd, grid::unit(.2, "null"), position = "top")
p4 =  insert_xaxis_grob(p2, x16, grid::unit(.3, "null"), position = "bottom")
p6 = ggdraw(p4)
p6
# save 
ggsave(p6, filename = paste0(figpath, "monocyte_merged_plot.pdf"), width = 6, height = 8)
ggsave(p6, filename = paste0(figpath, "monocyte_merged_plot.png"), width = 6, height = 8)






