suppressMessages(library(tidyverse))
suppressMessages(library(here))

figpath =here("mid_res/ru/figures/"); dir.create(figpath)
datapath =here("mid_res/ru/generated_data/"); dir.create(datapath, recursive = TRUE)


# define adjuvant subjects 
met = data.table::fread(here("data/CHI_H5N1_data/clinical_info_adj.txt"))

# define adjuvant subjects 
adj.subjects = met %>%
  select(`Subject ID`, Adjuvant) %>%
  filter(Adjuvant == 'Adj') %>% 
  select(subjectid = `Subject ID`, Adjuvant) 

# load SPR data 
ru = read_delim(file = "data/CHI_H5N1_data/MN_abbinding/MPMEDIT_CHI_H5N1_AS03_SPR_data_2017_Khurana_SK.txt",delim = '\t')
CITE = c("H5N1-011", "H5N1-017", "H5N1-021", "H5N1-031", "H5N1-038", "H5N1-043")
ru.sub = 
  ru %>%
  separate(Sera,into = c('sx','time'),sep = '-') %>% 
  filter(time == ' D42') %>% 
  filter(subjectid %in% adj.subjects$subjectid) %>% 
  mutate(cite = ifelse(subjectid %in% CITE, '1', '0')) %>% 
saveRDS(ru.sub,file = paste0(datapath,'ru.sub.rds'))


p = 
  ggplot(data  = ru.sub, aes(x = Indo_RU_HA1, y = Viet_RU_HA1 )) + 
  theme_bw() +
  geom_smooth(method = "lm", color = "black") + 
  geom_point(shape = 21 , size = 3, fill = "deepskyblue3") + 
  theme(legend.position = "top") + 
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10, color = 'black')) +
  # ggrepel::geom_text_repel(data = ru.sub %>% filter( cite == 1), size = 2.8, aes(label = subjectid)) +
  xlab("Day 42 Antibody Binding (RU) \n Heterologous H5N1 strain (Vietnam)") + 
  ylab("Day 42 Antibody Binding (RU) \n vaccine H5N1 strain (Indonesia)") + 
  ggpubr::stat_cor(method = "pearson")
ggsave(p, filename = paste0(figpath,"RU_plot.pdf" ), width = 3.5, height = 3.5)

