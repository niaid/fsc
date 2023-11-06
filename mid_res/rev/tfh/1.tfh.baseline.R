suppressMessages(library(tidyverse))
suppressMessages(library(here))

source(here('functions/MattPMutils.r'))

figpath = here('mid_res/rev/tfh/figures/'); dir.create(figpath)
datapath = here('mid_res/rev/tfh/generated_data/'); dir.create(datapath)

titer = read.csv(file = here('data/CHI_H1N1_data/titer/titer_processed.txt'),header = TRUE,sep = '\t')
titer$Subject = as.character(titer$Subject)

tester= read.csv(file = here('mid_res/rev/tfh/T4.exported.txt'),header= TRUE,sep = '\t')
colnames(tester)

# function provided by Mani 
parse.T4 <- function(day, controls=F, counts=F)
{
  hdr = read.delim(here("mid_res/rev/tfh/T4.exported.txt"), stringsAsFactors=F, header=F, nrows=1)
  a = read.delim(here("mid_res/rev/tfh/T4.exported.txt"), stringsAsFactors=F)
  stopifnot(grepl(".*T4  H1N1-", a$Sample), length(hdr)==ncol(a))
  if (!controls) {  a = a[!grepl("CHI-007",a$Sample),];  } else {  a$Sample = gsub("CHI-007","CHI007",a$Sample);  }
  
  # get the Tfh1 and Tfh1/Tfh17 cell popns, along with their CD278+PD1+ subsets.
  cidx = c(IDhlg = which(hdr == "Lymphocyte/single cells/Alive cells,Count"), #higher-level gate in the hierarchy
           ID242 = grep("Tfh1,Count$",hdr),
           ID244 = grep("Tfh1\\/.*CD278\\+.*PD\\-1\\+,Count$",hdr),
           ID247 = grep("Tfh1 Tfh17,Count$",hdr),
           ID249 = grep("Tfh1 Tfh17\\/.*CD278\\+.*PD\\-1\\+,Count$",hdr))
  stopifnot(length(cidx)==5)
  print(cbind(names(cidx), t(hdr[cidx])))
  stopifnot(make.names(hdr)[cidx] == colnames(a)[cidx])
  
  a$subject = do.call(rbind, strsplit(a$Sample, "-"))[,2]
  a$day = gsub("day ", "day", do.call(rbind, strsplit(a$Sample, "-"))[,3])
  print(table(a$day))
  if (!controls) {  a = a[a$day==day,];  } else {  a = a[a$subject=="CHI007",];  }
  print(table(a$day))
  print(table(a$subject))
  
  if (counts) {  a[,cidx['IDhlg']]=1;  }
  stopifnot(names(cidx)[1]=="IDhlg")
  b = matrix(NA, nrow(a), length(cidx)-1, dimnames=list(paste0("X",a$subject,ifelse(grepl("^day",a$day),"",paste0(".",a$day))), names(cidx)[-1]))
  for (i in names(cidx)[-1])
  {
    b[,i] = a[,cidx[i]]/a[,cidx['IDhlg']]
  }
  stopifnot(b >= 0, all(b <= 1) || counts)
  b = t(b)
  #write.table(file='tmp.T4.day0.txt', b, quote=F, sep="\t")
  b
}
d = parse.T4(day = 'day0') 
d2 = d %>% 
  as.data.frame() %>%  
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('subject') %>% 
  mutate(subject = as.character(str_sub(subject, 2,4))) %>% 
  filter(subject %in% titer$Subject) %>% 
  mutate(adjmfc = plyr::mapvalues(x = subject, from = titer$Subject,to = titer$adjMFC_class))

set.seed(1992)
d3 = d2 %>% 
  filter(!adjmfc=='NaN') %>%
  filter(adjmfc %in% c(0,2)) %>%  
  rename('AdjMFC' = adjmfc)
d3$AdjMFC[d3$AdjMFC == 0] <- 'low'
d3$AdjMFC[d3$AdjMFC == 2] <- 'high'
d3$AdjMFC= factor(d3$AdjMFC, levels = c('low','high'))
cu =c(col.alpha('dodgerblue',0.8),col.alpha('red', 0.8))

# TFH1 
p = ggplot(d3, aes(x = AdjMFC, y = ID242, fill = AdjMFC) )+
  theme_bw() + 
  geom_violin(show.legend = FALSE, trim = FALSE) + 
  geom_jitter(show.legend = FALSE, width = 0.3) + 
  ggpubr::stat_compare_means(method = 'wilcox',paired = FALSE) + 
  scale_fill_manual(values =cu) + 
  ylab('TFH percent of lymphocytes')

# TFHPd1 
p2 = ggplot(d3, aes(x = AdjMFC, y = ID244, fill = AdjMFC) )+
  theme_bw() + 
  geom_violin(show.legend = FALSE,trim = FALSE) + 
  geom_jitter(show.legend = FALSE, width = 0.3) + 
  ggpubr::stat_compare_means(method = 'wilcox',paired = FALSE) + 
  scale_fill_manual(values =cu) + 
  ylab('PD1+ TFH percent of lymphocytes')

# combine 
p3 = cowplot::plot_grid(p,p2, nrow = 1)

# save 
ggsave(p3,filename = paste0(figpath,'tfh.baseline.pdf'), width = 7, height = 3)

# save data 
saveRDS(object = d,file = paste0(datapath,'d.rds'))

