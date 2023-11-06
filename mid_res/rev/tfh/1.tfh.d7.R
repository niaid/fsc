suppressMessages(library(tidyverse))
suppressMessages(library(here))
source(here('functions/MattPMutils.r'))

figpath = here('mid_res/rev/tfh/figures_7/'); dir.create(figpath)
datapath = here('mid_res/rev/tfh/generated_data_7/'); dir.create(datapath)

titer = read.csv(file = here('data/CHI_H1N1_data/titer/titer_processed.txt'),header = TRUE,sep = '\t')
titer$Subject = as.character(titer$Subject)

tester= read.csv(file = here('mid_res/rev/tfh/T4.exported.txt'),header= TRUE,sep = '\t')
colnames(tester)

# function provided by Mani. 
# modify in steps below to load multi timepoints

# parse.T4 <- function(day, controls=F, counts=F)
# {
#   hdr = read.delim(here("mid_res/rev/tfh/T4.exported.txt"), stringsAsFactors=F, header=F, nrows=1)
#   a = read.delim(here("mid_res/rev/tfh/T4.exported.txt"), stringsAsFactors=F)
#   stopifnot(grepl(".*T4  H1N1-", a$Sample), length(hdr)==ncol(a))
#   if (!controls) {  a = a[!grepl("CHI-007",a$Sample),];  } else {  a$Sample = gsub("CHI-007","CHI007",a$Sample);  }
#   
#   # get the Tfh1 and Tfh1/Tfh17 cell popns, along with their CD278+PD1+ subsets.
#   cidx = c(IDhlg = which(hdr == "Lymphocyte/single cells/Alive cells,Count"), #higher-level gate in the hierarchy
#            ID242 = grep("Tfh1,Count$",hdr),
#            ID244 = grep("Tfh1\\/.*CD278\\+.*PD\\-1\\+,Count$",hdr),
#            ID247 = grep("Tfh1 Tfh17,Count$",hdr),
#            ID249 = grep("Tfh1 Tfh17\\/.*CD278\\+.*PD\\-1\\+,Count$",hdr))
#   stopifnot(length(cidx)==5)
#   print(cbind(names(cidx), t(hdr[cidx])))
#   stopifnot(make.names(hdr)[cidx] == colnames(a)[cidx])
#   
#   a$subject = do.call(rbind, strsplit(a$Sample, "-"))[,2]
#   a$day = gsub("day ", "day", do.call(rbind, strsplit(a$Sample, "-"))[,3])
#   print(table(a$day))
#   if (!controls) {  a = a[a$day==day,];  } else {  a = a[a$subject=="CHI007",];  }
#   print(table(a$day))
#   print(table(a$subject))
#   
#   if (counts) {  a[,cidx['IDhlg']]=1;  }
#   stopifnot(names(cidx)[1]=="IDhlg")
#   b = matrix(NA, nrow(a), length(cidx)-1, 
#              dimnames= list(
#                paste0("X",a$subject,ifelse(grepl("^day",a$day),"",paste0(".",a$day))),
#                names(cidx)[-1]
#                )
#              )
#   for (i in names(cidx)[-1])
#   {
#     b[,i] = a[,cidx[i]]/a[,cidx['IDhlg']]
#   }
#   stopifnot(b >= 0, all(b <= 1) || counts)
#   b = t(b)
#   #write.table(file='tmp.T4.day0.txt', b, quote=F, sep="\t")
#   b
# }


# modify to look at day 0 and 7 to calc fc
#parse.T4 <- function(day, controls=F, counts=F){
hdr = read.delim(
  here("mid_res/rev/tfh/T4.exported.txt"),
  stringsAsFactors = F,
  header = F,
  nrows = 1
)
a = read.delim(here("mid_res/rev/tfh/T4.exported.txt"), stringsAsFactors = FALSE)
stopifnot(grepl(".*T4  H1N1-", a$Sample), length(hdr) == ncol(a))
# if (!controls) {
a = a[!grepl("CHI-007", a$Sample), ]
# } else {
# a$Sample = gsub("CHI-007", "CHI007", a$Sample)
# }

# get the Tfh1 and Tfh1/Tfh17 cell popns, along with their CD278+PD1+ subsets.
cidx = c(
  IDhlg = which(hdr == "Lymphocyte/single cells/Alive cells,Count"),
  #higher-level gate in the hierarchy
  ID242 = grep("Tfh1,Count$", hdr),
  ID244 = grep("Tfh1\\/.*CD278\\+.*PD\\-1\\+,Count$", hdr),
  ID247 = grep("Tfh1 Tfh17,Count$", hdr),
  ID249 = grep("Tfh1 Tfh17\\/.*CD278\\+.*PD\\-1\\+,Count$", hdr)
)
stopifnot(length(cidx) == 5)
print(cbind(names(cidx), t(hdr[cidx])))
stopifnot(make.names(hdr)[cidx] == colnames(a)[cidx])

a$subject = do.call(rbind, strsplit(a$Sample, "-"))[, 2]
a$day = gsub("day ", "day", do.call(rbind, strsplit(a$Sample, "-"))[, 3])
print(table(a$day))

#if (!controls) {
#  a = a[a$day == day, ]
a = a[a$day %in% c('day0', 'day7'), ]
# } else {
#   a = a[a$subject == "CHI007", ]
# }
print(table(a$day))
print(table(a$subject))

# if (counts) {
#   a[, cidx['IDhlg']] = 1
# }
stopifnot(names(cidx)[1] == "IDhlg")

# Append days onto sample names 
b = matrix(NA, nrow(a), length(cidx)-1, 
           dimnames= list(
             #paste0("X",a$subject,ifelse(grepl("^day",a$day),"",paste0(".",a$day))),
             paste0(a$subject,a$day),
             names(cidx)[-1]
           )
)

# calc fc matrix as in original function 
for (i in names(cidx)[-1]){
  b[, i] = a[, cidx[i]] / a[, cidx['IDhlg']]
}
stopifnot(b >= 0, all(b <= 1) || counts)

# make dataframe 
bd = as.data.frame(b)
bd = bd %>%  rownames_to_column('sample')
fc = bd %>% 
  pivot_longer(!sample, names_to = 'population', values_to = 'freq') %>%  
  mutate(day = str_sub(sample, -4,-1)) %>% 
  arrange(population) %>% 
  mutate(fold_change = lead(freq) / freq) %>% # non log FC Data. 
  filter(day == 'day0') %>%  
  select(-day) %>%  
  mutate(subject = str_sub(sample,1,3)) %>%  
  filter(subject %in% titer$Subject) %>% 
  mutate(AdjMFC = plyr::mapvalues(x = subject, from = titer$Subject,to = titer$adjMFC_class)) 

fc$AdjMFC[fc$AdjMFC == '0'] <- 'low'
fc$AdjMFC[fc$AdjMFC == '2'] <- 'high'
fc$AdjMFC = factor(fc$AdjMFC, levels = c('low', 'high'))

# 
set.seed(1992)
cu = c(col.alpha('dodgerblue',0.8),col.alpha('red', 0.8))


p = ggplot(fc %>%  filter(population == 'ID242') %>%  filter(AdjMFC %in% c('low','high')), aes(x = AdjMFC, y = fold_change, fill = AdjMFC) )+
  theme_bw() + 
  geom_violin(show.legend = FALSE,trim = FALSE) + 
  geom_jitter(show.legend = FALSE, width = 0.3) + 
  ggpubr::stat_compare_means(method = 'wilcox',paired = FALSE) + 
  scale_fill_manual(values =cu) + 
  ylab('TFH day 7 vs baseline fold change')

p2 = ggplot(fc %>%  filter(population == 'ID244') %>%  filter(AdjMFC %in% c('low','high')), aes(x = AdjMFC, y = fold_change, fill = AdjMFC) )+
  theme_bw() + 
  geom_violin(show.legend = FALSE,trim = FALSE) + 
  geom_jitter(show.legend = FALSE, width = 0.3) + 
  ggpubr::stat_compare_means(method = 'wilcox',paired = FALSE) + 
  scale_fill_manual(values =cu) + 
  ylab('PD1+ TFH day 7 vs baseline fold change')

# combine 
p3 = cowplot::plot_grid(p,p2, nrow = 1)
ggsave(p3,filename = paste0(figpath,'tfh.d7fc.pdf'), width = 7, height = 3)

# what is the actual relative frquency of the cells 
apply(bd[ ,2:5], 1, median) %>% median
#  0.001618353   # 0.1% across samples both timepoints 

# load adjuvant signatures 
# d7fc vs nat adj sig 
mono.validated = readRDS(file = here('mid_res/nat_adj/generated_data/V4/mono.as03.sig.validated.rds'))
dc.validated = readRDS(file = here('mid_res/nat_adj/generated_data/V4/dc.as03.sig.validated.rds'))
adjuvant.signatures = c(list('Mono_AS03_validated' = mono.validated), list('DC_AS03_validated' = dc.validated))
# load average baseline high  and low responders from unadjuvanted cohort 
av0 = readRDS(file = here('mid_res/nat_adj/generated_data/V4/av0.rds'))
mono.sig.av2 = av0$CD14_Mono %>% 
  filter(gene %in% adjuvant.signatures$Mono_AS03_validated) %>% 
  group_by(sample, response) %>% 
  summarize(meansig = mean(count)) %>% 
  mutate(subject = str_sub(sample, 1,3)) 

# combine with tFH data 
d4 = full_join(mono.sig.av2, fc, by = 'subject')
d4 = d4 %>% filter(response %in% c('Low', "High"))

# tfh
p1 = ggplot(d4 %>%  filter(population == 'ID242'), aes(x = fold_change , y = meansig)) + 
  theme_bw() +
  geom_point() +
  stat_smooth(method = 'lm', color = 'black', show.legend = FALSE) + 
  ggpubr::stat_cor(show.legend = FALSE) +
  scale_color_manual(values = cu) + 
  xlab('TFH cell day 7 fold change') + 
  ylab('Baseline CD14 Monocyte \n  Natural Adjuvant signature') + 
  ggtitle('high outlier included')
p1
ggsave(p1,filename = paste0(figpath,'tfh.d7fc.NATADJ_sig_combined.pdf'),width = 3, height = 3)

# outlier exc
p1 = ggplot(d4 %>%  filter(population == 'ID242')%>%  filter(!subject == '209'), aes(x = fold_change , y = meansig)) + 
  theme_bw() +
  geom_point() +
  stat_smooth(method = 'lm', color = 'black', show.legend = FALSE) + 
  ggpubr::stat_cor(show.legend = FALSE) +
  scale_color_manual(values = cu) + 
  xlab('TFH cell day 7 fold change') + 
  ylab('Baseline CD14 Monocyte \n  Natural Adjuvant signature') + 
  ggtitle('high outlier excluded')
p1
ggsave(p1,filename = paste0(figpath,'tfh.d7fc.NATADJ_sig_combined_outlierexc.pdf'),width = 3, height = 3)



# pd1 +
p1 = ggplot(d4 %>%  filter(population == 'ID244'), aes(x = fold_change , y = meansig)) + 
  theme_bw() +
  geom_point() +
  stat_smooth(method = 'lm', color = 'black', show.legend = FALSE) + 
  ggpubr::stat_cor(show.legend = FALSE) +
  scale_color_manual(values = cu) + 
  xlab('PD1+ TFH cell day 7 fold change') + 
  ylab('Baseline CD14 Monocyte \n  Natural Adjuvant signature') + 
  ggtitle('high outlier included')
p1
ggsave(p1,filename = paste0(figpath,'tfh.d7fc.NATADJ_sig_combined_PD1.pdf'),width = 3, height = 3)

# tfh
p1 = ggplot(d4 %>%  filter(population == 'ID242'), aes(x = fold_change , y = meansig, color = response)) + 
  theme_bw() +
  geom_point() +
  facet_wrap(~response) + 
  stat_smooth(method = 'lm', color = 'black', show.legend = FALSE) + 
  ggpubr::stat_cor(show.legend = FALSE) +
  scale_color_manual(values = cu) + 
  xlab('TFH cell day 7 fold change') + 
  ylab('Baseline CD14 Monocyte \n  Natural Adjuvant signature') + 
  ggtitle('high responder outlier included')
p1
ggsave(p1,filename = paste0(figpath,'tfh.d7fc.NATADJ_sig_thfmain.pdf'),width = 7, height = 4)


p1 = ggplot(d4 %>%  filter(population == 'ID242'), aes(x = fold_change , y = meansig, color = response)) + 
  theme_bw() +
  geom_point() +
  stat_smooth(method = 'lm', color = 'black', show.legend = FALSE) + 
  ggpubr::stat_cor(aes(x = fold_change , y = meansig), inherit.aes = FALSE, show.legend = FALSE) +
  scale_color_manual(values = cu) + 
  xlab('TFH cell day 7 fold change') + 
  ylab('Baseline CD14 Monocyte \n  Natural Adjuvant signature')  
p1
ggsave(p1,filename = paste0(figpath,'tfh.d7fc.NATADJ_sig_thfmain.pdf'), width = 5, height = 4)

## take away the high value outlier and view trend
p1 = ggplot(d4 %>%  filter(population == 'ID242') %>%  filter(!subject == '209'), aes(x = fold_change , y = meansig, color = response)) + 
  theme_bw() +
  geom_point() +
  #facet_wrap(~response) + 
  stat_smooth(method = 'lm', color = 'black', show.legend = FALSE) + 
  ggpubr::stat_cor(show.legend = FALSE) +
  scale_color_manual(values = cu) + 
  xlab('TFH cell day 7 fold change') + 
  ylab('Baseline CD14 Monocyte \n  Natural Adjuvant signature')  +
  ggtitle('high responder outlier excluded')
p1
ggsave(p1,filename = paste0(figpath,'tfh.d7fc.NATADJ_sig_thfmain_outlierRM.pdf'), width = 7, height = 4)


mdc= av0$mDC %>% 
  filter(gene %in% adjuvant.signatures$DC_AS03_validated) %>% 
  group_by(sample, response) %>% 
  summarize(meansig = mean(count)) %>% 
  mutate(subject = str_sub(sample, 1,3)) 

# combine with tFH data 
d4 = full_join(mdc, fc, by = 'subject')
d4 = d4 %>% filter(response %in% c('Low', "High"))

# tfh
p1 = ggplot(d4 %>%  filter(population == 'ID242'), aes(x = fold_change , y = meansig)) + 
  theme_bw() +
  geom_point() +
  stat_smooth(method = 'lm', color = 'black', show.legend = FALSE) + 
  ggpubr::stat_cor(show.legend = FALSE) +
  scale_color_manual(values = cu) + 
  xlab('TFH cell day 7 fold change') + 
  ylab('Baseline CD14 mDC \n  Natural Adjuvant signature') + 
  ggtitle('high outlier included')
p1
ggsave(p1,filename = paste0(figpath,'tfh.d7fc.NATADJ_sig_MDC_combined.pdf'),width = 3, height = 3)

####################################

