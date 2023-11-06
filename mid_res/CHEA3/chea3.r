library(tidyverse)
source(here('functions/scglmmr.functions.R'))
source(here('functions/MattPMutils.r'))

figpath = here('mid_res/CHEA3/figures/'); dir.create(figpath, recursive = TRUE)
datapath = here('mid_res/CHEA3/generated_data/'); dir.create(datapath, recursive = TRUE)

# write wrapper for chea3 
RunChea3 = function(gene.list){
  require(httr)
  require(jsonlite)
  get.res.1 <- function() {
    payload = list(query_name = names(gene.list), gene_set = gene.list[[1]])
    response = POST(url = "https://maayanlab.cloud/chea3/api/enrich/",
                    body = payload, 
                    encode = "json")
    return(response)
  }
  get.response = function() {
    httr::with_config(
      config = httr::config(ssl_verifypeer = FALSE),
      get.res.1()
    )
  }
  response2 = get.response()
  json = content(response2, "text")
  results = fromJSON(json)
  return(results)
}


g0 = readRDS(file = here('mid_res/baseline_response/dataV3/g0.sub.rds'))
g0 = lapply(g0, function(x) x %>%  filter(NES>0))
li = LeadingEdgeIndexed(gsea.result.list = g0,padj.threshold = 0.05)
li.joint = lapply(li, function(x) { x %>%  unlist %>% unique })
#mdc.genes = li$mDC %>% unlist %>%  unique 

# rm length 0 elements 
li.joint = Filter(li.joint, f = length)

# add other adjuvant signatures with various levels of pruning 
as03.sig.list = readRDS(file = here('mid_res/nat_adj/generated_data/V4/as03.sig.list.rds'))
mono.as03.sig.validated = readRDS(file = here('mid_res/nat_adj/generated_data/V4/mono.as03.sig.validated.rds'))
dc.as03.sig.validated = readRDS(file = here('mid_res/nat_adj/generated_data/V4/dc.as03.sig.validated.rds'))
li.joint = c(li.joint, as03.sig.list)
li.joint[['mono.as03.sig.validated']] = mono.as03.sig.validated
li.joint[['dc.as03.sig.validated']] = dc.as03.sig.validated

# natural adjuvant high responder baseline 
na.mono = read_csv(file = here('mid_res/nat_adj/generated_data/V4/mono.na.gsea.leadingEdge.txt'), col_names = FALSE)
na.dc = read_csv(file = here('mid_res/nat_adj/generated_data/V4/mdc.na.gsea.leadingEdge.txt'), col_names = FALSE)
li.joint[['Natural.Adjuvant.Mono']] = na.mono$X1
li.joint[['Natural.Adjuvant.DC']] = na.dc$X1
names(li.joint)


res = list()
for (i in 1:length(li.joint)) {
  gene.use = li.joint[i]
  res[[i]] = RunChea3(gene.list = li.joint[i])
}
names(res) = names(li.joint)
saveRDS(res,file = paste0(datapath, 'res.rds'))



# NA mono
test = 
  res$Natural.Adjuvant.Mono$`Integrated--meanRank` %>% 
  select(TF, Rank, Score) %>% 
  mutate(Rank = as.numeric(Rank),
         Score = as.numeric(Score)) %>% 
  mutate(enrscore = 1/Score)


test %>% head(21)
# TF Rank Score   enrscore
# 1         NFE4    1 10.00 0.10000000
# 2         TFEC    2 16.00 0.06250000
# 3         SPI1    3 16.50 0.06060606
# 4         MTF1    4 17.00 0.05882353
# 5        BATF2    5 22.67 0.04411116
# 6       ZNF467    6 26.00 0.03846154
# 7       ZNF438    7 27.67 0.03614022
# 8         MXD1    8 28.00 0.03571429
# 9        SP110    9 29.33 0.03409478
# 10        NFE2   10 31.00 0.03225806
# 11       STAT1   11 33.33 0.03000300
# 12       CEBPB   12 33.67 0.02970003
# 13      ZNF267   13 36.67 0.02727025
# 14         HLX   14 37.67 0.02654632
# 15        IRF7   15 38.67 0.02585984
# 16        TET2   16 40.50 0.02469136
# 17        IRF1   17 40.83 0.02449180
# 18      PLSCR1   18 43.00 0.02325581
# 19       TIGD3   19 43.33 0.02307870
# 20        IRF8   20 45.50 0.02197802
# 21       CEBPE   21 47.67 0.02097755


p=
  ggplot(test, aes(x = Rank, y = enrscore)) + 
  theme_minimal() +
  geom_point(size= 0.5, alpha = 1/5) + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  ggrepel::geom_text_repel(data = test %>% 
                             filter(enrscore > 0.02),
                           aes(x = Rank, y = enrscore, label = TF),
                           size = 2, max.overlaps = 50, nudge_x = 150, segment.size = 0.2) + 
  ylab('Predicted TF Enrichment Score') + 
  xlab('Predicted TF Enrichment Rank') + 
  ggtitle('CD14 Monocytes Natural adjuvant')
ggsave(p, filename = paste0(figpath,'tfrankmono.pdf'), width = 4 , height = 4)



# NA mDC
test = 
  res$Natural.Adjuvant.DC$`Integrated--meanRank` %>% 
  select(TF, Rank, Score) %>% 
  mutate(Rank = as.numeric(Rank),
         Score = as.numeric(Score)) %>% 
  mutate(enrscore = 1/Score)


p=
  ggplot(test, aes(x = Rank, y = enrscore)) + 
  theme_minimal() +
  geom_point(size= 0.5, alpha = 1/5) + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  ggrepel::geom_text_repel(data = test %>% 
                             filter(enrscore > 0.02),
                           aes(x = Rank, y = enrscore, label = TF),
                           size = 2, max.overlaps = 50, nudge_x = 150, segment.size = 0.2) + 
  ylab('Predicted TF Enrichment Score') + 
  xlab('Predicted TF Enrichment Rank') + 
  ggtitle('mDCs Natural adjuvant')
ggsave(p, filename = paste0(figpath,'tfrank.mdc.pdf'), width = 4 , height = 4)


