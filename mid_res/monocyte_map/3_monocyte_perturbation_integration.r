# R 3.5 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
suppressMessages(library(here))
library(monocle)
source("functions/analysis_functions.R")
source("functions/MattPMutils.r")
figpath = here("mid_res/monocyte_map/figures/"); dir.create(figpath, recursive = TRUE)
datapath = here("mid_res/monocyte_map/generated_data/"); dir.create(datapath, recursive = TRUE)

# load monocle object 
sm = readRDS(file = here("mid_res/monocyte_map/generated_data/sm_cd14_cd16_d1_monocle_object.rds"))

# load day 1 enrichment from monocytes 
g1c = readRDS(file = here('mid_res/1_H1N1_pseudobulk_DE/dataV4/gsea/g1c.rds'))
g1 = g1c$CD14_Mono %>% 
  filter(NES > 0) %>% 
  filter(padj < 0.05)
monole = g1$leadingEdge
names(monole) = g1$pathway
cgene2 = unique(unlist(monole))
saveRDS(cgene2,file = paste0(datapath, 'cgene2.rds'))

# define branch dependent genes
de.branch <- BEAM(sm[cgene2, ], branch_point = 1, cores = 4)
saveRDS(de.branch,file = paste0(datapath,'de_branch.rds'))
de.branch = readRDS(file = here('mid_res/monocyte_map/generated_data/de_branch.rds'))
de.branch.sub = de.branch %>% filter(qval < 0.05)
branch.genes = as.character(de.branch.sub$gene_short_name)


# Visualization of branch genes 
led = exprs(sm)[branch.genes, ] %>% 
  t() %>% 
  as.matrix() %>% 
  as.data.frame()
d = cbind(led, sm@phenoData@data)

dd = d %>% select(OASL, CCL2, IFITM2, FCER1G, TNFSF10, FCGR1B)
dd$Pseudotime = d$Pseudotime
dd$timepoint = d$timepoint

# mono act 
dd = dd %>%  filter(Pseudotime > 5)
dd = dd %>% gather(gene, value, OASL:FCGR1B)

# vis theme 
mtheme =  list(
  theme_bw(), 
  geom_smooth(size = 2), 
  theme(axis.title = element_text(size = 18)), 
  ggsci::scale_color_jama(alpha = 0.8),
  theme(legend.position = c(0.2, 0.8)),
  xlab('Pseudotime') 
)

# example category 1 gene 
p1 = ggplot(dd %>% filter(gene %in% c('CCL2')), 
            aes(x = Pseudotime, y = value,  color = timepoint)) + 
  ylab('CCL2 Expression') +
  mtheme 
ggsave(p1, filename = paste0(figpath,'Category1_1_CCL2.pdf'), width = 3.5, height = 3.5)

# example category 2 gene 
p1 = ggplot(dd %>% filter(gene %in% c('TNFSF10')), 
            aes(x = Pseudotime, y = value,  color = timepoint)) + 
  ylab('TNFSF10 Expression') +
  mtheme 
ggsave(p1, filename = paste0(figpath,'Category1_2_TNFSF10.pdf'), width = 3.5, height = 3.5)

# example category 2 gene
p1 = ggplot(dd %>% filter(gene %in% c('FCER1G')), 
            aes(x = Pseudotime, y = value,  color = timepoint)) + 
  ylab('FCER1G Expression') +
  mtheme 
ggsave(p1, filename = paste0(figpath,'Category2_1_FCER1G.pdf'), width = 3.5, height = 3.5)

# example category 2 gene
p1 = ggplot(dd %>% filter(gene %in% c('IFITM2')), 
            aes(x = Pseudotime, y = value,  color = timepoint)) + 
  mtheme + 
  ylab('IFITM2 Expression')
ggsave(p1, filename = paste0(figpath,'Category2__1_IFITM2.pdf'), width = 3.5, height = 3.5)



#########################
# leading edge signatures analysis 

# Interferon 
sig = intersect(branch.genes, monole$`reactome interferon signaling`)
pdf(file = paste0(figpath, 'IFN_branch_de.pdf'),width = 5, height = 5)
plot_genes_branched_heatmap(sm[sig, ],
                            branch_point = 1,
                            cores = 1,
                            num_clusters = 3,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()

# set categ 
cat2 = c('IFITM2', 'PTPN1', 'EIF4E2', 'IFITM3', 'HLA-C')
cat1 = sig[!sig %in% cat2]
# get data 
led = exprs(sm)[sig, ] %>% 
  t() %>% 
  as.matrix() %>% 
  as.data.frame()
d = cbind(led, sm@phenoData@data)
dd = d %>% select(sig)
dd = apply( dd, 2, scale.simple) %>%  as.data.frame()
dd$Pseudotime = d$Pseudotime
dd$timepoint = d$timepoint
index1 = sig[1]
index2 = sig[length(sig)]
d3 = dd %>% gather(gene, value, index1:index2 ) 
d3$cat = ifelse(d3$gene %in% cat1, '1', '2')



p1 = ggplot(data = d3 %>%  
              filter(Pseudotime > 5 & cat ==1 )) + 
  geom_smooth(data = d3 %>% 
                filter(Pseudotime > 5 & timepoint == 'd0' & cat ==1 ), 
              mapping = aes(x  = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(1), se = FALSE) + 
  geom_smooth(data = d3 %>% filter(Pseudotime > 5 & timepoint == 'd1' & cat ==1 ),
              mapping =  aes(x  = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(2)[2],
              se = FALSE, size = 2) + 
  theme_bw() + 
  geom_vline(xintercept = 9.5, linetype = 'dashed') + 
  theme(strip.background = element_blank()) + 
  ylab('reactome interferon \n Category 1 genes ') + 
  theme(axis.title = element_text(size = 18))
p1
ggsave(p1,filename = paste0(figpath, 'IFNcat1.pdf'), width = 3.7, height = 3.5)


p2 = ggplot(data = d3 %>%  
              filter(Pseudotime > 5 & cat ==2 )) + 
  geom_smooth(data = d3 %>% 
                filter(Pseudotime > 5 & timepoint == 'd0' & cat ==2 ), 
              mapping = aes(x  = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(1), se = FALSE) + 
  geom_smooth(data = d3 %>% 
                filter(Pseudotime > 5 & timepoint == 'd1' & cat ==2),
              mapping =  aes(x  = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(2)[2], 
              se = FALSE, size = 2) + 
  theme_bw() + 
  geom_vline(xintercept = 9.5, linetype = 'dashed') + 
  theme(strip.background = element_blank()) + 
  ylab('reactome interferon \n Category 2 genes ') + 
  theme(axis.title = element_text(size = 18))
p2
ggsave(p2,filename = paste0(figpath, 'IFNcat2.pdf'), width = 3.7, height = 3.5)


# mtor hypoxia 
sig = intersect(branch.genes,
                c(monole$`HALLMARK hypoxia`, monole$`HALLMARK MTORC1 signaling`)
                )

pdf(file = paste0(figpath, 'mtorhypoxia_branch_de.pdf'),width = 5, height = 5)
plot_genes_branched_heatmap(sm[sig, ],
                            branch_point = 1,
                            cores = 1,
                            num_clusters = 3,
                            use_gene_short_name = T,
                            show_rownames = T)

dev.off()

#set categ 
cat2 = c('CTSC', 'PFKL', 'ACTR3', 'CITED2', 'PGK1', 'INSIG1', 'CHST2')
cat1 = sig[!sig %in% cat2]

# get data 
led = exprs(sm)[sig, ] %>% 
  t() %>% 
  as.matrix() %>% 
  as.data.frame()
d = cbind(led, sm@phenoData@data)
dd = d %>% select(sig)
dd = apply( dd, 2, scale.simple) %>%  as.data.frame()
dd$Pseudotime = d$Pseudotime
dd$timepoint = d$timepoint
index1 = sig[1]
index2 = sig[length(sig)]
d3 = dd %>% gather(gene, value, index1:index2 ) 
d3$cat = ifelse(d3$gene %in% cat1, '1', '2')


p1 = ggplot(data = d3 %>% 
              filter(Pseudotime > 5 & cat ==1 )) + 
  geom_smooth(data = d3 %>% 
                filter(Pseudotime > 5 & timepoint == 'd0' & cat ==1 ), 
              mapping = aes(x  = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(1), se = FALSE) + 
  geom_smooth(data = d3 %>% 
                filter(Pseudotime > 5 & timepoint == 'd1' & cat ==1 ),
              mapping =  aes(x = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(2)[2], 
              se = FALSE, size = 2) + 
  theme_bw() + 
  geom_vline(xintercept = 9.5, linetype = 'dashed') + 
  theme(strip.background = element_blank()) + 
  ylab('MTORC1 and Hypoxia\n Category 1 genes ') + 
  theme(axis.title = element_text(size = 18))
p1
ggsave(p1,filename = paste0(figpath, 'mtorcat1.pdf'), width = 3.7, height = 3.5)


p2 = ggplot(data = d3 %>% 
              filter(Pseudotime > 5 & cat ==2 )) + 
  geom_smooth(data = d3 %>% 
                filter(Pseudotime > 5 & timepoint == 'd0' & cat ==2 ), 
              mapping = aes(x  = Pseudotime, y = value, group = gene),
              color = ggsci::pal_jama(alpha = 0.2)(1), se = FALSE) + 
  geom_smooth(data = d3 %>% filter(Pseudotime > 5 & timepoint == 'd1' & cat ==2),
              mapping =  aes(x  = Pseudotime, y = value, group = gene), 
              color = ggsci::pal_jama(alpha = 0.2)(2)[2], 
              se = FALSE, size = 2) + 
  theme_bw() + 
  geom_vline(xintercept = 9.5, linetype = 'dashed') + 
  theme(strip.background = element_blank()) + 
  ylab('MTORC1 and Hypoxia\n Category 2 genes ') + 
  theme(axis.title = element_text(size = 18))
ggsave(p2,filename = paste0(figpath, 'mtorcat2.pdf'), width = 3.7, height = 3.5)



