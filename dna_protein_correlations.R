library(tidyverse)
library(data.table)
library(ggpubr)
library(patchwork)
library(broom)

rna <- fread('rna_cptac.txt', header = T)
prot <- fread('tumours_prot.txt', header = T) %>% filter(type_tumour != 'HCC')


### With Absolute CNs ###
cns_abs <- fread('cptac_gene_CNs_cut.bed')

sig_abs <- rna %>%
  left_join(cns_abs, by = c('Gene' = 'V1', 'Sample' = 'V2')) %>%
  left_join(prot, by = c('Gene', 'Sample')) %>%
  filter(!is.na(Abundance_tumour), !is.na(Expr), !is.na(V3), Expr != 0) %>%
  add_count(Gene, Cancer) %>%
  filter(n >= 30) 

cors_abs <- sig_abs %>%
  group_by(Gene, Cancer) %>%
  summarize(codr = cor(V3, Expr),
            corp = cor(Expr, Abundance_tumour),
            codp = cor(V3, Abundance_tumour)) 

p1 <- cors_abs %>%
  ggplot(aes(x = codr*corp, y = codp)) + geom_point() + stat_cor() + facet_wrap(~Cancer) +
  xlab('DNA-RNA * RNA-Protein') + ylab('DNA-Protein Correlation') + labs(subtitle = 'Absolute CNs')

p1_waterfall <- cors_abs %>%
  group_by(Gene) %>%
  summarize(med_codr = median(codr),
            med_corp = median(corp),
            med_codp = median(codp)) %>%
  mutate(DNA = 1, RNA = DNA*med_codr, Protein = RNA*med_corp) %>%
  sample_n(1000) %>%
  pivot_longer(DNA:Protein, names_to = 'Stage', values_to = 'Abundance') %>%
  ggplot(aes(x = factor(Stage, levels = c('DNA', 'RNA', 'Protein')), y = Abundance)) +
  geom_line(aes(group = Gene), alpha = 0.1) + geom_boxplot(width = 0.2) + xlab('') +
  labs(subtitle = 'Absolute CNs')


### With Relative CNs ###
cancers <- c('COAD', 'HNSC', 'LUAD', 'OV', 'RCC', 'UCEC')
lis <- vector(mode = 'list', length = length(cancers))
names(lis) <- cancers
for (cancer in cancers) {
  lis[[cancer]] <- fread(paste0('relcn/', cancer, '.cct'), check.names = F) %>% pivot_longer(-1, names_to = 'Sample', values_to = 'CN')
  names(lis[[cancer]])[1] <- 'Gene'
}
cns_rel <- bind_rows(lis, .id = 'Cancer')
cns_rel$Sample <- sub('\\.', '-', cns_rel$Sample)

sig_rel <- rna %>%
  left_join(cns_rel, by = c('Gene', 'Sample', 'Cancer')) %>%
  left_join(prot, by = c('Gene', 'Sample')) %>%
  filter(!is.na(Abundance_tumour), !is.na(Expr), !is.na(CN), Expr != 0)  %>%
  add_count(Gene, Cancer) %>%
  filter(n >= 30)

cors_rel <- sig_rel %>%
  group_by(Gene, Cancer) %>%
  summarize(codr = cor(CN, Expr),
            corp = cor(Expr, Abundance_tumour),
            codp = cor(CN, Abundance_tumour)) 

p2 <- cors_rel %>%
  ggplot(aes(x = codr*corp, y = codp)) + geom_point() + stat_cor() + facet_wrap(~Cancer) +
  xlab('DNA-RNA * RNA-Protein') + ylab('DNA-Protein Correlation') + labs(subtitle = 'Relative CNs')

p2_waterfall <- cors_rel %>%
  group_by(Gene) %>%
  summarize(med_codr = median(codr),
            med_corp = median(corp),
            med_codp = median(codp)) %>%
  mutate(DNA = 1, RNA = DNA*med_codr, Protein = RNA*med_corp) %>%
  sample_n(1000) %>%
  pivot_longer(DNA:Protein, names_to = 'Stage', values_to = 'Abundance') %>%
  ggplot(aes(x = factor(Stage, levels = c('DNA', 'RNA', 'Protein')), y = Abundance)) +
  geom_line(aes(group = Gene), alpha = 0.1) + geom_boxplot(width = 0.2) + xlab('') +
  labs(subtitle = 'Relative CNs')

### Linear Regression for each Gene ###

trep <- read.csv('~/Downloads/transcriptomic_reproducibility.csv') %>% rename(Gene = X)
repro <- read.csv('protein_repro.csv') %>% rename(Gene = X)
trepro <- repro %>% left_join(trep, by = c('Gene'))

sig_abs_lm <- sig_abs %>%
  left_join(trepro, by = c('Gene')) %>%
  filter(!is.na(Aggregated.Reproducibility.Rank), !is.na(CCLE.Klijn)) %>%
  group_by(Gene) %>%
  mutate(cancers = n_distinct(Cancer)) %>%
  filter(cancers > 1)

# reg_abs <- vector(mode = 'list', length = NROW(unique(sig_abs_lm$Gene)))
# names(reg_abs) <- unique(sig_abs_lm$Gene)
# for (gene in unique(sig_abs_lm$Gene)) {
#    print(gene)
#    reg_abs[[gene]] <- sig_abs_lm %>%
#      lm(Abundance_tumour ~ V3 + CCLE.Klijn + Aggregated.Reproducibility.Rank + Cancer, subset = Gene == gene, data = .) %>%
#      summary() %>% 
#      tidy() 
# }

genes <- read.table('genes_hg38.tsv', header = T, sep = '\t')
complex <- read.table('corum_complexes.txt', sep = '\t', header = T)
complex_genes <- complex %>% separate_rows(components_genesymbols, sep = '_') %>% pull(components_genesymbols) 
genes$complex <- ifelse(genes$external_gene_name %in% complex_genes, 'Complex', 'Non-Complex')

lr_rp_abs <- cors_abs %>% left_join(trepro, by = c('Gene')) %>%
  filter(!is.na(Aggregated.Reproducibility.Rank), !is.na(CCLE.Klijn)) %>%
  group_by(Gene) %>%
  mutate(cancers = n_distinct(Cancer)) %>%
  filter(cancers > 1) %>%
  left_join(genes, by = c('Gene' = 'external_gene_name')) %>%
  filter(!is.na(copyconserved)) %>%
  lm(corp ~ DBO + CNV + Ohnolog + copyconserved + complex + Aggregated.Reproducibility.Rank + Cancer, data = .) %>%
  summary() %>% tidy()

lr_dp_abs <- cors_abs %>% left_join(trepro, by = c('Gene')) %>%
  filter(!is.na(Aggregated.Reproducibility.Rank), !is.na(CCLE.Klijn)) %>%
  group_by(Gene) %>%
  mutate(cancers = n_distinct(Cancer)) %>%
  filter(cancers > 1) %>%
  left_join(genes, by = c('Gene' = 'external_gene_name')) %>%
  filter(!is.na(copyconserved)) %>%
  lm(codp ~ DBO + CNV + Ohnolog + copyconserved + complex + CCLE.Klijn + Aggregated.Reproducibility.Rank, data = .) %>%
  summary() %>% tidy()

lr_rp_rel <- cors_rel %>% left_join(trepro, by = c('Gene')) %>%
  filter(!is.na(Aggregated.Reproducibility.Rank), !is.na(CCLE.Klijn)) %>%
  group_by(Gene) %>%
  mutate(cancers = n_distinct(Cancer)) %>%
  filter(cancers > 1) %>%
  left_join(genes, by = c('Gene' = 'external_gene_name')) %>%
  filter(!is.na(copyconserved)) %>%
  lm(corp ~ DBO + CNV + Ohnolog + copyconserved + complex + Aggregated.Reproducibility.Rank + Cancer, data = .) %>%
  summary() %>% tidy()

lr_dp_rel <- cors_rel %>% left_join(trepro, by = c('Gene')) %>%
  filter(!is.na(Aggregated.Reproducibility.Rank), !is.na(CCLE.Klijn)) %>%
  group_by(Gene) %>%
  mutate(cancers = n_distinct(Cancer)) %>%
  filter(cancers > 1) %>%
  left_join(genes, by = c('Gene' = 'external_gene_name')) %>%
  filter(!is.na(copyconserved)) %>%
  lm(codp ~ DBO + CNV + Ohnolog + copyconserved + complex + CCLE.Klijn + Aggregated.Reproducibility.Rank, data = .) %>%
  summary() %>% tidy()


reg_plot <- function(df, co, cns) {
df %>%
  mutate(significant = ifelse(p.value < 0.05/4, 'Yes', 'No')) %>%
  filter(term != '(Intercept)') %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate,
                         ymin = estimate - std.error, ymax = estimate + std.error,
                         colour = significant)) +
  geom_point() +
  geom_errorbar() +
  geom_hline(colour = 'black', linetype = 'dashed', yintercept = 0) +
  xlab('') + ylab('Regression Estimate') +
  coord_flip() +
  theme(legend.position = 'bottom') +
  labs(subtitle = cns) +
  ggtitle(co)

}

r1 <- reg_plot(lr_rp_abs, 'RNA-Protein Correlation', 'Absolute CNs')
r2 <- reg_plot(lr_dp_abs, 'DNA-Protein Correlation', 'Absolute CNs')
r3 <- reg_plot(lr_rp_rel, 'RNA-Protein Correlation', 'Relative CNs')
r4 <- reg_plot(lr_dp_rel, 'DNA-Protein Correlation', 'Relative CNs')

cag <- cors_abs %>% left_join(genes, by = c('Gene' = 'external_gene_name')) %>% filter(!is.na(DBO))
crg <- cors_rel %>% left_join(genes, by = c('Gene' = 'external_gene_name')) %>% filter(!is.na(DBO))

dsg_rp_abs <- cag %>% ggplot(aes(x = complex, y = corp)) + geom_violin(aes(fill = complex)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('RNA-Protein Correlation') | 
    cag %>% ggplot(aes(x = DBO, y = corp)) + geom_violin(aes(fill = DBO)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('RNA-Protein Correlation') | 
    cag %>% ggplot(aes(x = copyconserved, y = corp)) + geom_violin(aes(fill = copyconserved)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('RNA-Protein Correlation') | 
    cag %>% filter(CNV %in% c('Class P', 'Class B')) %>% ggplot(aes(x = CNV, y = corp)) + geom_violin(aes(fill = CNV)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('RNA-Protein Correlation')  

dsg_dp_abs <- cag %>% ggplot(aes(x = complex, y = codp)) + geom_violin(aes(fill = complex)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('DNA-Protein Correlation') | 
  cag %>% ggplot(aes(x = DBO, y = codp)) + geom_violin(aes(fill = DBO)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('DNA-Protein Correlation') | 
  cag %>% ggplot(aes(x = copyconserved, y = codp)) + geom_violin(aes(fill = copyconserved)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('DNA-Protein Correlation') | 
  cag %>% filter(CNV %in% c('Class P', 'Class B')) %>% ggplot(aes(x = CNV, y = codp)) + geom_violin(aes(fill = CNV)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('DNA-Protein Correlation')  

dsg_rp_rel <- crg %>% ggplot(aes(x = complex, y = corp)) + geom_violin(aes(fill = complex)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('RNA-Protein Correlation') | 
  crg %>% ggplot(aes(x = DBO, y = corp)) + geom_violin(aes(fill = DBO)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('RNA-Protein Correlation') | 
  crg %>% ggplot(aes(x = copyconserved, y = corp)) + geom_violin(aes(fill = copyconserved)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('RNA-Protein Correlation') | 
  crg %>% filter(CNV %in% c('Class P', 'Class B')) %>% ggplot(aes(x = CNV, y = corp)) + geom_violin(aes(fill = CNV)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('RNA-Protein Correlation')  

dsg_dp_rel <- crg %>% ggplot(aes(x = complex, y = codp)) + geom_violin(aes(fill = complex)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('DNA-Protein Correlation') | 
  crg %>% ggplot(aes(x = DBO, y = codp)) + geom_violin(aes(fill = DBO)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('DNA-Protein Correlation') | 
  crg %>% ggplot(aes(x = copyconserved, y = codp)) + geom_violin(aes(fill = copyconserved)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('DNA-Protein Correlation') | 
  crg %>% filter(CNV %in% c('Class P', 'Class B')) %>% ggplot(aes(x = CNV, y = codp)) + geom_violin(aes(fill = CNV)) + geom_boxplot(alpha = 0.2) + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('DNA-Protein Correlation')  

