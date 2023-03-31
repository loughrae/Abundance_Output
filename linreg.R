library(tidyverse)
library(textshape)
library(broom)
library(openxlsx)
library(TCGAutils)
library(ggrepel) 

OV_cn <- read.table('CN/OV.cct', header = T, check.names = F)
HNSC_cn <- read.table('CN/HNSC.cct', header = T, check.names = F)
RCC_cn <- read.table('CN/RCC.cct', header = T, check.names = F)
LUAD_cn <- read.table('CN/LUAD.cct', header = T, check.names = F)
LUAD_cn <- LUAD_cn %>% column_to_rownames(loc = 1) %>% t() %>% as.data.frame() %>% rownames_to_column(var = 'gene') #reformat LUAD CN to match the others
UCEC_cn <- read.table('CN/UCEC.cct', header = T, check.names = F)
COAD_cn <- read.table('CN/COAD.cct', header = T, check.names = F)

cn_list <- list(UCEC = UCEC_cn, LUAD = LUAD_cn, HNSC = HNSC_cn, RCC = RCC_cn, OV = OV_cn, COAD = COAD_cn)

long_cn <- lapply(cn_list, pivot_longer, cols = -1, names_to = 'Sample', values_to = 'CN') %>% lapply(setNames,  c('Gene', 'Sample', 'CN')) %>% bind_rows(.id = 'Cancer')
aneuploidy <- long_cn %>% group_by(Sample, Cancer) %>% summarize(std = sd(CN, na.rm = T))

complex <- read.table('processed_complexes.txt', sep = '\t', header = T) #import complexes from CORUM database
deviations <- read.table('deviations_cptac.txt', sep = '\t', header = T)

ids <- read.xlsx('RNA/HCC_IDs.xlsx', startRow = 2)
colnames(ids) <- c('Prot_T', 'Prot_N', 'RNA_T', 'RNA_N')
ids$matcho <- paste0('X', ids$RNA_T)

vip <- read.table('viper_6_pleiotropy_scale_ofinterest.txt', sep = '\t', header = T) %>%
  left_join(ids, by = c('condition' = 'matcho')) %>% 
  mutate(condition = case_when(cancer == 'HCC' ~ Prot_T, TRUE ~ condition))

tog <- complex %>%
  separate_rows(components_genesymbols) %>%
  left_join(vip, by = c('components_genesymbols' = 'source')) %>%
  filter(!is.na(score)) %>%
  add_count(name, condition) %>% 
  filter(n >= 3) %>% 
  group_by(name, condition) %>%
  summarize(avg_score = median(score)) %>%
  left_join(deviations, by = c('name', 'condition' = 'Sample')) %>%
  filter(!is.na(cv)) %>% 
  left_join(aneuploidy, by = c('condition' = 'Sample')) %>%
  filter(!is.na(std)) %>%
  group_by(name) %>%
  mutate(types = n_distinct(Cancer)) %>%
  filter(types > 1)

head(tog)
NROW(unique(tog$name)) 

complex_names <- as.list(unique(tog$name))
names(complex_names) <- complex_names
length(complex_names)

regs <- lapply(complex_names, function(complex) lm(avg_score ~ cv + std + Cancer, data = tog[tog$name == complex,]))
regs_df <- regs %>% lapply(summary) %>% lapply(tidy) %>% bind_rows(.id = 'Complex') 
regs_df %>% write.table('complex_regs.txt', sep = '\t', quote = F, row.names = F, col.names = T)

regs_df %>% filter(term == 'cv')  %>% ggplot(aes(x = reorder(Complex, estimate), y = estimate, colour = p.value < 0.05 / nrow(regs_df))) + 
geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), colour = 'black') + 
geom_point()  + 
labs(colour = 'Significant') + 
ylab('Regression Estimate') + xlab('Complex') + 
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
#geom_text_repel(aes(label = substr(Complex, start = 1, stop = 10)))

ggsave(filename = 'complex_regression.png', width = 20, height = 15)

#### Gene-Level Analyses ####
d <- read.table('../cn_cptac/segments.txt', header = T)
names(d)[ncol(d)] <- 'filename'
d <- d %>% filter(Copy_Number != 'Copy_Number') %>% 
mutate(fl = gsub("\\..*","",filename)) %>% 
separate(fl, into = c('folder', 'file'), sep = '_') %>%
mutate(Copy_Number = as.numeric(Copy_Number), Start = as.numeric(Start), End = as.numeric(End)) %>% 
mutate(seg_len = End - Start) %>% 
mutate(prod = Copy_Number * seg_len) %>%
group_by(GDC_Aliquot, folder) %>%
summarize(ploidy = sum(prod) / sum(seg_len)) 

case_ids <- UUIDtoUUID(d$folder, to = 'case_id')
sample_names <- UUIDtoBarcode(unique(case_ids$cases.case_id), from_type = 'case_id')
d2 <- d %>% left_join(case_ids, by = c('folder' = 'file_id')) %>% left_join(sample_names, by = c('cases.case_id' = 'case_id'))
write.table(d2, 'cptac_ploidy.txt', sep = '\t', col.names = T, row.names = F, quote = F)

prot <- read.table('tumours_prot.txt', header = T, sep = '\t') 

gene_scores <- vip %>% left_join(d2, by = c('condition' = 'submitter_id')) %>%
	 filter(!is.na(ploidy)) %>%
	 left_join(prot, by = c('source' = 'Gene', 'condition' = 'Sample')) %>%
	 filter(!is.na(unlog)) %>% 
	 group_by(source) %>% 
	 mutate(types = n_distinct(cancer_type)) %>%
	 filter(types > 1)

gene_names <- as.list(unique(gene_scores$source))
names(gene_names) <- gene_names #name the elements of the list to be used for .id later
length(gene_names)

gene_regs <- lapply(gene_names, function(gene) lm(score ~ unlog + I(unlog^2) + ploidy + unlog*ploidy + cancer_type, data = gene_scores[gene_scores$source == gene,]))
gene_regs_df <- gene_regs %>% lapply(summary) %>% lapply(tidy) %>% bind_rows(.id = 'Gene')
gene_regs_df %>% write.table('gene_regs.txt', sep = '\t', quote = F, row.names = F, col.names = T)

ng <- length(gene_names)
gene_regs_df %>% filter(term == 'unlog')  %>%
mutate(significant = ifelse(p.value < 0.05 / ng, 'Yes', 'No')) %>%
mutate(lab = ifelse(significant == 'Yes', Gene, '')) %>%
ggplot(aes(x = reorder(Gene, estimate), y = estimate, colour = p.value < 0.05 / ng)) +
geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), colour = 'black') +
geom_point()  +
labs(colour = 'Significant') +
ylab('Regression Estimate') + xlab('Gene') +
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
#geom_text_repel(aes(label = lab))

ggsave(filename = 'gene_regression.png', width = 20, height = 15)
