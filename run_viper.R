library(tidyverse)
library(OmnipathR)
library(openxlsx)
library(janitor)
library(reticulate)
library(ineq)
library(litteR)
library(aracne.networks)
library(binilib)
library(decoupleR)
library(textshape)


complex <- read.table('processed_complexes.txt', sep = '\t', header = T) #import complexes from CORUM database
complex_genes <- unique(unlist(str_split(complex$components_genesymbols, pattern = '_')))

cosmic <- read.table('cosmic_genes.tsv', header = T, sep = '\t')
cos <- cosmic %>% filter(Tier == 1) %>% pull(Gene.Symbol)

interest <- union(complex_genes, cos)

cancers <-  c('UCEC', 'LUAD', 'HNSC', 'RCC', 'OV', 'COAD', 'HCC') 
regulons <- list(UCEC = regulonucec, LUAD = regulonluad, HCC = regulonlihc, HNSC = regulonhnsc, RCC = regulonkirc, OV = regulonov, COAD = reguloncoad)

tidy_regulons <- function(regulon) {
  regulon %>%
    reg2tibble() %>% 
    mutate(source = any2symbol(source), target = any2symbol(target)) %>% #convert Entrez ID to gene symbol
    filter(!is.na(target), !is.na(source))  #cases where Entrez ID fails to match a symbol
}

regulons_proc <- lapply(regulons, tidy_regulons) %>% lapply(function(x) dplyr::filter(x, source %in% interest))
lapply(regulons_proc, pull, source) %>% lapply(unique) %>% lapply(NROW) 

#Fix HCC RNA
#HCC RNA dataset is formatted differently. This code removes the extra column of NAs at the end + 
#selects only tumour samples (not normals) + 
#converts Ensembl gene-version IDs to gene symbols
#deals with duplicated gene names by keeping the row with most non-zero values
ids <- read.xlsx('RNA/HCC_IDs.xlsx', startRow = 2)
ids$matcho <- paste0('X', ids[,3])
hcr = read.table('RNA/HCC_all.cct', header = T, sep = '\t') %>% select(protein, ids$matcho)
hcr$zeros <- rowSums(hcr == 0)
hcc_rna <- hcr %>% arrange(zeros) %>%
  mutate(protein = sub("\\..*","", protein)) %>%
  mutate(protein = any2symbol(protein)) %>%
  filter(!is.na(protein)) %>%
  select(-zeros) %>%
  distinct(protein, .keep_all = T)
write.table(hcc_rna, file = 'RNA/HCC_Tumor.cct', col.names = T, row.names = F, quote = F)


#Read in expression matrices
print('PROCESSING EXPRESSION MATRICES')
exp_mats <- vector(mode = 'list')
for (cancer in cancers) {
    exp_mats[[paste(cancer, 'Tumor', sep = '_')]] <- read.table(paste0('RNA/', cancer, '_Tumor.cct'), header = T, check.names = F) %>%
      column_to_rownames(loc = 1) %>%
      as.matrix() %>%
      na.omit()
  
}

#Run VIPER
viper_res <- mapply(run_viper, exp_mats, regulons_proc, verbose = TRUE, method = 'scale', pleiotropy = TRUE, SIMPLIFY = FALSE) #the right default colnames are already there, and pleiotropy is true as default
vip <- bind_rows(viper_res, .id = 'cancer_type') %>% mutate(cancer = gsub('_Tumor', '', cancer_type))
write.table(vip, file = 'viper_7_pleiotropy_scale_ofinterest.txt', sep = '\t', col.names = T, row.names = T, quote = F)

#Run ULM
ulm_res <- mapply(run_ulm, exp_mats, regulons_proc, center = T, SIMPLIFY = FALSE) #the right default colnames are already there, and pleiotropy is true as def>
ulm <- bind_rows(ulm_res, .id = 'cancer_type') %>% mutate(cancer = gsub('_Tumor', '', cancer_type))
write.table(ulm, file = 'ulm_7.txt', sep = '\t', col.names = T, row.names = T, quote = F)

#Run Weighted Mean
wmean_res <- mapply(run_wmean, exp_mats, regulons_proc,  SIMPLIFY = FALSE) #the right default colnames are already there, and pleiotropy is true as def>
wmean <- bind_rows(wmean_res, .id = 'cancer_type') %>% mutate(cancer = gsub('_Tumor', '', cancer_type))
write.table(wmean, file = 'wmean_7.txt', sep = '\t', col.names = T, row.names = T, quote = F)
