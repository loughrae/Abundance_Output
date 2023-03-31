library(tidyverse)
library(OmnipathR)
library(openxlsx)
library(janitor)
library(reticulate)
library(ineq)
library(litteR)
library(aracne.networks)
library(textshape)
library(binilib)
library(decoupleR)

#### COMPLEXES ####
print('PROCESSING COMPLEXES')
complex <- read.table('corum_complexes.txt', sep = '\t', header = T) %>%
  mutate(name_separations = str_count(components_genesymbols, '_')) %>% 
  mutate(n_members = name_separations + 1) %>%
  filter(name_separations > 1) %>%
  mutate(zero_before = str_count(stoichiometry, '0:')) %>%
  mutate(zero_after = str_count(stoichiometry, ':0')) %>% 
  filter(zero_before + zero_after == 0) %>% #ensure stoichiometry fully known
  filter(name != 'Molybdopterin synthase') %>% 
  arrange(desc(n_members)) %>% 
  distinct(name, .keep_all = T) #a few complexes have same name and slightly different makeup; keep the one with most members

write.table(complex, file = 'processed_complexes.txt', sep = '\t', col.names = T, row.names = F, quote = F)


#### PROTEOMICS ——— CALCULATING STOICHIOMETRIC IMBALANCE ####
#read in CPTAC tumour files
print('READING IN CPTAC')
coad <- read.table('proteomics/COAD_tumor.cct', sep = '\t', header = T, check.names = F) %>%
  pivot_longer(cols = -attrib_name, names_to = 'Sample', values_to = 'Abundance') %>%
  mutate(type = 'COAD') %>%
  mutate(Gene = attrib_name) %>%
  select(-attrib_name)

ov <- read.table('proteomics/ovarian_tumor.cct', sep = '\t', header = T, check.names = F) %>%
  pivot_longer(cols = -gene, names_to = 'Sample', values_to = 'Abundance') %>%
  mutate(type = 'OV') %>%
  mutate(Gene = gene) %>%
  select(-gene) %>%
  filter(Gene != '') #this dataset has some rows with no Gene listed


ucec <- read.table('proteomics/UCEC_tumor.cct', sep = '\t', header = T, check.names = F) %>%
  pivot_longer(cols = -idx, names_to = 'Sample', values_to = 'Abundance') %>%
  mutate(type = 'UCEC') %>%
  mutate(Gene = idx) %>%
  select(-idx)

rcc <- read.table('proteomics/RCC_tumor.cct', sep = '\t', header = T, check.names = F) %>%
  pivot_longer(cols = -gene, names_to = 'Sample', values_to = 'Abundance') %>%
  mutate(type = 'RCC') %>%
  mutate(Gene = gene) %>%
  select(-gene)

hnsc <- read.table('proteomics/HNSC_tumor.cct', sep = '\t', header = T, check.names = F) %>%
  pivot_longer(cols = -Index, names_to = 'Sample', values_to = 'Abundance') %>%
  mutate(type = 'HNSC') %>%
  mutate(Gene = Index) %>%
  select(-Index)

#these two are a mix of tumour and normal
hcc <- read.xlsx('proteomics/HCC.xlsx', sheet = 4, check.names = F) %>%
  select(starts_with('T'), Sample.ID.Protein) %>%
  pivot_longer(cols = -Sample.ID.Protein, names_to = 'Sample', values_to = 'Abundance') %>%
  mutate(type = 'HCC') %>%
  mutate(Gene = Sample.ID.Protein) %>%
  select(-Sample.ID.Protein)

luad <- read.table('proteomics/luad_proteomics.tsv', sep = '\t', header = T, check.names = F, na.strings = c(""))
luad <- luad[!luad$Name %in% c('Patient_ID', 'Database_ID'),] %>%
	column_to_rownames(loc = 1) %>% 
	t() %>%
	as.data.frame() %>% 
	rownames_to_column(var = 'Gene') %>% 
	pivot_longer(cols = -Gene, names_to = 'Sample', values_to = 'Abundance', values_drop_na = TRUE) %>%
	filter(!grepl('N$', Sample)) %>%
	mutate(Abundance = as.numeric(Abundance)) %>%
	group_by(Gene, Sample) %>% 
	summarize(Abundance = median(Abundance)) %>% #there should be no NAs
	mutate(type = 'LUAD') 

#read in CPTAC normals
n_coad <- read.table('proteomics/COAD_normal.cct', sep = '\t', header = T, check.names = F) %>%
  pivot_longer(cols = -attrib_name, names_to = 'Sample', values_to = 'Abundance') %>%
  mutate(type = 'COAD') %>%
  mutate(Gene = attrib_name) %>%
  select(-attrib_name)

n_ov <- read.table('proteomics/ovarian_normal.cct', sep = '\t', header = T, check.names = F) %>%
  pivot_longer(cols = -gene, names_to = 'Sample', values_to = 'Abundance') %>%
  mutate(type = 'OV') %>%
  mutate(Gene = gene) %>% 
  select(-gene) %>%
  filter(Gene != '')

n_ucec <- read.table('proteomics/UCEC_normal.cct', sep = '\t', header = T, check.names = F) %>%
  pivot_longer(cols = -idx, names_to = 'Sample', values_to = 'Abundance') %>%
  mutate(type = 'UCEC') %>%
  mutate(Gene = idx) %>%
  select(-idx)

n_rcc <- read.table('proteomics/RCC_normal.cct', sep = '\t', header = T, check.names = F) %>%
  pivot_longer(cols = -gene, names_to = 'Sample', values_to = 'Abundance') %>%
  mutate(type = 'RCC') %>%
  mutate(Gene = gene) %>%
  select(-gene)

n_hnsc <- read.table('proteomics/HNSC_normal.cct', sep = '\t', header = T, check.names = F) %>%
  pivot_longer(cols = -Index, names_to = 'Sample', values_to = 'Abundance') %>%
  mutate(type = 'HNSC') %>%
  mutate(Gene = Index) %>%
  select(-Index)

#these two are a mix of tumour and normal
n_hcc <- read.xlsx('proteomics/HCC.xlsx', sheet = 4, check.names = F) %>%
  select(starts_with('N'), Sample.ID.Protein) %>% #extract normals
  pivot_longer(cols = -Sample.ID.Protein, names_to = 'Sample', values_to = 'Abundance') %>%
  mutate(type = 'HCC') %>%
  mutate(Gene = Sample.ID.Protein) %>%
  select(-Sample.ID.Protein)

n_luad <- read.table('proteomics/luad_proteomics.tsv', sep = '\t', header = T, check.names = F, na.strings = c(""))    
n_luad <- n_luad[!n_luad$Name %in% c('Patient_ID', 'Database_ID'),] %>%
        column_to_rownames(loc = 1) %>% 
        t() %>%
        as.data.frame() %>%
        rownames_to_column(var = 'Gene') %>%
        pivot_longer(cols = -Gene, names_to = 'Sample', values_to = 'Abundance', values_drop_na = TRUE) %>%
        filter(!grepl('N$', Sample)) %>%
        mutate(Abundance = as.numeric(Abundance)) %>% 
        group_by(Gene, Sample) %>% 
        summarize(Abundance = median(Abundance)) %>% #there should be no NAs
        mutate(type = 'LUAD') 

#Read in CCLE proteomics
print('READ IN AND PROCESS CCLE')
ccle_orig = read.csv(gzfile('proteomics/CCLE_protein_quant.tsv.gz'), header = T) #file is gzipped
ccle_orig$count_na <- rowSums(is.na(ccle_orig)) #measure missingness. This will also count any NAs in metadata columns

samp <- read.xlsx('proteomics/CCLE_samples.xlsx', sheet = 2) #CCLE sample annotations
bridge_lines <- samp %>% filter(Notes == 'Bridge line') %>% pull(CCLE.Code)
samp <- samp %>% distinct(CCLE.Code, .keep_all = T) #deduplicate sample annotations by Line

ccle <- ccle_orig %>%
  filter(Gene_Symbol != '') %>% 
  arrange(count_na) %>% #deal with duplicate gene measurements
  distinct(Gene_Symbol, .keep_all = T) %>%
  select(-HCT15_LARGE_INTESTINE_TenPx18, -CAL120_BREAST_TenPx02, -SW948_LARGE_INTESTINE_TenPx11) %>% #remove duplicate samples (measured in different plexes)
  setNames(gsub("_TenPx.*","",names(.))) %>% #remove the TenPx ID from the end of Lines
  select(!starts_with('TenPx')) %>% #remove unnecessary columns
  select(-Group_ID, -Uniprot, -Uniprot_Acc, -Description, -Protein_Id, -Uniprot) %>%
  pivot_longer(cols = -Gene_Symbol, names_to = 'Line', values_to = 'Abundance') %>%
  filter(Line %in% samp$CCLE.Code) %>%
  left_join(samp, by = c('Line' = 'CCLE.Code')) %>% #join to the annotations to get tissue type
  filter(!Line %in% bridge_lines) %>% #remove bridge line samples (both copies)
  filter(!is.na(Abundance)) #remove NAs

#get pooled normals and compute logFCs for CPTAC
print('GET POOLED NORMALS AND LOGFCS')
tumours <- bind_rows(coad, ov, ucec, rcc, hnsc, hcc, luad) %>% filter(!is.na(Abundance))
normals <- bind_rows(n_coad, n_ov, n_ucec, n_rcc, n_hnsc, n_hcc, n_luad) %>% filter(!is.na(Abundance))
pooled_normal <- normals %>% group_by(type, Gene) %>% summarize(avg_abundance = median(Abundance))

matched <- tumours %>% 
  mutate(match_normals = case_when(type == 'LUAD' ~ paste0(Sample, '_n'), type == 'HCC' ~ gsub('T', 'N', Sample), type == 'COAD' ~ paste0(Sample, 'N'), TRUE ~ Sample)) %>%
  left_join(normals, by = c('match_normals' = 'Sample', 'Gene'), suffix = c('_tumour', '_normal')) %>%
  left_join(pooled_normal, by = c('type_tumour' = 'type', 'Gene')) %>%
  filter(!is.na(avg_abundance)) %>% #genes not covered in any sample in a given tissue type will be missing from pooled_normal and show up as NA here
  mutate(pooled_log2FC = Abundance_tumour - avg_abundance) %>%
  mutate(unlog = 2**pooled_log2FC)

write.table(matched, file = 'tumours_prot.txt', sep = '\t', col.names = T, row.names = F, quote = F)
#make pooled normals and compute log2FCs for CCLE
arms <- read.xlsx('CN/CCLE_aneuploidy.xlsx', sheet = 2) #arm-level aneuploidy calls to get ploidy

near_dips <- arms %>% #identify near-diploids to use as pseudonormals
  filter(CCLE_ID %in% ccle$Line) %>%
  filter(ploidy > 1.9, ploidy < 2.1) #78 samples

near_dips_lung <- arms %>% #near diploids specifically from lung tissue
  filter(CCLE_ID %in% ccle$Line) %>%
  filter(ploidy > 1.9, ploidy < 2.1) %>%
  filter(grepl('LUNG', CCLE_ID))  #8 samples 

normal_ccle <- ccle %>% 
  filter(Line %in% near_dips$CCLE_ID) %>% 
  group_by(Gene_Symbol) %>%
  summarize(avg_abundance = median(Abundance))

normal_lung <- ccle %>% 
  filter(Line %in% near_dips_lung$CCLE_ID) %>% 
  group_by(Gene_Symbol) %>%
  summarize(avg_abundance = median(Abundance))

ccle_df <- ccle %>% #calculate log2FCs against pooled data -- used to compare tissue-specific v non-specific pooled normals
  filter(!Line %in% near_dips$CCLE_ID) %>% #might be aneuploidy but point is they're not in the controls
  left_join(normal_ccle, by = c('Gene_Symbol')) %>%
  mutate(log2FC_pooled = Abundance - avg_abundance) %>%
  left_join(normal_lung, by = c('Gene_Symbol'), suffix = c('_all', '_lung')) %>%
  mutate(log2FC_lung = Abundance - avg_abundance_lung)

write.table(ccle_df, file = 'ccle_abundances.txt', sep = '\t', col.names = T, row.names = F, quote = F)

print('CALCULATE STOICHIOMETRIC DEVIATIONS')
deviations_cptac <- complex %>%
  select(components_genesymbols, name, n_members) %>%
  separate_rows(components_genesymbols, sep = '_') %>%
  left_join(matched, by = c('components_genesymbols' = 'Gene')) %>%
  filter(!is.na(Sample)) %>%
  add_count(name, Sample) %>%
  mutate(perc = n / n_members) %>%
  filter(perc >= 0.75) %>%
  #filter(!is.na(log2FC)) %>% #use this only if using a non-pooled normal
  group_by(name, Sample, n_members, perc, type_tumour) %>%
  summarize(cv = sd(unlog)/mean(unlog))


deviations_ccle <- complex %>%
  select(components_genesymbols, name, n_members) %>%
  separate_rows(components_genesymbols, sep = '_') %>%
  left_join(ccle_df, by = c('components_genesymbols' = 'Gene_Symbol')) %>% 
  filter(!is.na(Line)) %>%
  add_count(name, Line) %>%
  mutate(perc = n / n_members) %>%
  filter(perc >= 0.75) %>%
  mutate(unlog = 2**log2FC_pooled) %>%
  group_by(name, Line, n_members, perc) %>%
  summarize(cv = sd(unlog)/mean(unlog))

write.table(deviations_cptac, file = 'deviations_cptac.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(deviations_ccle, file = 'deviations_ccle.txt', sep = '\t', col.names = T, row.names = F, quote = F)
#summarize(rang = max(unlog) - min(unlog), ks = ks.test(unlog, punif(length(unlog), min = min(unlog), max = max(unlog)))$statistic, std = sd(unlog), me = mean(unlog), entro = entropy(unlog), gini = ineq(unlog, type = 'Gini'), RMAD = rmad(unlog)) %>%
#mutate(cv = std/me)

#chi-square not valid and gives warnings but it looked like "chi = chisq.test(unlog)$statistic"
