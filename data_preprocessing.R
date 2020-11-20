library(jsonlite)
library(xlsx)
library(dplyr)

########### INITIALIZATION ############

# Set working directory
setwd("~/Documents/Alzheimer-Disease")

# Import settings
settings <- jsonlite::read_json("settings.json")

# Create comand
`%notin%` <- Negate(`%in%`)

# Import covariates and create new variable
covars.df <- read.csv(settings$file$HLA_covars, sep = '\t')
covars.df$sample.id.file <- covars.df$FID %>% lapply(function(x) gsub('[[:digit:]]+', '', x %>% strsplit('_') %>% unlist() %>% .[1])) %>% unlist()
covars.df$sample.id.file <- covars.df$sample.id.file %>% paste0(rep('_', nrow(covars.df)), covars.df$IID)

# Import relatedness and create new variables
relat.df <- read.csv(settings$file$relatedness_df)
relat.df <- relat.df[,-1]
relat.df$sample.id.file1 <- relat.df$FID1 %>% lapply(function(x) gsub('[[:digit:]]+', '', x %>% strsplit('_') %>% unlist() %>% .[1])) %>% unlist()
relat.df$sample.id.file1 <- relat.df$sample.id.file1 %>% paste0(rep('_', nrow(relat.df)), relat.df$ID1)
relat.df$sample.id.file2 <- relat.df$FID2 %>% lapply(function(x) gsub('[[:digit:]]+', '', x %>% strsplit('_') %>% unlist() %>% .[1])) %>% unlist()
relat.df$sample.id.file2 <- relat.df$sample.id.file2 %>% paste0(rep('_', nrow(relat.df)), relat.df$ID2)



############ REMOVE DUPLICATE IDS AND MERGE ###############

### Remove duplicate Ids
# Load HLA file 
HLACalls_names <- list.files(settings$folder$HLACalls)
HLA.df <- read.csv(paste0(settings$folder$HLACalls,HLACalls_names[1])); 

# Create new variables
HLA.df$file.origin <- HLA.df$fileName %>% lapply(function(x) gsub('[[:digit:]]+', '', x %>% strsplit('_') %>% unlist() %>% .[1])) %>% unlist()
HLA.df$sample.id.file <- paste0(HLA.df$file.origin,'_', HLA.df$sample.id)

# Filter file to get subjects in covars file
HLA.df.filt <- HLA.df %>% filter(sample.id.file %in% covars.df$sample.id.file)

# Remove duplicates based on alleles, and sample.id.file
dup_ids <- HLA.df.filt[which(duplicated(HLA.df.filt$sample.id.file)),]
dup_ids$sample.id <-
HLA.df.filt <- HLA.df.filt %>% distinct(sample.id.file, .keep_all = TRUE)
HLA.df.filt <- HLA.df.filt[c(1,7,8,2,3)]

# Merge
for (file in HLACalls_names[2:length(HLACalls_names)]){
  
  print(paste0('Current iteration: ', file))
  
  # Load and create new variables
  df_loop <- read.csv(paste0(settings$folder$HLACalls, file)); 
  df_loop$file.origin <- df_loop$fileName %>% lapply(function(x) gsub('[[:digit:]]+', '', x %>% strsplit('_') %>% unlist() %>% .[1])) %>% unlist()
  df_loop$sample.id.file <- paste0(df_loop$file.origin,'_', df_loop$sample.id)
  
  # Remove duplcates
  df_loop <- df_loop %>% distinct(sample.id.file, .keep_all = TRUE)
  df_loop <- df_loop[c(1,7,8,2,3)]
  
  # Merge
  HLA.df.filt <- merge(HLA.df.filt, df_loop[-c(1:2)], by.x = 'sample.id.file', by.y = 'sample.id.file')
  
}

# Write 
write.csv(HLA.df.filt, file = settings$file$HLA_df, row.names = FALSE)

# Filter covars based on HLA.df.filt ids, remove duplicates 
covars.df.filt <- covars.df %>% filter(sample.id.file %in% HLA.df.filt$sample.id.file)
covars.df.filt <- covars.df.filt %>% distinct(sample.id.file, .keep_all = TRUE)

# Write 
write.csv(covars.df.filt, file = settings$file$HLA_covars_filt)


############# REMOVE DUPLICATE INDIVIDUALS ################

# Read pre-process HLA calls
HLA.df.filt <- read.csv(settings$file$HLA_df)

# Filter relatedness
related.df.filt <- relat.df %>% filter(sample.id.file1 %in% HLA.df.filt$sample.id.file & sample.id.file2 %in% HLA.df.filt$sample.id.file)

# Count duplicates
rel.dup.df <- related.df.filt %>% filter(InfType == 'Dup/MZ')

# Remove duplicates 
pairs.df <- data.frame(id1 = c(), id2 = c())
for (i in 1:nrow(rel.dup.df)){
  
  # Get ID1 and ID2 
  id1 <- rel.dup.df$sample.id.file1[i]; id2 <- rel.dup.df$sample.id.file2[i]
  
  # Filter 
  df <- HLA.df.filt %>% filter(sample.id.file %in% c(id1, id2))
  
  # Check if identical 
  if (any(duplicated(df[,c(4:23)]) == TRUE)){
    # Keep pair 
    pairs.df <- rbind(pairs.df, data.frame(id1 = id1, id2 = id2))
  }
  
}

# Sanity check with covars
check_idx <- c()
for (i in 1:nrow(pairs.df)){
  
  # Subset 
  cov_check <- covars.df %>% filter(sample.id.file %in% pairs.df[i,])
  
  if (!any(cov_check[,c(3,14,15)] %>% duplicated() ==TRUE)){
    check_idx <- c(check_idx, i)
  }
  
}

# Get duplicated subjects' ids and remove
dupSubj <- pairs.df$id1 %>% unique()
HLA.df.filt <- HLA.df.filt %>% filter(sample.id.file %notin% dupSubj)
covars.df.filt <- covars.df.filt %>% filter(sample.id.file %notin% dupSubj)

################### CHECK IDENTICAL HLA CALLS ############

# Check identical HLAs 
idx_iden <- which(duplicated(HLA.df.filt[,c(4:23)]) == TRUE)
IDs_iden <- HLA.df.filt$sample.id.file[idx_iden]
Dx_prop_iden <- covars.df.filt %>% filter(sample.id.file %in% IDs_iden) %>% .["pheno"] %>% table()

# Remove identical HLA calls and sanity checks
HLA.df.filt <- HLA.df.filt %>% filter(sample.id.file %notin% IDs_iden)
which(duplicated(HLA.df.filt[,c(4:23)]) == TRUE)
covars.df.filt <- covars.df.filt %>% filter(sample.id.file %notin% IDs_iden)
covars.df.filt$pheno %>% table() %>% prop.table()*100

# Write 
write.csv(HLA.df.filt, settings$file$HLA_df, row.names = FALSE)
write.csv(covars.df.filt, settings$file$HLA_covars_filt, row.names = FALSE)


#########SCRATCH###########
idx <- c()
for (i in 1:nrow(dup.df)){
  
  if (identical(c(dup.df[i,c(4:23)]) %>% unlist() %>%  unname(), tstVec)){
    idx <- c(idx, i)
  }
  
}
