library(jsonlite)
library(xlsx)
library(dplyr)

########### Initialization ############

# Set working directory
setwd("~/Documents/Alzheimer-Disease")

# Import settings
settings <- jsonlite::read_json("settings.json")

# Create comand
`%notin%` <- Negate(`%in%`)

# Import covariates
covars.df <- read.csv(settings$file$HLA_covars, sep = '\t')
covars.df$file.origin <- covars.df$FID %>% lapply(function(x) gsub('[[:digit:]]+', '', x %>% strsplit('_') %>% unlist() %>% .[1])) %>% unlist()
covars.df$sample.id.file <- paste0(covars.df$file.origin,'_', covars.df$IID)

# Import relatedness
relat.df <- read.csv(settings$file$relatedness_df)

############ REMOVE DUPLICATE IDS AND MERGE ###############

### Remove duplicate Ids
# Load HLA file 
HLACalls_names <- list.files(settings$folder$HLACalls)
HLA.df <- read.csv(paste0(settings$folder$HLACalls,HLACalls_names[1])); 

# Create new variables 
HLA.df$file.origin <- HLA.df$fileName %>% lapply(function(x) gsub('[[:digit:]]+', '', x %>% strsplit('_') %>% unlist() %>% .[1])) %>% unlist()
HLA.df$sample.id.file <- paste0(HLA.df$file.origin,'_', HLA.df$sample.id)
HLA.df <- HLA.df[c(1,7,8,2,3)]

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
  HLA.df <- merge(HLA.df, df_loop[-c(1:2)], by.x = 'sample.id.file', by.y = 'sample.id.file', all = TRUE)
  
}

# CHECK: NAs
missHLA <- HLA.df[which(is.na(HLA.df), arr.ind=TRUE),]

# Filter file to get subjects in covars file
HLA.df.filt <- HLA.df[which(HLA.df$sample.id.file %in% covars.df$sample.id.file),]

# Remove duplicates based on alleles, and sample.id.file
HLA.df.filt <- HLA.df.filt %>% distinct(sample.id.file, .keep_all = TRUE)
covars.df <- covars.df %>% distinct(sample.id.file, .keep_all = TRUE)

# Filter covars 
covars.df <- covars.df[which(covars.df$sample.id.file %in% HLA.df.filt$sample.id.file),]

############# REMOVE DUPLICATE INDIVIDUALS ################

# Filter relatedness
related.df.filt <- relat.df %>% filter(ID1 %in% HLA.df.filt$sample.id & ID2 %in% HLA.df.filt$sample.id)

# Count duplicates
rel.dup.df <- related.df.filt %>% filter(InfType == 'Dup/MZ')

# Remove duplicates 
pairs.df <- data.frame(id1 = c(), id2 = c())
for (i in 1:nrow(rel.dup.df)){
  
  # Get ID1 and ID2 
  id1 <- paste0(rel.dup.df$FID1 %>% lapply(function(x) gsub('[[:digit:]]+', '', x %>% strsplit('_') %>% unlist() %>% .[1])) %>% unlist() %>% .[i],
                "_",rel.dup.df$ID1[i]); 
  id2 <- paste0(rel.dup.df$FID2 %>% lapply(function(x) gsub('[[:digit:]]+', '', x %>% strsplit('_') %>% unlist() %>% .[1])) %>% unlist() %>% .[i]
                ,"_",rel.dup.df$ID2[i]); 
  
  # Filter 
  df <- HLA.df.filt %>% filter(sample.id.file %in% c(id1, id2))
  
  # Check if identical 
  if (any(duplicated(df[,c(4:23)]) == TRUE)){
    # Keep pair 
    pairs.df <- rbind(pairs.df, data.frame(id1 = id1, id2 = id2))
  }
  
}

# Get duplicated subjects' ids and remove
dupSubj <- pairs.df$id1 %>% unique()
HLA.df.filt <- HLA.df.filt %>% filter(sample.id %notin% dupSubj)
covars.df <- covars.df %>% filter(IID %notin% dupSubj)

# Write 
write.csv(HLA.df.filt, file = settings$file$HLA_df, row.names = FALSE)
write.csv(covars.df, file = settings$file$HLA_covars_filt, row.names = FALSE)

################# SKETCH #################
idx <- 145
for (i in 1:nrow(HLA.df.filt)){
  
  thr <- c()
  for (col in 4:23){
    if (HLA.df.filt[idx,col] == HLA.df.filt[i,col]){
      thr <- c(thr, TRUE)
    } else{
      thr <- c(thr, FALSE)
    }
  }
  
  if (all(thr)){
    x <- HLA.df.filt[i,]
  }
}
