# Import libraries
library(jsonlite)
library(tidyverse)
library(readr)
library(data.table)
library(ggrepel)
library(viridis)
library(hrbrthemes)
library(HIBAG)
library(parallel)
library(corrplot)
library(randomForest)
library(xlsx)
library(TreeBH)
library(zeallot)
library(epitools)

########## IMPORT #########

# Set working directory
setwd("~/Documents/Alzheimer-Disease")

# Import settings
settings <- jsonlite::read_json("settings.json")

# Set values to 0
settings$controlAlleles <- c()
settings$excludePosSubjsByAllele <- c()
settings$allele2Remove <- c()

# Create comand
`%notin%` <- Negate(`%in%`)

# Import HLA calls and covariates
HLA.df <- read.csv(settings$file$HLA_df)

# Import amino acid alignment 
AA_alignment <- read.table(settings$file$AA_alignment, header = TRUE, sep = ',')

# Remove negative positions 
full_prot = FALSE
if (!full_prot){
  AA_alignment[,3] <- AA_alignment %>% apply(MARGIN = 1, function(x) x[3] %>% substr(start = 30, stop = nchar(x[3])))
}

########### DATA TO OHE ###############

data2OHE = function(settings, data, allelesAA, covars.df, AA_alignment){
  
  # Get data with AAs
  data.AApos <- data %>% filter(get(A1) %in% allelesAA| get(A2) %in% allelesAA)
  data.AAneg <- data %>% filter(sample.id.file %notin% data.AApos$sample.id.file)
  
  # Create dataframe with OHE
  data.AA <- data.frame(sample.id.file = c(data.AApos$sample.id.file, data.AAneg$sample.id.file), AA = c(rep(1, nrow(data.AApos)), rep(0, nrow(data.AAneg))))
  
  # Merge with presence of controlled amino acid
  if (!settings$AA2control %>% is_empty()){
    data.AAcontrol <- controlAA(settings, data, AA_alignment)
    data.AA <- merge(data.AA, data.AAcontrol, by = "sample.id.file")
  }
  
  # Return 
  return(data.AA)
}

################## GET OHE ##############

# Loci
L <- "DRB1"

# Get locus ID
c(A1, A2) %<-% c(paste0(L,'.1'), paste0(L,'.2'))

# Get alignment subset, and counts subset
AA_locus <- AA_alignment %>% filter(locus == L)

# Get max sequence length 
maxLen <- AA_locus$sequence %>% lapply(nchar) %>% unlist() %>% max()

# Get unique AAs 
AAOHE <- data.frame("sample.id" = HLA.df$sample.id, "sample.id.file" = HLA.df$sample.id.file)
for (pos in 1:maxLen){
  
  # Get unique AAs
  AAs <- AA_locus$sequence %>% lapply(function(x, pos) substr(x, pos, pos), pos) %>% unlist() %>% unique()
  AAs <- AAs[which(AAs != '' & AAs!= '*')]
  
  # For each AA
  for (AA in AAs){
    
    # Get alleles with AAs
    allelesAA.OG <- AA_locus %>% apply(MARGIN = 1, function(x, pos, AA) if (substr(x[3],pos,pos) == AA) {return(x[2])}, pos, AA) %>% unlist()
    allelesAA <- allelesAA.OG %>% lapply(function(x) strsplit(x, split='\\*') %>% unlist() %>% 
                                           .[2] %>% strsplit(split=':') %>% unlist() %>% .[1:2] %>% paste(collapse=':')) %>% unlist() 
    
    # Crete one hot encoding dataframes
    data.AA <- data2OHE(settings, HLA.df, allelesAA, covars.df, AA_alignment); colnames(data.AA)[2] <- paste(L,pos,AA, sep = '_')
    
    # Merge
    AAOHE <- merge(AAOHE, data.AA, by ="sample.id.file")
    
  }

}

# Wrie 
write.table(AAOHE, file = "Outputs/OHE_AA/OHE_AA.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE) 
