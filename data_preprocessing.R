library(jsonlite)
library(xlsx)
library(dplyr)

########### Initialization ############

# Set working directory
setwd("~/Documents/Alzheimer-Disease")

# Import settings
settings <- jsonlite::read_json("settings.json")

### Merge HLA calls and write 

# Files name
HLACalls_names <- list.files(settings$folder$HLACalls)

# Initialize 
HLA.df <- read.csv(paste0(settings$folder$HLACalls,HLACalls_names[1])); HLA.df <- HLA.df[,1:3]
colnames(HLA.df) <- c('SampleID', 'A.1', 'A.2')

# Merge
for (file in HLACalls_names[2:length(HLACalls_names)]){
  
  # Load and merge
  df <- read.csv(paste0(settings$folder$HLACalls, file))
  df <- df[,1:3]
  HLA.df <- merge(HLA.df, df, by.x = 'SampleID', by.y = 'sample.id')
  
}

# Write 
write.csv(HLA.df, file = settings$file$HLA_df, row.names = FALSE)

###

############## PRE-PROCESSING #############

### Merge relatedness and write 

# Files name 
relatedness_names <- list.files(settings$folder$Relatedness)

# Initialize 
rtd.df <- read.csv(paste0(settings$folder$Relatedness, relatedness_names[1]), sep = '\t')

# Concat
for (file in relatedness_names){
  df <- read.csv(paste0(settings$folder$Relatedness, file), sep = '\t')
  rtd.df <- rbind(rtd.df, df)
}

# Write 
write.csv(rtd.df, file = settings$file$relatedness_df)



