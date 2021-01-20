## ---------------------------
##
## Script name: merge_HLA.R
##
## Purpose of script: Merge files from different datasets provided by Stanford
##
## Author: Vicente Peris Sempere
##
## Date Created: 12/21/2020
##
## Copyright (c) Vicente Peris Sempere, 2020
## Email: vipese@stanford.edu
##
## ---------------------------
##
## Notes: None
##   
##
## ---------------------------

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

########### MERGE FILES ############

# Get list of files
files = list.files(settings$folder$HLACalls_Individuals)

# Get loci
loci <- files %>% lapply(function(x) x %>% strsplit("_") %>% unlist() %>% tail(n=1) %>% 
                           strsplit("[.]") %>% unlist() %>% .[1]) %>% unlist() %>% unique()

# For each locus
for (locus in loci){
  
  # Initialize sub-loop
  locus.df <- data.frame()
  
  # For each file
  for (file in files){
    
    # If file has locus
    if (file %>% lapply(function(x) x %>% strsplit("_") %>% unlist() %>% tail(n=1) %>% 
                    strsplit("[.]") %>% unlist() %>% .[1]) %>% unlist() %>% unique() == locus){
      
      # Read data and bind row-wise
      data <- read.csv(paste0(settings$folder$HLACalls_Individuals, file))
      data$fileName <- rep(file, nrow(data))
      locus.df <- rbind(locus.df, data)
      
    }
  }
  
  # Write final result
  write.csv(x = locus.df, file = paste0(settings$folder$HLACalls, "HLA_", locus, ".csv"), row.names = FALSE)
  
}


