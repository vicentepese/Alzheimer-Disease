library(jsonlite)
library(xlsx)
library(dplyr)
library(gridExtra)
library(ggplot2)

########### Initialization ############

# Set working directory
setwd("~/Documents/Alzheimer-Disease")

# Import settings
settings <- jsonlite::read_json("settings.json")

# Create comand
`%notin%` <- Negate(`%in%`)

# Import covars
covars.df <- read.table(settings$file$covars_UKB, sep = '\t',header = TRUE)

# Import pheno
pheno.df <- read.table(file = settings$file$pheno_UKB, header = TRUE)
pheno.df$pheno <- with(pheno.df, ifelse(AD_proxy > 1, 1, 0))

# Import set of individuals
indivs.set <- read.table(file = settings$file$individuals_UBK, header = TRUE)

# Filter covars and pheno by set of individuals
covars.df <- covars.df %>% filter(IID %in% indivs.set$IID)
pheno.df <- pheno.df %>% filter(IID %in% indivs.set$IID)


########### MERGE HLA #############

# List of files to merge 
file.names <- list.files(settings$folder$HLACalls_UKB) 

# First pass
f <- file.names[1]
HLA.data.UK <- read.csv(paste0(settings$folder$HLACalls_UKB, f))

# Get probs and alleles
probs.UK <- HLA.data.UK[,c(1,4)]
HLA.data.UK <- HLA.data.UK[,1:3]

# Merge files
for (f in file.names[2:length(file.names)]){
  
  # Print for loop
  print(paste0("Current file: ", f))
  
  # Load file 
  data <- read.csv(paste0(settings$folder$HLACalls_UKB, f))
  
  # Merge by sample.id
  HLA.data.UK <- merge(HLA.data.UK, data[,1:3], by = "sample.id")
  
  # Merge probs by sample.id
  probs.UK <- merge(probs.UK, data[,c(1,4)], by = "sample.id")
  
}

# Compute average probabilities
probs.UK$mu_prob <- apply(probs.UK[,2:10], MARGIN = 1, FUN = function(x) mean(x))

# Filter by set of individuals 
HLA.data.UK <- HLA.data.UK %>% filter(sample.id %in% indivs.set$FID)
probs.UK <- probs.UK %>% filter(sample.id %in% indivs.set$FID)
  
# Filter individuals with probability higher than 0.5
if (settings$filterProb){
  probs.UK <- probs.UK %>% filter(DRB1.prob > 0.5)
  HLA.data.UK <- HLA.data.UK %>% filter(sample.id %in% probs.UK$sample.id)
}

# Merge pheno into covars
covars.df <- merge(pheno.df[,c("IID", "AD_proxy", "pheno")], covars.df[,2:ncol(covars.df)])
covars.df <- covars.df %>% filter(IID %in% HLA.data.UK$sample.id)

# Add variables for consistency in the pipeline 
HLA.data.UK$sample.id.file <- HLA.data.UK$sample.id
HLA.data.UK$file.origin <- HLA.data.UK$sample.id
covars.df$sample.id <- covars.df$IID
covars.df$sample.id.file <- covars.df$IID

# Write file
write.csv(HLA.data.UK, file = settings$file$HLA_df_UKB, row.names = FALSE)
write.csv(covars.df, file = settings$file$HLA_covars_UKB, row.names = FALSE)

############ PLOTS ############

library(ggplot2)
library(gridExtra)
library(rlist)
pl <- vector('list', ncol(probs.UK)-1)
idx <- 1
for (i in 2:ncol(probs.UK)){
  
  pl[[i-1]] <- local({
    i <- i
    p1 <- ggplot(probs.UK, aes(get(colnames(probs.UK)[i]))) +
      geom_histogram() +
      xlab(colnames(probs.UK)[i]) + xlim(0,1)
    print(p1)
  })
  
  
}
do.call(grid.arrange, pl)
