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
library(plyr)
library(TreeBH)

## DESCRIPTION ##
# This scripts uses one-hot encoding to predict with a logistic regression the 
# most significant alleles in each locus in patients where sample size is > 45 

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
covars.df <- read.csv(settings$file$HLA_covars_filt)
covars.df$pheno <- covars.df$pheno -1
covars.df$sex <- covars.df$sex -1


########### ONE HOT ENCODING FUNCTIONS ###############

# HLA Parse function
alleleFreqOHE=function(test_DF){
  test_DF = as.data.frame(test_DF)
  sample.id.file = rep(test_DF[,1], 2)
  make_HLA = c(as.character(test_DF[,2]), as.character(test_DF[,3]))
  make_DF = cbind.data.frame(sample.id.file, make_HLA)
  setDT(make_DF)
  dcast_HLA = dcast(make_DF, sample.id.file~make_HLA, fun.aggregate = length)
  dcast_HLA$sample.id.file = as.factor(dcast_HLA$sample.id.file)
  return(setDF(dcast_HLA))
}

# HLA Parse function
carrierFreqOHE=function(test_DF){
  test_DF = as.data.frame(test_DF)
  sample.id.file = rep(test_DF[,1], 2)
  make_HLA = c(as.character(test_DF[,2]), as.character(test_DF[,3]))
  make_DF = cbind.data.frame(sample.id.file, make_HLA)
  setDT(make_DF)
  dcast_HLA = dcast(make_DF, sample.id.file~make_HLA, fun.aggregate = length)
  dcast_vals <- dcast_HLA[,!c("sample.id.file")]; dcast_vals[dcast_vals > 1] <- 1
  dcast_HLA <- cbind(data.frame(sample.id.file =dcast_HLA$sample.id.file), dcast_vals)
  dcast_HLA$sample.id.file = as.factor(dcast_HLA$sample.id.file)
  return(setDF(dcast_HLA))
}

########### PRE-PROCES DATA ############
# Compute a logistic regression to test significance of carrier frequencies 

## Exclude subjects w/ positive allele 
excludeByAllele = function(data.filt, settings){
  
  # Alleles to control
  As2exclude = settings$excludePosSubjsByAllele
  
  exclude_list <- c()
  for (A in As2exclude){
    
    # Get locus and allele 
    locus <- A %>% strsplit('\\*') %>% unlist() %>% .[1]
    A1 <- paste('HLA', locus, '_A1', sep = ''); A2 <- paste('HLA', locus, '_A2', sep = '')
    allele2exclude = A %>% strsplit('\\*') %>% unlist() %>% .[2]
    
    # Get data 
    data.filt <- data.filt %>% filter(get(A1) != allele2exclude & get(A2) != allele2exclude)
    
  }
  return(data.filt)
}

preprocessData=function(cases, controls){
  
  # Import subject list 
  subjects.ID <- read.table(settings$file$patList, colClasses = c('character', 'NULL')) %>% unlist()
  changeIdFormat = function(GWASID){
    idpre <- GWASID %>% strsplit('_')  %>% unlist() 
    return(strsplit(idpre[length(idpre)], '\\.')[[1]][1])
  }
  subjects.GWASID <- sapply(subjects.ID, changeIdFormat, USE.NAMES = FALSE)
  
  # Filter out patients
  data.filt <- data[which(data$sample.id %in% subjects.ID),]
  
  # Add diagnosis to data
  Dx <- c()
  for (subj in as.character(data.filt$GWASID)){
    if (subj %in% controls){
      Dx <- c(Dx, 0)
    } else if (subj %in% cases) {
      Dx <- c(Dx, 1)
    } else{
      Dx <- c(Dx, -9)
    }
  }
  
  # Add to datframe 
  data.filt['Dx'] <- Dx
  
  # Exclude subjects with allele
  data.filt <- excludeByAllele(data.filt, settings)
  
  # Return
  return(data.filt)
  
}

########### ALLELE AND CARRIER FREQUENCY  #################


computeACFREQ = function(data, locus, Dx){
  # Compute heterozigous, homozigous, or absence count. 
  # Compute allele frequency, count and total.
  # Compute carrier frequency, count and total.
  
  # Carrier frequency
  A1 <- locus %>% paste0('.1'); A2 <- locus %>% paste0('.2')
  alleles <- list(data[, A1], 
                  data[, A2]) %>% unlist()
  carriers <- data[,c(A1, A2)]
  carriers.levels <- list(data[, A1], 
                          data[, A2]) %>% unlist() %>% levels()
  carriers.unique <- apply(carriers, 1, function(x) unique(x)) %>% unlist() %>% as.factor()
  carriers.count <- table (carriers.unique); carriers.count[c(carriers.levels %>% setdiff(carriers.count %>% names()))] <- 0
  carriers.freq <- carriers.count /nrow(data) * 100
  carrier.df <- data.frame(allele = carriers.count %>% names(), 
                           carrierCount = carriers.count %>% as.vector(),
                           carrierFreq = carriers.freq %>% as.vector(),
                           carrierTotal = nrow(data))
  
  # Heterozigous, homozigous, and absence count
  A0 <- c(); A1 <- c(); A2 <- c();
  for (A in levels(as.factor(alleles))){
    HH.data <- carriers[which(carriers[,1]==as.character(A) | carriers[,2]==as.character(A)),]
    HH.count <- HH.data %>% apply(1, function(x) x %>% unlist() %>% unique() %>% length()) %>% unlist() %>% table()
    A0 <- c(A0, nrow(data) - nrow(HH.data)); A1 <- c(A1,HH.count['2'] %>% unname()); A2 <- c(A2,HH.count['1'] %>% unname())
  }
  A1[is.na(A1)] <- 0; A2[is.na(A2)] <- 0; HH.data <- data.frame(allele = levels(as.factor(alleles)), A0 = A0, A1 = A1, A2 = A2)
  
  # Merge 
  ACFREQ.df <- merge(HH.data, carrier.df, by = 'allele')
  
  switch (Dx,
          'case' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2','carrierCount', 'carrierFreq', 'carrierTotal'),
                                                     rep('Case',4), sep = ''))
            return(ACFREQ.df)
          },
          'control' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2','carrierCount', 'carrierFreq', 'carrierTotal'),
                                                     rep('Control',4), sep = ''))
            return(ACFREQ.df)
          }
  )
}

########### CONTROL ALLELES ##########

controlAllele = function(as2control, HLA.df){
  
  # Get unique loci 
  lociControl = lapply(as2control, function(x) x %>% strsplit('\\*') %>% unlist() %>% .[1]) %>% unlist() %>% unique() 
  
  # For each allele 
  alleleControl.df = data.frame(sample.id.file = HLA.df$sample.id.file)
  for (A in as2control){
    
    # Get locus and allele 
    locus <- A %>% strsplit('\\*') %>% unlist() %>% .[1]
    allele2control = A %>% strsplit('\\*') %>% unlist() %>% .[2]
    allele1 <- paste(locus, '.1', sep = '')
    allele2 <- paste(locus, '.2', sep = '')
    
    # Get subjects
    alleleControl.df[A] <- as.logical(c(HLA.df[,c(allele1)] %>% as.character() == allele2control) + 
                                        c(HLA.df[,c(allele2)] %>% as.character() == allele2control)) %>% as.integer()
    
  }
  return(alleleControl.df)
}


############ REGRESSION MODEL ############

runLogisticRegression = function(locus, OHE.carrierFreq.data, covars.df, as2control = NULL){
  
  ## Allele Frequency 
  # Merge dataset to include PCs
  covars.df$sample.id.file <- as.factor(covars.df$sample.id.file)
  alleles.freq <- colnames(OHE.carrierFreq.data)[-c(1,(ncol(OHE.carrierFreq.data)-length(as2control)):ncol(OHE.carrierFreq.data))]
  OHE.carrierFreq.data <- merge(OHE.carrierFreq.data, covars.df[,-3], by = 'sample.id.file')
  
  # Remove alleles for control 
  locus.subset <- as2control[grepl(as2control, pattern = locus)]
  allele.subset <- locus.subset %>% lapply(function(x) x %>% strsplit('\\*') %>% unlist() %>% .[2]) %>% unlist()
  OHE.carrierFreq.data[allele.subset] <- NULL
  alleles.freq <- alleles.freq[!alleles.freq %in% allele.subset]
  
  # Run logistic regression on carrier frequency 
  Acarrier.model.df <- data.frame()
  for (allele in alleles.freq){
    if (!is.null(as2control)){
      control.alleles <- paste(' ', as2control %>% sapply(function (x) paste('`', x ,'`', sep = '')) %>% paste(collapse = ' + '), sep = '+ ')
    } else{
      control.alleles <- ''
    }
    glm.formula <- paste('pheno ~ `',allele, '` + PC1 + PC2 + PC3', control.alleles, sep = '')
    Acarrier.model <- glm(data = OHE.carrierFreq.data, 
                          formula = as.formula(glm.formula),
                          family = 'binomial', maxit = 100) %>% summary()
    Acarrier.model.df <- rbind(Acarrier.model.df, c(Acarrier.model$coefficients[2,1], 
                                                    Acarrier.model$coefficients[,dim(Acarrier.model$coefficients)[2]]))
  }
  colnames(Acarrier.model.df) <- c('allele.COEF.CARRIER', c('Incercept', 'allele', Acarrier.model$coefficients[-c(1,2),] %>% row.names()) %>%
                                     paste('.CARRIER.pval', sep = ''))
  Acarrier.model.df <- data.frame(allele=alleles.freq, Acarrier.model.df)
  
  # Merge 
  glm.data <- Acarrier.model.df
  
  # Return
  return(glm.data)
  
}

fitGLM = function(settings, locus, HLA.df, data.cases, data.controls, covars.df, as2control = NULL){
  
  # Control for allele
  controlAllele.df <- controlAllele(as2control, HLA.df)
  
  # Subset locus
  allele1 <- paste(locus, '.1',sep = '')
  allele2 <- paste(locus,'.2', sep =  '')
  data.locus <- HLA.df[,c('sample.id.file', allele1, allele2)]
  
  # Compute allele frequencies and coutns, and carrier frequencies and counts
  ACFREQ.cases <- computeACFREQ(data.cases, locus, 'case');
  carrierCases <- unique(ACFREQ.cases$carrierTotalCase)
  ACFREQ.controls <- computeACFREQ(data.controls, locus, 'control');
  carrierControls <- unique(ACFREQ.controls$carrierTotalControl)
  
  # Merge and clean
  ACFREQ.df <- merge(ACFREQ.cases[,!names(ACFREQ.cases) %in% c('A0','A1','A2')], ACFREQ.controls, by = 'allele', all = TRUE) 
  ACFREQ.df[is.na(ACFREQ.df)] <- 0
  ACFREQ.df$carrierTotalCase <- rep(carrierCases, nrow(ACFREQ.df)); ACFREQ.df$carrierTotalControl <- rep(carrierControls, nrow(ACFREQ.df))
  ACFREQ.df <- ACFREQ.df %>% filter(allele != '')
  
  # Parse one hot encoding and merge
  OHE.carrierFreq.data <- carrierFreqOHE(data.locus); if ('V1' %in% colnames(OHE.carrierFreq.data)){OHE.carrierFreq.data$V1 <- NULL}
  OHE.carrierFreq.data<-  merge(as.data.frame(OHE.carrierFreq.data), 
                                covars.df[c('pheno', 'sample.id.file')], by.x ='sample.id.file', by.y = 'sample.id.file') %>% 
    merge(controlAllele.df, by = 'sample.id.file')
  
  # Remove subjects thar are not controls or cases
  OHE.carrierFreq.data <- OHE.carrierFreq.data %>% filter(pheno != -9)
  
  # Run logistic regression model 
  glm.data <- runLogisticRegression(locus, OHE.carrierFreq.data, covars.df, as2control)
  
  # Create dataframes 
  HLA.GLM_carriers.df <-merge(glm.data[,c(1,which(grepl('CARRIER', colnames(glm.data))))],
                              ACFREQ.df[,c(1,which(grepl(paste(c('A0','A1','A2', 'carrier'), collapse = '|'), colnames(ACFREQ.df))))],
                              by = 'allele')
  
  return(HLA.GLM_carriers.df)
}

############### HLA ANALYSIS ###########


# Get cases and controls
cases.ids <- covars.df$sample.id.file[which(covars.df$pheno ==1)]
controls.ids <- covars.df$sample.id.file[which(covars.df$pheno ==0)]
data.cases <- HLA.df %>% filter(sample.id.file %in% cases.ids)
data.controls <- HLA.df %>% filter(sample.id.file %in% controls.ids)

# Initialize while lopp
pval <- 0; as2control <- c(); signAlleles <- list(); 

# HLA Loci
loci <- c('A','B','C','DQA1', 'DQB1', 'DPB1', 'DRB1','DRB3','DRB4','DRB5')

# While signifiant alleles
idx <- 1; iter <- 1
while (pval < 0.05){
  
  # Fit GLM 
  HLA.GLM_carriers.list <- list(); pvalTotal <- c()
  for (locus in loci){
    HLA.GLM_carriers.df <- fitGLM(settings, locus, HLA.df, data.cases, data.controls, covars.df, as2control)
    HLA.GLM_carriers.list[[locus]] <- HLA.GLM_carriers.df 
    for (A in HLA.GLM_carriers.df$allele){
      pvalA <- HLA.GLM_carriers.df$allele.CARRIER.pval[which(HLA.GLM_carriers.df$allele == A)]
      pvalTotal[paste(locus, '*', A, sep = '')] <- pvalA
    }
  }
  
  # Get minimum allele value 
  pval <- pvalTotal[which(pvalTotal == min(pvalTotal))][1]
  pvalMin <- pval %>% names(); pvalMinLocus <- pvalMin %>% strsplit('\\*') %>% unlist() %>% .[1]
  pvalMinAllele <- pvalMin %>% strsplit('\\*') %>% unlist() %>% .[2]
  
  # Add to outputs
  as2control <- c(as2control, pvalMin)
  preOut <- HLA.GLM_carriers.list[[ pvalMinLocus]] %>% filter(allele == pvalMinAllele)
  preOut$allele <- pvalMin  ; signAlleles[[idx]] <- preOut
  
  # Compute allele and locus group for correction
  alleleGroup <- pvalTotal %>% names(); locusGroup <- lapply(alleleGroup, function(x) x %>% strsplit('\\*') %>% unlist() %>% .[1]) %>% unlist()
  groups <- matrix(c(as.factor(locusGroup), as.factor(alleleGroup)), ncol = 2)
  
  # Compute TreeBH
  TreeBH.res <- get_TreeBH_selections(pvals = pvalTotal %>% unname(), groups = groups, q = rep(0.05, ncol(groups)))
  TreeBH.res.DF <- data.frame(allele = alleleGroup, pval = pvalTotal %>% unname()); 
  TreeBH.res.DF <- TreeBH.res.DF[which(TreeBH.res[,ncol(TreeBH.res)] == 1),]
  
  # Write TreeBH for each iteration, and append at the end of the sheet controlled alleles
  # if (nrow(TreeBH.res.DF) > 0){
  #   write.xlsx(x = TreeBH.res.DF,  paste(settings$folder$HLA_Output_GLM_Iter, 'HLA_TreeBH_iter','.xlsx', sep = ''), 
  #              sheetName = paste('Iteration', as.character(iter), sep = ''), col.names = TRUE, row.names = FALSE )
  #   # write.xlsx(x = '\n', paste(settings$folder$HLA_Output_GLM_Iter, 'HLA_TreeBH_iter','.xlsx', sep = ''), 
  #   #            sheetName = paste('Iteration', as.character(iter), sep = ''), append = TRUE);
  #   for (A in as2control){
  #     write.xlsx(x = as.character(A), paste(settings$folder$HLA_Output_GLM_Iter, 'HLA_TreeBH_iter','.xlsx', sep = ''), 
  #                sheetName = paste('Iteration', as.character(iter), sep = ''), append = TRUE);
  #   }
  # }
  
  # Update
  iter <- iter + 1
  idx <- idx +1
  
}

# Format output
allele <- as.data.frame(signAlleles[[1]]); 
for (idx in 2:length(signAlleles)){
  allele <- rbind.fill(allele, signAlleles[[idx]] %>% as.data.frame())
}

# Write
for (idx in 1:length(HLA.GLM_carriers.list)){
  HLA.GLM_carriers.df <- HLA.GLM_carriers.list[[idx]]; locus <- HLA.GLM_carriers.list %>% names() %>% .[idx]
  write.xlsx(x = HLA.GLM_carriers.df, file = paste(settings$folder$HLA_Output_GLM_Iter, 'HLA_GLM_Carriers','.xlsx', sep = ''), sheetName = locus, 
             col.names = TRUE, row.names = FALSE, append = TRUE)
}
write.xlsx(x = allele, file = paste(settings$folder$HLA_Output_GLM_Iter, 'HLA_GLM_Carriers','.xlsx', sep = ''), sheetName = 'Significant_alleles', 
           col.names = TRUE, row.names = FALSE, append = TRUE)



