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


########### INITIALIZATION ############

# Set working directory
setwd("~/Documents/Alzheimer-Disease")

# Import settings
settings <- jsonlite::read_json("settings.json")

# Create comand
`%notin%` <- Negate(`%in%`)

# Import HLA calls and covariates
if(settings$dataset == "Stanford"){
  HLA.df <- read.csv(settings$file$HLA_df)
  covars.df <- read.csv(settings$file$HLA_covars_filt)
} else if (settings$dataset == "UKB"){
  HLA.df <- read.csv(settings$file$HLA_df_UKB)
  covars.df <- read.csv(settings$file$HLA_covars_UKB)
  covars.df$pheno <- covars.df$pheno +1
  
}


############ COMPUTE ALLELE/CARRIER COUNT/FREQUENCIES ############

computeACFREQ = function(data, locus, Dx){
  # Compute heterozigous, homozigous, or absence count. 
  # Compute allele frequency, count and total.
  # Compute carrier frequency, count and total.
  
  # Allele frequency 
  A1 <- locus %>% paste0('.1'); A2 <- locus %>% paste0('.2')
  alleles <- list(data[, A1], 
                  data[, A2]) %>% unlist()
  alleles.count <- table(alleles)
  alleles.freq <- alleles.count %>% prop.table() * 100
  alleles.df <- data.frame(allele = alleles.count %>% names(), 
                           alleleCount = alleles.count %>% as.vector(),
                           alleleFreq = alleles.freq %>% as.vector(), 
                           alleleTotal = nrow(data)*2)
  
  # Carrier frequency
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
  ACFREQ.df <- merge(HH.data, carrier.df, by = 'allele') %>% merge(alleles.df, by = 'allele')
  
  switch (Dx,
          'case' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2','carrierCount', 'carrierFreq', 'carrierTotal',
                                                       'alleleCount', 'alleleFreq','alleleTotal'),
                                                     rep('Case',7), sep = ''))
            return(ACFREQ.df)
          },
          'control' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2','carrierCount', 'carrierFreq', 'carrierTotal',
                                                       'alleleCount', 'alleleFreq','alleleTotal'),
                                                     rep('Control',7), sep = ''))
            return(ACFREQ.df)
          }
  )
}


############ COMPUTE CHI SQUARE ##############

computeChi2 = function(ACFREQ.df){
  
  # Initialize 
  chi2.pval.alleles <- c() ; chi2.pval.carriers <- c();
  chi2.OR.alleles <- c(); chi2.OR.carriers <- c();
  fishers.pval.alleles <- c(); fishers.pval.carriers <- c()
  fishers.OR.alleles <- c(); fishers.OR.carriers <- c()
  fishers.UI.alleles <- c(); fishers.UI.carriers <- c()
  fishers.LI.alleles <- c(); fishers.LI.carriers <- c()
  OR.carrier <- c(); OR.allele <- c()
  
  # For each allele in the locus 
  for (A in ACFREQ.df$allele){
    ## Allele frequency 
    # Create contingency table  and run chi2 test
    allele.data <- ACFREQ.df %>% filter(allele == A)
    cont.table.allele <- matrix(c(allele.data$alleleCountCase,
                                  allele.data$alleleCountControl,
                                  allele.data$alleleTotalCase - allele.data$alleleCountCase,
                                  allele.data$alleleTotalControl - allele.data$alleleCountControl),
                                nrow = 2)
    # Tests
    chi2.allele.res <- chisq.test(x = cont.table.allele)
    fishers.allele.res <- fisher.test(cont.table.allele)
    # Chi2
    chi2.pval.alleles <- c(chi2.pval.alleles, chi2.allele.res$p.value) ;
    chi2.OR.alleles <- c(chi2.OR.alleles, chi2.allele.res$statistic %>% unname())
    # Fishers' exact test 
    fishers.pval.alleles <- c(fishers.pval.alleles, fishers.allele.res$p.value);
    fishers.OR.alleles <- c(fishers.OR.alleles, fishers.allele.res$estimate %>% unname())
    fishers.LI.alleles <- c(fishers.LI.alleles, fishers.allele.res$conf.int[1])
    fishers.UI.alleles <- c(fishers.UI.alleles, fishers.allele.res$conf.int[2])
    # OR 
    OR.allele <- c(OR.allele, (cont.table.allele[1]*cont.table.allele[4])/(cont.table.allele[2]*cont.table.allele[3]))
    
    
    ## Carrier Frequency
    # Create contingency table  and run chi2 test
    cont.table.carrier <- matrix(c(allele.data$carrierCountCase,
                                   allele.data$carrierCountControl,
                                   allele.data$carrierTotalCase - allele.data$carrierCountCase,
                                   allele.data$carrierTotalControl - allele.data$carrierCountControl),
                                 nrow = 2)
    # Tests
    chi2.carrier.res <- chisq.test(x = cont.table.carrier)
    fishers.carrier.res <- fisher.test(cont.table.carrier)
    # Chi2
    chi2.pval.carriers <- c(chi2.pval.carriers, chi2.carrier.res$p.value) ; 
    chi2.OR.carriers <- c(chi2.OR.carriers, chi2.carrier.res$statistic %>% unname())
    # Fishers' exact test 
    fishers.pval.carriers <- c(fishers.pval.carriers, fishers.carrier.res$p.value);
    fishers.OR.carriers <- c(fishers.OR.carriers, fishers.carrier.res$estimate %>% unname())
    fishers.LI.carriers <- c(fishers.LI.carriers, fishers.carrier.res$conf.int[1])
    fishers.UI.carriers <- c(fishers.UI.carriers, fishers.carrier.res$conf.int[2])
    # OR 
    OR.carrier <- c(OR.carrier, (cont.table.carrier[1]*cont.table.carrier[4])/(cont.table.carrier[2]*cont.table.carrier[3]))
    
  }
  
  # Create dataframe 
  stats.data <- data.frame(allele = ACFREQ.df$allele, FishersCarrierPVAL = fishers.pval.carriers, 
                           FishersCarrierOR = fishers.OR.carriers, FishersCarrierLI = fishers.LI.carriers,
                           FishersCarrierUI = fishers.UI.carriers,
                           ChiCarrierPVAL = chi2.pval.carriers, ChiCarrierEST = chi2.OR.carriers, ORCarrier = OR.carrier,
                           FishersAllelePVAL = fishers.pval.alleles, FishersAlleleOR = fishers.OR.alleles, 
                           FishersAlleleLI = fishers.LI.alleles, FishersAlleleUI = fishers.UI.alleles,
                           ChiAllelePVAL = chi2.pval.alleles, ChiAlleleEST = chi2.OR.alleles,
                           ORAllele = OR.allele)
  
  # Remove NAs
  stats.data[is.na(stats.data)] <- 0
  
  # Return 
  return(stats.data)
  
}


############# MAIN LOOP ##############

# Get cases and controls and separate datasets
cases.ids <- covars.df$sample.id.file[which(covars.df$pheno ==2)]
controls.ids <- covars.df$sample.id.file[which(covars.df$pheno ==1)]
data.cases <- HLA.df %>% filter(sample.id.file %in% cases.ids)
data.controls <- HLA.df %>% filter(sample.id.file %in% controls.ids)

# Initialize multiple-test correction 
pvals <- c(); l1group <- c(); l2group <- c(); l2locus <- c()

# For each locus 
loci <- c("A","B","C","DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5")

idx <- 1
for (locus in loci){
  
  print(paste0("Current locus: ", locus))
  
  # Compute allele and carrier counts and frequencies
  ACFREQ.cases <- computeACFREQ(data.cases, locus, 'case');
  totalCases <-unique(ACFREQ.cases$alleleTotalCase); carrierCases <- unique(ACFREQ.cases$carrierTotalCase)
  ACFREQ.controls <- computeACFREQ(data.controls, locus, 'control');
  totalControls <-unique(ACFREQ.controls$alleleTotalControl); carrierControls <- unique(ACFREQ.controls$carrierTotalControl)
  
  # Merge and clean 
  ACFREQ.df <- merge(ACFREQ.cases[,!names(ACFREQ.cases) %in% c('A0','A1','A2')], ACFREQ.controls, by = 'allele', all = TRUE) 
  ACFREQ.df[is.na(ACFREQ.df)] <- 0
  ACFREQ.df$alleleTotalCase <- rep(totalCases, nrow(ACFREQ.df)); ACFREQ.df$alleleTotalControl <- rep(totalControls, nrow(ACFREQ.df));
  ACFREQ.df$carrierTotalCase <- rep(carrierCases, nrow(ACFREQ.df)); ACFREQ.df$carrierTotalControl <- rep(carrierControls, nrow(ACFREQ.df))
  ACFREQ.df <- ACFREQ.df %>% filter(allele != '')
  
  # Compute Chi2
  stats.data <- computeChi2(ACFREQ.df)
  stats.data[stats.data==Inf] = 'Inf'
  
  # Correction for TreeHB (see paper)
  pvals <- c(pvals, stats.data$FishersCarrierPVAL)
  l1group <- c(l1group, paste(rep(locus, nrow(stats.data)), rep('*', nrow(stats.data)), stats.data$allele, sep = ''));
  l2group <- c(l2group, rep(idx, nrow(stats.data))); l2locus <-  c(l2locus, rep(locus, nrow(stats.data)));
  idx <- idx +1
  
  # Create dataframes 
  HLA.alleles.df <- merge(stats.data[,c(1,which(grepl('Allele', colnames(stats.data))))],
                          ACFREQ.df[,c(which(grepl(paste(c('A0','A1','A2', 'allele'), collapse = '|'), colnames(ACFREQ.df))))],
                          by = 'allele')
  HLA.carriers.df <-  merge(stats.data[,c(1,which(grepl('Carrier', colnames(stats.data))))],
                            ACFREQ.df[,c(1,which(grepl(paste(c('A0','A1','A2', 'carrier'), collapse = '|'), colnames(ACFREQ.df))))],
                            by = 'allele')
  
  # Write to excel output
  if (settings$dataset == "Stanford"){
    write.xlsx(x = HLA.alleles.df, file = paste(settings$folder$HLA_Output_Chi2, 'HLA_AnalysisAlleles','.xlsx', sep = ''), sheetName = locus,
               col.names = TRUE, row.names = FALSE, append = TRUE)
    write.xlsx(x = HLA.carriers.df, file = paste(settings$folder$HLA_Output_Chi2, 'HLA_AnalysisCarriers','.xlsx', sep = ''), sheetName = locus,
               col.names = TRUE, row.names = FALSE, append = TRUE)
  } else if (settings$dataset == "UKB"){
    write.xlsx(x = HLA.alleles.df, file = paste(settings$folder$HLA_Output_Chi2, 'HLA_AnalysisAlleles_UKB','.xlsx', sep = ''), sheetName = locus,
               col.names = TRUE, row.names = FALSE, append = TRUE)
    write.xlsx(x = HLA.carriers.df, file = paste(settings$folder$HLA_Output_Chi2, 'HLA_AnalysisCarriers_UKB','.xlsx', sep = ''), sheetName = locus,
               col.names = TRUE, row.names = FALSE, append = TRUE)
  }


}


# Apply correction and apply TreeHB
groups <- cbind(as.factor(l2group),as.factor(1:length(l2group)))
resSec <- get_TreeBH_selections(pvals = pvals, groups = groups, q = rep(0.05, 2))
allelesSec <- l1group[which(resSec[,2] == 1)]
locusSec <- l2locus[which(resSec[,1] == 1)]
totalSec <-l1group[which(resSec[,2] == 1 & resSec[,1] == 1)]

# dataframe results 
HBTree <- data.frame(pval = pvals, locus=l2locus, allele = l1group); HBTree <- HBTree[which(resSec[,2] == 1),]

# Write results 
write.csv(x = HBTree, file = paste(settings$folder$HLA_Output_Chi2, 'HLA_HBTree_Carriers','.xlsx', sep = ''))


