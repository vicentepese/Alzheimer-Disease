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
covars.df <- read.csv(settings$file$HLA_covars_filt)
covars.df$pheno <- covars.df$pheno -1
covars.df$sex <- covars.df$sex -1

# Add phenotype to HLA dataframee 
HLA.df <- merge(HLA.df, covars.df[,c("sample.id.file", "pheno")], by = "sample.id.file")

# Import amino acid alignment 
AA_alignment <- read.table(settings$file$AA_alignment, header = TRUE, sep = ',')

# Remove negative positions 
full_prot = FALSE
if (!full_prot){
  AA_alignment[,3] <- AA_alignment %>% apply(MARGIN = 1, function(x) x[3] %>% substr(start = 30, stop = nchar(x[3])))
}


########## EXTRA FUNCTIONS ############

carrierCount = function(settings, data){
  
  # Get loci, number of cases, number of controls
  loci <- colnames(data)[which(grepl('HLA', colnames(data)))] %>%
    lapply(function(x) x %>% strsplit('_') %>% unlist() %>% .[1] %>% strsplit('HLA') %>% unlist() %>% .[2]) %>% unlist() %>% unique()
  c(Ncases, Ncontrols) %<-% c(data$Dx %>% table() %>% .['1'], data$Dx %>% table() %>% .['0'])
  
  # Initialize loop 
  alleles <- c(); locusVec <- c(); cases.Count <- c(); controls.Count <- c()
  
  # For each locus, count carrier in cases and controls
  for (locus in loci){
    
    # Get alleles 
    c(A1, A2) %<-% c(paste('HLA', locus, '_A1', sep = ''), paste('HLA', locus, '_A2', sep = ''))
    data.locus <- data[, c('GWASID', A1, A2, 'Dx')]
    
    # Get carriers and count 
    carriers.locus <- c(data.locus[,A1], data.locus[,A2]) %>% unique()
    carriers.locus <- carriers.locus[which(carriers.locus != "")]
    for(carrier in carriers.locus){
      data.carrier <- data.locus %>% filter(get(A1) == carrier | get(A2) == carrier)
      c(C.cases, C.controls ) %<-% c(table(data.carrier$Dx)['1'], table(data.carrier$Dx)['0']);
      if(C.cases %>% is.na()){C.cases <- 0}; if (C.controls %>% is.na()) {C.controls <- 0}
      c(F.cases, F.controls) %<-% c(C.cases/Ncases, C.controls/Ncontrols)
      
      # If frequency is smaller than a threshold, next. Else add to vector 
      if (is.na(C.cases) | is.na(C.controls)){
        next 
      } else if (F.cases < settings$min_carrierFreq | F.controls < settings$min_carrierFreq){
        next
      } else{
        locusVec <- c(locusVec, locus)
        A <- paste( locus, '*', carrier, ':01', sep = ''); alleles <- c(alleles, A)
        cases.Count <- c(cases.Count, C.cases); controls.Count <- c(controls.Count, C.controls)
      }
    }
  }
  
  # Create data frame 
  counts.data <- data.frame(locus = locusVec, allele = alleles, Ncases = cases.Count, Ncontrols = controls.Count)
  
  # Return 
  return (counts.data)
  
}


countAA2control = function(settings, data, AA_alignment){
  
  # Get locus, position and AA to control 
  AA2control <- settings$AA2control %>% strsplit('_') %>% unlist()
  L2control <- AA2control[1]; posControl <- AA2control[2]; AAcontrol <- AA2control[3]
  c(L2controlA1, L2controlA2) %<-% c(paste0('HLA',L2control,'_A1'), paste0('HLA',L2control,'_A2'))
  
  # Get sequences
  AA2control_locus <- AA_alignment %>% filter(locus == AA2control[1])
  
  # Get alleles 
  allelesAA.OG <- AA2control_locus %>% apply(MARGIN = 1, function(x, posControl, AAcontrol)
    if (substr(x[3],posControl,posControl) == AAcontrol) {return(x[2])}, posControl, AAcontrol) %>% unlist()
  allelesAA <- allelesAA.OG %>% lapply(function(x) strsplit(x, split='\\*') %>% unlist() %>% 
                                         .[2] %>% strsplit(split=':') %>% unlist() %>% .[1:2] %>% paste(collapse=':')) %>% unlist()
  
  # Get data and count cases and controls
  data.AA <- data %>% filter(get(L2controlA1) %in% allelesAA| get(L2controlA2) %in% allelesAA)
  c(AA2cntrl.cases, AA2cntrl.controls) %<-% c(data.AA$Dx %>% table() %>% .['1'], data.AA$Dx %>% table() %>% .['0']);
  
  return(c(AA2cntrl.cases, AA2cntrl.controls))
}


########## AMINO ACID ANALYSIS ##########

# Get loci, number of cases and controls
loci <- unique(AA_alignment$locus)
c(Ncases, Ncontrols) %<-% c(HLA.df$pheno %>% table %>%.['1'], HLA.df$pheno %>% table %>%.['0'])

# Get AA2control counts
if (settings$AA2control %>% length() > 0){
  c(AA2cntrl.cases, AA2cntrl.controls) %<-% countAA2control(settings,HLA.df, AA_alignment)
} else{
  AA2cntrl.cases <- 0; AA2cntrl.controls <- 0
}

# Initialize loop
AA.total <- c(); pval <- c(); OR <- c(); pos.total <- c(); locus <- c();
cases.Count <- c(); controls.Count <- c(); alleles <- c()

# For each locus 
for (L in loci){
  
  # Get locus ID
  c(A1, A2) %<-% c(paste0(L,'.1'), paste0(L,'.2'))
  
  # Get alignment subset, and counts subset
  AA_locus <- AA_alignment %>% filter(locus == L)
  
  # Get max sequence length 
  maxLen <- AA_locus$sequence %>% lapply(nchar) %>% unlist() %>% max()
  
  # Get unique AAs 
  for (pos in 1:maxLen){
    
    # Get unique AAs
    AAs <- AA_locus$sequence %>% lapply(function(x, pos) substr(x, pos, pos), pos) %>% unlist() %>% unique()
    
    # If more than one, count 
    if (length(AAs) >1){
      
      for (AA in AAs){
        # Get alleles with AAs
        allelesAA.OG <- AA_locus %>% apply(MARGIN = 1, function(x, pos, AA) if (substr(x[3],pos,pos) == AA) {return(x[2])}, pos, AA) %>% unlist()
        allelesAA <- allelesAA.OG %>% lapply(function(x) strsplit(x, split='\\*') %>% unlist() %>% 
                                               .[2] %>% strsplit(split=':') %>% unlist() %>% .[1:2] %>% paste(collapse=':')) %>% unlist() 
        
        # Get data with such alleles
        data.AA <- HLA.df %>% filter(get(A1) %in% allelesAA| get(A2) %in% allelesAA)
        
        # Get counts 
        c(AA.cases, AA.controls) %<-% c(data.AA$pheno %>% table() %>% .['1'], data.AA$pheno %>% table() %>% .['0']);
        if(is.na(AA.cases)) {AA.cases <- 0}; if(is.na(AA.controls)){AA.controls <- 0}
        
        # Contingency table 
        cont.table <- matrix(c(AA.cases, AA.controls, 
                               Ncases - AA.cases, Ncontrols - AA.controls), nrow = 2)
        
        # Compute chiSq test
        chi.sq.res <- fisher.test(cont.table)
        
        # Append to vectors
        AA.total <- c(AA.total, AA); pval <- c(pval, chi.sq.res$p.value);
        OR <- c(OR, (cont.table[1]*cont.table[4])/(cont.table[2]*cont.table[3])); 
        pos.total <- c(pos.total, pos); locus <- c(locus, L);
        cases.Count <- c(cases.Count, AA.cases); controls.Count <- c(controls.Count, AA.controls)
        alleles <- c(alleles, paste(allelesAA.OG, collapse = ', '))
      }
    }
  }
}

# Create dataframe
AA.analysis.results <- data.frame(locus= locus, AA = AA.total, pos = pos.total, Ncases = cases.Count, Ncontrol = controls.Count,
                                  pval = pval, OR = OR, alleles = alleles)
# Save dataframe 
write.csv(AA.analysis.results, file = paste0(settings$folder$HLA_Output_AA, "HLA_AA.csv"), row.names = FALSE)


