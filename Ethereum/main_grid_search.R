library(igraph)
library(TDA)
library(tidyverse)
library(depthTools) # to compute MBD
library(randomForest) # to fit RF
library(caret)
library(fda.usc)
library(pROC)
library(lubridate)
library(here)
options(dplyr.summarise.inform = FALSE)
#library(fda.usc) # to compute modal depth

#############
# MAIN BODY #
#############

#mainDir <- "/Users/iumar/Box Sync/Research/EthereumCurves-new/"
mainDir <- here::here("Ethereum")

folders <- c('data','tokenPrice','betti','graph','depth','merge','model','pd','results')
for (f in folders){
  assign(paste0(f,'Dir'),file.path(mainDir,f))
  if (!file.exists(file.path(mainDir,f))) 
    dir.create(file.path(mainDir,f,'/'),showWarnings = FALSE)
}

# load all functions
source(file.path(mainDir,'functions.R'))

# filtration type
filtration <- 'sublevel'

# load all token data
allTokens <- readRDS(file=file.path(dataDir,'allTokens.rds'))

# select days on which at least 5 tokens are traded
selectedDays <- allTokens %>% group_by(time) %>%
  summarise(n=length(unique(name))) %>%
  filter(n>=5) %>% .[,1] %>% unlist()
allTokens %>% filter(time %in% selectedDays) -> networkDF
networkDF$time <- ymd(networkDF$time)

# normalize "value" 
minMaxValues <- networkDF %>% group_by(name) %>% 
  summarise(minValue=min(value),maxValue=max(value)) 
networkDF <- inner_join(networkDF,minMaxValues,by="name") %>% 
  mutate(value=(log10(maxValue)-log10(value))/(log10(maxValue)-log10(minValue))) %>% .[,1:4]

# compute graph features
featureGraph() 

resultsSummary <- tibble()
# grid search
for (topRank in c(150,200,250,300))
  for (threshold in c(0.04,0.05,0.06)){
      # compute PDs
      computePD(topRank)
      
      # compute Betti sequences
      scale_seq=c(seq(0,1,by=0.01),2) # sequences of scale values
      computeBetti()
      
      # compute rolling depth
      rollDepth()
      
      # merge data
      # day is anomalous if price return changes by more than this threshold in the next h days. 
      dataMerge(threshold)
      
      # fit RF model
      RFmodel(threshold,repNum=5)
      
      # gather all results
      gatherResults()
      
      resultsCurrent <- read_csv(file.path(resultsDir,paste0('AUC_',filtration,'.csv')))[,-1]
      resultsCurrent <- resultsCurrent %>% mutate(Comb=paste0(topRank,'+',threshold))
      resultsSummary <- rbind(resultsSummary,resultsCurrent)
}

saveRDS(resultsSummary,file=file.path(resultsDir,'resultsSummary.rds'))






