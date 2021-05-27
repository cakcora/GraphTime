library(igraph)
library(TDA)
library(tidyverse)
library(depthTools) # to compute MBD
library(randomForest) # to fit RF
library(caret)
library(fda.usc)
library(pROC)
library(lubridate)
options(dplyr.summarise.inform = FALSE)
#library(fda.usc) # to compute modal depth

#############
# MAIN BODY #
#############
mainDir <- "Ethereum"

folders <- c('data','tokenPrice','betti','graph','depth','merge','model','pd','results')
for (f in folders){
  assign(paste0(f,'Dir'),paste0(mainDir,f,'/'))
  if (!file.exists(paste0(mainDir,f))) 
    dir.create(paste0(mainDir,f,'/'),showWarnings = FALSE)
}

# load all functions
source(paste0(mainDir,'functions.R'))

# filtration type
filtration <- 'sublevel'

# load all token data
allTokens <- readRDS(file=paste0(dataDir,'allTokens.rds'))

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

# compute PDs
topRank <- 250
computePD(topRank)

# compute Betti sequences
scale_seq=c(seq(0,1,by=0.01),2) # sequences of scale values
computeBetti()

# compute rolling depth
rollDepth()

# merge data
threshold <- 0.05 # day is anomalous if price return changes by more than this threshold in the next h days. 
dataMerge(threshold)

# fit RF model
RFmodel(threshold,repNum=5)

# gather all results
gatherResults()






