library(igraph)
library(TDA)
library(tidyverse)
library(depthTools) # to compute MBD
library(randomForest) # to fit RF
library(caret)
library(fda.usc)
library(ModelMetrics)
library(lubridate)
options(dplyr.summarise.inform = FALSE)
#library(fda.usc) # to compute modal depth

#############
# MAIN BODY #
#############
mainDir <- "/Users/iumar/Box Sync/Research/EthereumCurves-new/"

folders <- c('data','tokenPrice','betti','graph','depth','merge','model','results')
for (f in folders){
  assign(paste0(f,'Dir'),paste0(mainDir,f,'/'))
  if (!file.exists(paste0(mainDir,f))) 
    dir.create(paste0(mainDir,f,'/'),showWarnings = FALSE)
}

# load all functions
source(paste0(mainDir,'functions.R'))

# sequences of scale parameter values
scale_seq=seq(0,1,by=0.01)  

# filtration type
filtration <- factor('sublevel',levels = c('sublevel','superlevel')) # or 'superlevel'

# load all token data
allTokens <- readRDS(file=paste0(dataDir,'allTokens.rds'))

# select days on which at least 5 tokens are traded
selectedDays <- allTokens %>% group_by(time) %>% 
                summarise(n=length(unique(name))) %>% 
                filter(n>=5) %>% .[,1] %>% unlist()
allTokens %>% filter(time %in% selectedDays) -> networkDF
networkDF$time <- ymd(networkDF$time) 

# drop extreme observations
dropExtreme()

# normalize "value" to lie between 0.1 and 1
minMaxValues <- networkDF %>% group_by(name) %>% 
                summarise(minValue=min(value),maxValue=max(value)) 
networkDF <- inner_join(networkDF,minMaxValues,by="name") %>% 
             mutate(value=1/(1+9*(log10(value)-log10(minValue))/(log10(maxValue)-log10(minValue)))) %>% .[,1:4]

# compute graph features
featureGraph() 

# compute Betti sequences
computeBetti(toprankFrom = 250,toprankTo = 250)

# compute rolling depth
rollDepth()

# merge data
threshold <- 0.15 # day is anomalous if price return changes by more than this threshold in the next h days. 
dataMerge(threshold)

# fit RF model
RFmodel(threshold,repNum=10)

# gather all results
gatherResults()






