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
library(tictoc)
options(dplyr.summarise.inform = FALSE)
#library(fda.usc) # to compute modal depth

#############
# MAIN BODY #
#############

mainDir <- here::here("Ethereum")

folders <- c('data','tokenPrice','betti','pl','pi','graph','depth','merge','model','pd','results')
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

networkDF <- allTokens %>% filter(time %in% selectedDays)
networkDF$time <- ymd(networkDF$time)

# normalize "value" 
networkDF <- networkDF %>% group_by(name) %>% 
  mutate(value=log(value)/log(max(value))) %>% .[,1:4]          
             
# compute graph features
featureGraph() 

# compute PDs
topRank <- 150
computePD(topRank)

# compute Betti sequences
scale_seq=seq(0,1.01,length.out = 100) # sequences of scale values
tic()
computeBetti(scale_seq)
toc()

# compute Persistence Landscapes (PL)
scale_seq=seq(0,1.01,by=0.01) # sequences of scale values
tic()
computePL(scale_seq)
toc()

# compute Persistence Images (PI)
res <- 10 # PI resolution
tic()
computePI(res)
toc()

# compute rolling depth
topoSignature <- 'betti' # choose from 'betti', pl','pi'
rollDepth(topoSignature)

# merge data
threshold <- 0.05 # day is anomalous if price return changes by more than this threshold in the next h days. 
dataMerge(threshold)

# fit RF model
set.seed(1)
RFmodel(threshold,topoSignature,repNum=10)

# gather all results
gatherResults(topoSignature)






