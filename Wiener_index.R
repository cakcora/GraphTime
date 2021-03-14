## function which compute Wiener index of graph
wiener_index = function(input_graph){
  W = as_adjacency_matrix(input_graph, type = c("both"),
                          edges = FALSE, names = TRUE,attr="weight")
  W = as.matrix(W)
  A = as_adjacency_matrix(input_graph, type = c("both"),
                          edges = FALSE, names = TRUE,attr=NULL)
  A = as.matrix(A)
  t = dim(A)[1]
  for(i in 1:t){
    for(j in 1:t){
      if(A[i,j]>1){A[i,j]=1}
    }
  }
  
  D = A
  
  for(i in 1:t){
    for(j in 1:t){
      if(D[i,j]== 0 & i!=j){D[i,j]=t}
    }
  }
  
  
  for(k in 1:t){
    for(i in 1:t){
      if(i!=k & D[i,k]<t){
        for(j in 1:t){
          if(j!=k & D[i,k]+D[k,j]<D[i,j]){
            D[i,j]=D[i,k]+D[k,j]
            W[i,j]=W[i,k]+W[k,j]
          }
        }
      }
    }
  }
  
  Wiener = sum(as.matrix(W))/2
  Wiener
}


######################################################
#Wiener index of Aragon token

library(anytime)
library(igraph)
library(TDA)
library(tidyverse)
library(depthTools) # to compute MBD
library(randomForest) # to fit RF
library(caret)

## Needed data

## Path
mainDir <-"C:/Users/o_kho/Desktop/2020_Winter/DTA/EthereumCurves_new/"
depthDir = paste0(mainDir,'depth','/')
mergeDir = paste0(mainDir,'merge','/')
modelDir = paste0(mainDir,'model','/')
graphDir = paste0(mainDir,'graph','/')
dataDir = paste0(mainDir,'data','/')
tokenDir = paste0(mainDir,'tokenPrice','/')

depthList <- list.files(depthDir,pattern = "csv.rds$")
mergeList <- list.files(mergeDir,pattern = "rds$")
modelList <- list.files(modelDir,pattern = "rds$")
graphList <- list.files(graphDir,pattern = "rds$")
dataList <- list.files(dataDir,pattern = "txt$")
tokenList <- list.files(tokenDir,pattern = "txt$")


x1 = readRDS(paste0(depthDir,'rd_aragon_B0.csv.rds'))

z1 = read_rds(paste0(mergeDir,"df_aragon_abs25.rds"))

w1 = read_rds(paste0(graphDir,"graphFeature_aragon.rds"))

######################################################

## Creating the aragon netwrok

tokenName = dataList[1] %>% str_remove(.,"network") %>% str_remove(.,"TX.txt")
fileLoc = paste0(dataDir,dataList[1])
networkDF <- read_delim(fileLoc," ",col_names = FALSE,
                        col_types =cols(.default = col_double())) %>% 
  setNames(.,c("from","to","time","value"))
timeFormat = "%Y%m%d"
networkDF$time <- networkDF$time %>% anydate() %>% format(.,timeFormat) %>% as.factor()

# transform transaction amounts
networkDF$value2 <- log(1+networkDF$value)
# OR
#networkDF$value2 <- 1/(1+networkDF$value)

# compute graph features
networkDF$time <- networkDF$time %>% anydate()
networkDF <- networkDF %>% arrange(time)

net_time = unique(networkDF$time)
tmp = networkDF %>% filter(time==net_time[1])

# create graph
rawGraph <- graph.data.frame(tmp[,c("from","to")],directed = F)
E(rawGraph)$weight = tmp$value2

#connected components of aragon network
com = decompose(rawGraph, mode = c("weak", "strong"), max.comps = NA,
                min.vertices = NA)

############################################################

## Wiener index of aragon token with reporting the computation time

start.time <- Sys.time()
W_vector = rep(0,length(net_time))

for(i in 1:length(W_vector)){
  tmp = networkDF %>% filter(time==net_time[i])
  
  # create graph
  rawGraph <- graph.data.frame(tmp[,c("from","to")],directed = F)
  E(rawGraph)$weight = tmp$value2
  
  com = decompose(rawGraph, mode = c("weak", "strong"), max.comps = NA,
                  min.vertices = NA)
  Wgraph = rep(0,length(com))
  for(j in 1:length(com)){
    Wgraph[j] = wiener_index(com[[j]])
  }
  W_vector[i] = sum(Wgraph)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken







