########
featureGraph = function(){
  
  # calculate graph features per day
  featurePerDay = function(inputDay){
    
    periodDF = networkDF %>% filter(time==inputDay)
    
    # construct graph
    G <- graph.data.frame(periodDF[,c("from","to")],directed = F)
    
    # the return result
    tibble(
      day = inputDay,
      vertexNum = vcount(G),
      edgeNum = ecount(G),
      clusterCoef = transitivity(G)
    )
  }
  # body of featureGraph()
  periodList = networkDF$time %>% unique()
  res = map_dfr(periodList,featurePerDay)
  write_rds(res,paste0(graphDir,"graphFeature_allTokens.rds"))
  
}

#################
computePD = function(topRank){
  
  selectTop <- function(dat,topRank){ 
    frqFrom <- table(dat$from) %>% sort(decreasing = T) %>% names %>% .[1:topRank]
    dat <- subset(dat,from %in% frqFrom)
    frqTo <- table(dat$to) %>% sort(decreasing = T) %>% names %>% setdiff(frqFrom) %>% .[1:topRank]
    topList <- union(frqFrom,frqTo)
    dat <- subset(dat,(to %in% topList) & (from %in% topList))
  } 
  
  # body of computePD()
  periodList = networkDF$time %>% unique() 
  nDays <- length(periodList)
  
  for (i in 1:nDays){
    periodDF <- networkDF %>% filter(time==periodList[i])
    periodDF <- selectTop(periodDF,topRank)
    periodDF <- periodDF %>% group_by(from,to) %>% summarise(value=mean(value)) 
    
    G <- graph.data.frame(periodDF[,c("from","to")],directed = F)
    #print(vcount(G))
    
    # compute all simplices of dimension 0,1,2
    cmplx <- cliques(G,min=1,max=3) # 
    
    # sublevel filtration
    vertValues<-rep(1,length(V(G)))
    names(vertValues)<-unique(c(periodDF$from,periodDF$to))
    
    # compute values of the feature function (here sum()) for each vertex
    valuesFrom <- tapply(periodDF$value,periodDF$from,mean) 
    vertValues[names(valuesFrom)] <- valuesFrom
    
    #valuesTo <- tapply(periodDF$value,periodDF$to,mean) 
    #vertValues[names(valuesTo)] <- vertValues[names(valuesTo)] + valuesTo
    
    Flt <- funFiltration(FUNvalues = vertValues, 
                         cmplx = cmplx, 
                         sublevel = (filtration=='sublevel'))
    PD <- filtrationDiag(filtration = Flt,maxdimension = 1,library = 'Dionysus',location = T)$diagram
    
    # replacing infinity in PD with finite value of 2
    PD[PD[,3]==Inf,3]=2
    
    saveRDS(PD,file=paste0(pdDir,'PD_',filtration,'_',periodList[i],'.rds'))
    
    if (i %% 30 == 0) print(paste(i,'out of',nDays,'days processed'))
    
  }
}

#################
computeBetti = function(){
  
  extractBetti=function(D){
    x <- D[,1] # birth times
    y <- D[,2] # death times 
    delta <- diff(scale_seq) # deltas 
    n <- length(delta)
    betti<-numeric(length = n)
    
    for (k in 1:n){
      b <- pmin(scale_seq[k+1],y)-pmax(scale_seq[k],x)
      betti[k] <- sum(pmax(0,b))/delta[k] 
    }
    betti # returned object
  }
  
  # body of computeBetti()
  periodList = networkDF$time %>% unique() 
  nDays <- length(periodList)

  filt_len<-length(scale_seq)-1
  B0 = B1 = matrix(0,nrow = nDays,ncol = filt_len)
  
  for (i in 1:nDays){
    # read in PD  
    PD <- readRDS(file=paste0(pdDir,'PD_',filtration,'_',periodList[i],'.rds')) 
    
    # compute Betti sequences
    B0[i,]<-extractBetti(PD[PD[,1]==0,2:3,drop=F])
    B1[i,]<-extractBetti(PD[PD[,1]==1,2:3,drop=F])
    
    if (i %% 30 == 0) print(paste(i,'out of',nDays,'days processed'))
    
  }
  
  colnames(B0)<-rep(scale_seq[-1],times=1)
  rowNames<-data.frame(Time=paste0('B0',periodList))
  saveRDS(cbind(rowNames,B0),file=paste0(bettiDir,'allTokens_B0.rds'))
  
  colnames(B1)<-rep(scale_seq[-1],times=1)
  rowNames<-data.frame(Time=paste0('B1',periodList))
  saveRDS(cbind(rowNames,B1),file=paste0(bettiDir,'allTokens_B1.rds'))
  
}

########
rollDepth = function(){
  
  rollDepthPerFile = function(f,rollSize=7){
    #message(f)
    # read in Betti sequences
    df <- readRDS(paste0(bettiDir,f))
    # normalize Betti sequences
    df[,-1] <- df[,-1]/apply(df[,-1],1,function(x) max(1,max(x)))
    rollDepth = NULL
    
    for (i in rollSize:nrow(df)) {
      betti_roll = df[(i-rollSize+1):i,-1]  #df[(i-(rollSize-1)):i,-1]  
      # MBD depth
      rollDepth[i] = MBD(df[(i-1):i,-1],betti_roll,plotting = F)$MBD[2]
    }
    res <- tibble(
      day = str_sub(df$Time,-10) %>% ymd,
      rollDepth = rollDepth,
      type = str_sub(df$Time,0,-11)
    )
    write_rds(res,paste0(depthDir,"rd_",f))
  }
  
  # body of rollDepth()
  fileList = list.files(bettiDir,paste0("(allTokens).*")) #list.files(bettiDir)
  invisible(map_dfr(fileList,rollDepthPerFile))
}

################# 
dataMerge <- function(threshold){
  addGraphFeature = function(keyword="graphFeature"){
    inputFile = list.files(graphDir,paste0("(allTokens).*"))
    df = read_rds(paste0(graphDir,inputFile))
  }
  #
  addRollDepth = function(keyword="rd"){
    # 
    inputFiles = list.files(depthDir,paste0("(allTokens).*"))
    df = map_dfr(paste0(depthDir,inputFiles),function(inputFile){
      read_rds(inputFile)
    })
    
    df = df %>% spread(type,rollDepth)
    df = df %>% arrange(day)
  }
  #
  addPrice = function(){
    inputFile = list.files(tokenPriceDir,paste0("Eth_price.*"))
    df = read_csv(paste0(tokenPriceDir,inputFile),col_types = cols(Date = col_date(format = "%Y-%m-%d")))
    df = rename(df,Open=`24h Open (USD)`)
  }
  #
  procPrice = function(){
    df = priceEth %>% dplyr::select(Date,Open)
    price = df$Open
    priceDiff = diff(price)
    priceReturn = priceDiff / price[1:length(priceDiff)]
    df = df %>% mutate(openNorm = Open/max(Open),priceReturn = c(NA,priceReturn))
  }
  
  # flag under all horizon
  flagPeriod = function(flagRange = 1:7){
    # flag under one horizon
    flagDaysFun = function(inputVector,ndays=3){
      # flag 1:7 days
      absReturn = abs(inputVector)
      l = length(inputVector)
      flagDays = vector(length = l)
      for (i in 1:(l-1)) {
        maxChange = max(absReturn[(i+1):(i+ndays)], na.rm=TRUE)
        # assign as factor
        flagDays[i] = ifelse(maxChange>=threshold,
                             TRUE,
                             FALSE)
      }
      return(flagDays)
    }
    df = price
    for (x in flagRange) {
      df = df%>% mutate(!!(paste0("flag",x)):=flagDaysFun(priceReturn,x))
    }
    df
  }
  #
  mergeData = function(){
    df = inner_join(ghFeature,rollDepth,by="day")
    df = inner_join(df,price,by=c("day"="Date"))
  }
  # output the final data set
  outputRDS = function(){
    fileName = paste0(mergeDir,"df_allTokens_abs",100*threshold,".rds")
    saveRDS(df,fileName)
  }
  # body of dataMerge()
  ghFeature = addGraphFeature()
  rollDepth = addRollDepth()
  priceEth = addPrice()
  price = procPrice()
  price = flagPeriod()
  df = mergeData()
  outputRDS()
}

###############
RFmodel <- function(threshold,repNum){
  
  RFpred = function(repNum = 100){
    
    # iterate from flag1 to flag7
    res = map_dfr(flags,function(flagi){
      set.seed(1)
      #message(flagi)
      flag_resRuns = map_dfr(1:repNum,function(repID){
        # iterate all formula
        flag_res = map2_dfr(inputfmls,names(inputfmls),function(fmlstr,fmlname){
          
          rf = randomForest::randomForest(as.formula(paste0(flagi,fmlstr)),
                                          dftrain,importance = TRUE, na.action = na.omit)
          
          cf = caret::confusionMatrix(
            predict(rf,dftest),
            dftest %>% pull(flagi),
            positive = "TRUE"
          )
          cfmatrix = rbind(as.matrix(cf,what="overall"),
                           as.matrix(cf,what="classes"))
          cftb = cfmatrix %>% as.data.frame()  %>% rownames_to_column() %>% as_tibble()
          names(cftb) = c("measure","value")
          cftb[nrow(cftb)+1,"measure"] <- "AUC"
          cftb[nrow(cftb),"value"] <- pROC::auc(unlist(dftest[,flagi]),predict(rf,dftest,type = 'prob')[,"TRUE"],quiet=T) %>% as.numeric()
          cftb$horizon = flagi
          cftb$modelType = fmlname
          return(cftb)
        })
        flag_res$id = repID
        return(flag_res)
        
      })
      
    })
    res
    
  }
  outputRDS = function(prefix = "rf_"){
    saveRDS(res,paste0(modelDir,prefix,"allTokens.rds"))
  }
  
  # body of RFmodel()
  inputfmls = c(M1="~ vertexNum + edgeNum + clusterCoef + openNorm",
                M2="~ vertexNum + edgeNum + clusterCoef + openNorm + B0",
                M3="~ vertexNum + edgeNum + clusterCoef + openNorm + B0 + B1")
  file = paste0(mergeDir,"df_allTokens_abs",100*threshold,".rds")
  df = read_rds(file)
  flags = df %>% dplyr::select(contains("flag")) %>% colnames()
  df = df %>% mutate_if(is.logical,factor,levels=c("TRUE","FALSE"))
  df$clusterCoef = replace(df$clusterCoef,is.nan(df$clusterCoef),0)
  dateRange = range(df$day)
  spDate = dateRange[1] + diff(dateRange)/3*2
  dftrain = df[df$day<spDate,]
  dftest = df[!df$day<spDate,]
  
  # Model fit 
  res = RFpred(repNum)
  outputRDS()
  
}

#################
gatherResults <- function(){
  files = list.files(modelDir)
  
  gain = function(x,ben){
    return(1-ben/x)
  }
  
  for (measureStr in c('Accuracy','Sensitivity','Precision','AUC')){
    df = map_dfr(files,function(f){
      read_rds(paste0(modelDir,f)) %>% filter(measure == measureStr)
    })

    dfs = df %>% group_by(horizon,modelType) %>% summarise(avg = mean(value,na.rm = TRUE)) %>% ungroup()
    
    dfs = dfs %>% spread(modelType,avg)
    
    write.csv(dfs,paste0(resultsDir,measureStr,"_",filtration,".csv"))
    # dfs
    #  formula: 1 - M1/M4
    dfplot = map_dfc(paste0("M",2:3),function(x){
      res_tb = tibble(gain(dfs %>% pull(x),
                           dfs$M1))
      names(res_tb) = paste0(x,"_gain")
      res_tb
    })
    dfplot = bind_cols(dfs %>% dplyr::select(horizon),dfplot)
    dfplot = dfplot %>% gather(M2_gain:M3_gain,key="type",value='value')
    dfplot2 = dfplot %>% mutate(horizon=str_remove(horizon,"flag"),type=str_remove(type,"_gain"))
    
    g = ggplot(data = dfplot2,aes(x=horizon,y=value,fill = type)) + geom_bar(colour="black", position="dodge",stat = "identity") + 
      scale_fill_brewer() + xlab("Prediction horizon (unit: day)") + ylab(paste("Gain in:",measureStr)) + 
      theme_minimal() + theme(text = element_text(size=22),legend.position="bottom")
    
    ggsave(paste0(resultsDir,measureStr,"_",filtration,".png"),g,width = 8,height = 6)
  }
  
}