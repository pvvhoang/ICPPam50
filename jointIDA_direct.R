jointIDA_direct <- function(datacsv, cause, effect, method = c("min","max","median"), pcmethod = "stable", alpha, num.cores = 1, mem.efficient = FALSE, technique = c("RRC","MCD")){
  if(is.character(datacsv)){
    data <- read.csv(datacsv)
    #data=data[,-1] # if the dataset have the sample names column, otherwise comment this out.
  } else {
    data <- datacsv #Assume there is no samplenames column and this is a data.frame.
  }      				#To allow both .csv data input or a matrix in R. This will help the IDAbootstrap, as IDA can run on sampling matrices.
  data <- scale(data) #standardise the data
  #print(data[1:5,])
  allnames <- colnames(data)
  causenames <- allnames[cause]
  effectnames <- allnames[effect]
  
  multiset <- character(0)
  result <- matrix(nrow=length(effect), ncol = length(cause))
  suffStat <- list(C=cor(data), n = nrow(data))
  indepTest <- gaussCItest
  
  start_total_jointida <- proc.time()
  print(" Start runing pc...")
  if ((pcmethod == "stable") || (pcmethod == "original")){
    pcFit <- pc_stable(suffStat, indepTest, p = ncol(data), alpha = alpha, skel.method = pcmethod)
  }else {
    pcFit <- pc_parallel(suffStat, indepTest = gaussCItest, p = ncol(data), skel.method = pcmethod, alpha = alpha, num.cores = num.cores, mem.efficient = mem.efficient)
  }
  print("Finished. Now we calculate the joint causal effects...")
  #pcFit<-pc(suffStat, indepTest=gaussCItest, p=ncol(data),alpha=alpha)
  #return(jointIda(cause,effect,cov(data),pcFit@graph,technique=technique))
  
  # get number of cores to run
  cl<-makeCluster(num.cores)
  registerDoParallel(cl)
  
  temp = rep(NA, length(cause))
  sink("log.txt")
  result.tmp <- foreach(k = 1:length(effect)) %dopar% {
    ##=====================
    #sink("log.txt", append=TRUE)
    #cat("calculate the effects on the ", k, "th mRNA", "\n")
    ##=====================
    caef <- pcalg::jointIda(cause, effect[k], cov(data), pcFit@graph, technique=technique, type="cpdag")
    ##====================
    #cat("Finished jointIDA, now get the unique causal effect", "\n")
    ##====================
    caefabs <- abs(caef)
    
    mat.tmp = temp
    
    for(l in 1:length(cause)){
      
      if(method=="min"||method=="max"){
        if(method=="min"){
          index <- which(caefabs==min(caefabs[l, ], na.rm = TRUE), arr.ind = TRUE)
        }else{
          index <- which(caefabs==max(caefabs[l, ], na.rm = TRUE), arr.ind = TRUE)
        }
        # cat("index",index,"\n")
        # cat("index[1,2]",index[1,2],"\n")
        
        pos <- index[1, 2]
        mat.tmp[l] <- caef[l, pos]
      }else if(method=="median"){
        mat.tmp[l] <- median(caef[l, ], na.rm = TRUE)
      }
    }
    mat.tmp
  }
  
  # shut down the workers
  stopCluster(cl)
  stopImplicitCluster()
  
  # create a matrix from a list
  result = matrix(unlist(result.tmp), ncol = length(cause), byrow = T)
  
  total_t_jointida <- proc.time()-start_total_jointida
  #cat('Total Time jointida=', total_t_jointida[3], '\n', sep=" ")
  colnames(result) <- causenames
  rownames(result) <- effectnames
  return(result)
  
}#jointIDA
