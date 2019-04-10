#=====================================
# This function is to add new_death and death_event to clinical
# Input parameters:
#   clinical: data of clinical
# The output data is the clinical with new_death and death_event
#=====================================
processClinical <- function(clinical){
  clinical <- cbind(clinical, new_death = 0, death_event = 0)
  clinical[, "new_death"] <- ifelse(is.na(as.numeric(as.character(clinical[, "days_to_death"]))),
                                  as.numeric(as.character(clinical[, "days_to_last_follow_up"])),
                                  as.numeric(as.character(clinical[, "days_to_death"])))
  clinical[, "death_event"] <- ifelse(clinical[, "vital_status"] == "alive", 0, 1)
  return(clinical)
}

#' A data getting function
#'
#' This function allows you to get data of miRs and mRs.
#' @param matchedData Matched data including miRs, mRNAs and clinical.
#' @param cutoff_miR Value to cut off miRs. Defaults to 0.05.
#' @param cutoff_mR value to cut off mRs. Defaults to 0.05.
#' @return A list containing samples of miRs and mRs. The list of elements are:
#' \item{d} {The data with rows being samples and columns being miRs and mRs.}
#' \item{miRs} {The names of miRs.}
#' \item{mRs} {The names of mRs.}
#' @keywords 
#' @export
#' @examples
#' load(paste(directoryPath, "01BRCA/BRCA_matchedData.RData", sep=""))
#' cutoff_miR <- 0.05
#' cutoff_mR <- 0.0005
#' l <- getData(matchedData, cutoff_miR, cutoff_mR)
getData <- function(matchedData, cutoff_miR = 0.05, cutoff_mR = 0.05) {
  
  # Add new_death and death_event to clinical
  clinical <- processClinical(matchedData$clinical)
  
  # Identify significant miRNAs by using function FSbyCox in CancerSubtypes package
  miRNAsData = FSbyCox(t(matchedData$miRs), as.numeric(clinical[, "new_death"]),
                       as.numeric(clinical[, "death_event"]), cutoff=cutoff_miR)
  miRNAsData <- t(miRNAsData)
  miRs <- colnames(miRNAsData)
  
  # Identify significant mRNAs by using function FSbyCox in CancerSubtypes package
  mRNAsData = FSbyCox(t(matchedData$mRNAs), as.numeric(clinical[, "new_death"]),
                      as.numeric(clinical[, "death_event"]), cutoff=cutoff_mR)
  mRNAsData <- t(mRNAsData)
  # Remove duplicated data
  mRNAsData <- mRNAsData[,-which(duplicated(colnames(mRNAsData)))]
  
  mRs <- colnames(mRNAsData)
  
  # Combine data
  d <- cbind(miRNAsData, mRNAsData)
  
  l = list(d = d, miRs = miRs, mRs = mRs)
  
  return(l)
}

#' A data getting function
#'
#' This function allows you to get data of miRs and mRs.
#' @param matchedData Matched data including miRs, mRNAs.
#' @param topk_miR Value to cut off miRs.
#' @param topk_mR value to cut off mRs.
#' @return A list containing samples of miRs and mRs. The list of elements are:
#' \item{d} {The data with rows being samples and columns being miRs and mRs.}
#' \item{miRs} {The names of miRs.}
#' \item{mRs} {The names of mRs.}
#' @keywords 
#' @export
#' @examples
#' load(paste(directoryPath, "01BRCA/BRCA_matchedData.RData", sep=""))
#' topk_miR <- 30
#' topk_mR <- 1500
#' l <- getData(matchedData, topk_miR, topk_mR)
getDatabyMAD <- function(matchedData, topk_miR, topk_mR) {
  
  # Identify significant miRNAs by using function FSbyMAD in CancerSubtypes package
  miRNAsData = FSbyMAD(t(matchedData$miRs), value = topk_miR)
  miRNAsData <- t(miRNAsData)
  miRs <- colnames(miRNAsData)
  
  # Identify significant mRNAs by using function FSbyMAD in CancerSubtypes package
  mRNAsData = FSbyMAD(t(matchedData$mRNAs), value=topk_mR)
  mRNAsData <- t(mRNAsData)
  # Remove duplicated data
  mRNAsData <- mRNAsData[,!duplicated(colnames(mRNAsData))]
  
  mRs <- colnames(mRNAsData)
  
  # Combine data
  d <- cbind(miRNAsData, mRNAsData)
  
  l = list(d = d, miRs = miRs, mRs = mRs)
  
  return(l)
}

#' A function for getting mapping information of miRs and mRNAs
#'
#' This function allows you to convert miR names to version 21 and get mapping information of miRs and mRs.
#' @param l A list including samples of miRs and mRs as well as miR names and mR names.
#' @param mappingTable Mapping table of miRs and mRs.
#' @return A list containing samples of miRs and mRs as well as mapping list. The list of elements are:
#' \item{d} {The data with rows being samples and columns being miRs and mRs.}
#' \item{miRs} {The names of miRs.}
#' \item{mRs} {The names of mRs.}
#' \item{mapingList} {The mapping list of miRs and mRs.}
#' @keywords 
#' @export
#' @examples
#' load(paste(directoryPath, "01BRCA/BRCA_matchedData.RData", sep=""))
#' cutoff_miR <- 0.05
#' cutoff_mR <- 0.0005
#' l <- getData(matchedData, cutoff_miR, cutoff_mR)
#' mappingTable = read.csv(paste(directoryPath, "miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv", sep=""), header = FALSE)
#' l <- getMappingList(l, mappingTable)
getMappingList <- function(l, mappingTable) {
  # Get version of miRs
  version=checkMiRNAVersion(l$miRs,verbose=FALSE)
  
  # Convert non mature miRNAs' names to mature names
  miRMature=miRNA_PrecursorToMature(l$miRs, version=version)
  
  # Convert to version 21
  miRNameList_Accession = miRNA_NameToAccession(miRMature[, 2], version=version)
  miRNameList_Converted = miRNA_AccessionToName(miRNameList_Accession[, 2], targetVersion = "v21")
  
  # Process miRNAs' names
  # Update the converted items to the miR list
  # if a converted item is "NA", get the corresponding value from the mature list
  l$miRs <- miRNameList_Converted[,2]
  naList <- which(is.na(miRNameList_Converted[,2]))
  for(i in 1:length(naList)) {
    k <- naList[i]
    l$miRs[k] <- miRMature[k,2]
  }
  
  # Update miR names of the data set
  colnames(l$d)[1:length(l$miRs)] <- l$miRs
  
  # Find mapping information between miRNA and mRNA
  curmiR = NULL;
  curmR = NULL;
  mappingList = NULL;
  for(i in 1:nrow(mappingTable)) {
    curmiR <- as.character(mappingTable$V1[i])
    curmR <- as.character(mappingTable$V2[i])
    
    if(curmiR %in% l$miRs & curmR %in% l$mRs) {
      mappingList=rbind(mappingList,c(curmiR, curmR))
    }
  }
  
  # Sort mappingList based on mRNAs
  mappingList <- mappingList[order(mappingList[,2]),]
  
  # Remove duplicated values
  mappingListFinal <- unique(mappingList)
  
  # Create the returned object
  l = list(d = l$d, miRs = l$miRs, mRs = l$mRs, mappingList = mappingListFinal)
  
  return(l)
}

#' A function for converting miRs to version 21
#'
#' This function allows you to convert miR names to version 21
#' @param miRs miR names.
#' @return converted miR names
#' @keywords 
#' @export
#' @examples
convertmiRs <- function(miRs) {
  # Get version of miRs
  version=checkMiRNAVersion(miRs,verbose=FALSE)
  
  # Convert non mature miRNAs' names to mature names
  miRMature=miRNA_PrecursorToMature(miRs, version=version)
  
  # Convert to version 21
  miRNameList_Accession = miRNA_NameToAccession(miRMature[, 2], version=version)
  miRNameList_Converted = miRNA_AccessionToName(miRNameList_Accession[, 2], targetVersion = "v21")
  
  # Process miRNAs' names
  # Update the converted items to the miR list
  # if a converted item is "NA", get the corresponding value from the mature list
  miRs <- miRNameList_Converted[,2]
  naList <- which(is.na(miRNameList_Converted[,2]))
  for(i in 1:length(naList)) {
    k <- naList[i]
    miRs[k] <- miRMature[k,2]
  }
  
  return(miRs)
}

#' A function for validating the causal effects between miRs and mRNAs
#'
#' This function allows you to validate the causal effects between miRs and mRNAs
#' @param d p values of causal effects between miRs and mRNAs.
#' @param datacsv the true causal effects.
#' @param pVal value to filter causal effects.
#' @return the list of confirmed pairs of miRs and mRNAs which have causal effects
#' @keywords 
#' @export
#' @examples
validate=function(d, datacsv, pVal=0.05){
  
  nmiR <- ncol(d)
  nmR <- nrow(d)
  result = NULL
  for (i in 1:nmiR) {
    for (j in 1:nmR) {
      if(d[j,i] <= pVal) {
        result <- rbind(result, c(colnames(d)[i], row.names(d)[j], d[j,i]))
      }
    }
  }
  
  mappingTable = read.csv(datacsv, header = FALSE)
  l <- nrow(result)
  curmiR = NULL;
  curmR = NULL;
  r = NULL;
  for(i in 1:l) {
    curmiR <- result[i,1]
    curmR <- result[i,2]
    
    if(curmiR %in% as.character(mappingTable$V1) & curmR %in% as.character(mappingTable$V2)) {
      r=rbind(r,result[i,])
    }
  }
  
  return(r)
}

#' A function for converting and processing duplicated miR names
#'
#' This function allows you to convert and process duplicated miR names
#' @param miRs a matrix of expression of miRs with rows being samples and columns being miR names.
#' @return miRs with miR names converted and duplicated miR names removed
#' @keywords 
#' @export
#' @examples
processmiRName=function(miRs){
  
  miRNameList <- colnames(miRs)
  
  # Convert miR names in miRNameList
  # Convert the version of miRNAs to 21
  # Check miRNAs' version
  version=checkMiRNAVersion(miRNameList,verbose=FALSE)
  # Convert non mature miRNAs' names to mature names
  miRMature=miRNA_PrecursorToMature(miRNameList, version=version)
  # Convert to version 21
  miRNameList_Accession = miRNA_NameToAccession(miRMature[, 2], version=version)
  miRNameList_Converted = miRNA_AccessionToName(miRNameList_Accession[, 2], targetVersion = "v21") 
  # Process miRNAs' names
  # Update the converted items to the miR list
  # if a converted item is "NA", get the corresponding value from the mature list
  miRNameList <- miRNameList_Converted[,2]
  naList <- which(is.na(miRNameList_Converted[,2]))
  for(i in 1:length(naList)) {
    k <- naList[i]
    miRNameList[k] <- miRMature[k,2]
  }
  
  # Update miR names in miRs and remove duplicated miR names in miRs
  colnames(miRs) <- miRNameList
  uniquemiRNameList <- unique(miRNameList)
  noOfmiR <- length(uniquemiRNameList)
  temp <- matrix(0, nrow = nrow(miRs), ncol = noOfmiR)
  rownames(temp) <- rownames(miRs)
  colnames(temp) <- uniquemiRNameList
  for(column in 1:noOfmiR){
    miRList <- which(miRNameList == colnames(temp)[column])
    for(r in 1:nrow(temp)) {
      temp[r, column] <- mean(miRs[r, miRList])  
    }
  }
  
  return(temp)
}

#' A function for processing duplicated mRNA names
#'
#' This function allows you to process duplicated mRNA names
#' @param mRNAs a matrix of expression of mRNAs with rows being samples and columns being mRNA names.
#' @return mRNAs with duplicated mRNA names removed
#' @keywords 
#' @export
#' @examples
processmRNAName=function(mRNAs){
  
  mRNANameList <- colnames(mRNAs)
  uniquemRNANameList <- unique(mRNANameList)
  noOfmRNA <- length(uniquemRNANameList)
  temp <- matrix(0, nrow = nrow(mRNAs), ncol = noOfmRNA)
  rownames(temp) <- rownames(mRNAs)
  colnames(temp) <- uniquemRNANameList
  for(column in 1:noOfmRNA){
    mRNAList <- which(mRNANameList == colnames(temp)[column])
    for(r in 1:nrow(temp)) {
      temp[r, column] <- mean(mRNAs[r, mRNAList])  
    }
  }
  
  return(temp)
}

#' A function for preparing data to classify cancer subtypes
#'
#' This function allows you to prepare data for classifying cancer subtypes
#' @param d a matrix of expression of 50 mRNAs in Pam50 with rows being samples and columns being mRNA names.
#' @param directoryPath the directory path to save result file.
#' @param dataDir the directory for saving the result data.
#' @param inputFileName the file name of the input data which is used to classify cancer subtypes.
#' @return none, just writing the resullt to file
#' @keywords 
#' @export
#' @examples
prepareData=function(d, directoryPath, dataDir = "bioclassifier_data", inputFileName = "inputFile.txt"){
  
  d <- t(d)
  noOfRow <- nrow(d)
  noOfCol <- ncol(d)
  temp <- matrix(0, nrow = (noOfRow + 1), ncol = noOfCol)
  row.names(temp) <- c("T", row.names(d))
  colnames(temp) <- colnames(d)
  temp[2:(noOfRow + 1), 1:noOfCol] <- d[,]
  write.table(temp, paste(directoryPath, "/", dataDir, "/", inputFileName, sep = ""), col.names=NA, sep = "\t")
}

#' Validate the targets of miRNAs using transfection data
#'
#' Given the predicted target of miRNAs, the function returns a list of targets that are  confirmed
#' based on the curated transfection data. Users need to download the file logFC.imputed.rda from nugget.unisa.edu.au/Thuc/miRLAB/ and place it in the working directory (this file is obtained from the TargetScoreData package)
#' @param topkList a matrix with the following columns. The first column is the miRNA name, the second contains the target mRNAs, and the third contains the correlation values/ causal effects/ scores
#' @param LFC the log fold change threshold. The targets that have the absolute value of log fold change greater than the LFC
#' will be regarded as the confirmed targets.
#' @return a matrix in the same format of the input matrix but only contains the confirmed interactions.
#' @keywords 
#' @export
#' @examples
ValidationTPlus=function(topkList, LFC){
  
  uni <- unique(topkList[,1])
  l <- length(uni)
  result = NULL
  for(i in 1:l) {
    curmiR <- uni[i]
    curPairs <- topkList[which(topkList[,1] == curmiR),]
    curPairs <- matrix(curPairs, ncol=4);
    curConfirmed = ValidationT(curPairs, LFC)
    result <- rbind(result, curConfirmed[[1]])
  }
  return(result)
}

ValidationTPlusForAll=function(directoryPath, infile, outfile1, outfile2) {
  # Validate with transfection data
  # Validate the results of the topk targets of miRNAs predicted by the method
  
  setwd(directoryPath)
  
  # Get the result
  temp <- read.csv(file = paste(directoryPath, infile, sep = ""), row.names = 1)
  temp <- data.matrix(temp)
  nmiR <- ncol(temp)
  
  for(topk in c(500, 1000, 1500, 2000)) {
    top = Extopk(temp, topk)
    topConfirmed = ValidationTPlus(top, 0.3)
    write.csv(topConfirmed, file = paste(directoryPath, outfile1, topk, ".csv", sep = ""))
  }
  
  ##======================
  # Validate for each miRNA
  for(k in c(50, 100, 150, 200)) {
    topConfirmed = NULL
    for(i in 1:nmiR) {
      miRiTopk = bRank(temp, i, k, FALSE)
      topConfirmedi = ValidationT(miRiTopk, 0.3)
      if(!is.null(topConfirmedi[[1]])) {
        topConfirmed <- rbind(topConfirmed, topConfirmedi[[1]])
      }
    }
    write.table(topConfirmed, file = paste(directoryPath, outfile2, k, ".csv", sep = ""), sep=",", row.names = FALSE, col.names = FALSE)
  }
  ##======================
}

checkBRCAInteraction=function(topkList, BRCAmiRs, BRCAmRNAs){
  
  l <- nrow(topkList)
  result = NULL
  for(i in 1:l) {
    curmiR <- topkList[i,1]
    curmR <- topkList[i,2]
    if(curmiR %in% BRCAmiRs$V1 & curmR %in% BRCAmRNAs$V1) {
      result <- rbind(result, c(curmiR, curmR))
    }
  }
  return(result)
}

checkBRCAInteractionForAll=function(directoryPath, infile, outfile1, outfile2) {
  # Validate with BRCA interactions
  # Validate the results of the topk targets of miRNAs predicted by the method
  
  # Get the result
  temp <- read.csv(file = paste(directoryPath, infile, sep = ""), row.names = 1)
  temp <- data.matrix(temp)
  nmiR <- ncol(temp)
  BRCAmiRs <- read.csv(file = paste(directoryPath, "BRCA_miRNA.csv", sep = ""), header = FALSE)
  BRCAmRNAs <- read.csv(file = paste(directoryPath, "BRCA_GENE.csv", sep = ""), header = FALSE)
  
  for(topk in c(500, 1000, 1500, 2000)) {
    top = Extopk(temp, topk)
    topConfirmed = checkBRCAInteraction(top, BRCAmiRs, BRCAmRNAs)
    write.csv(topConfirmed, file = paste(directoryPath, outfile1, topk, ".csv", sep = ""))
  }
  
  ##======================
  # Validate for each miRNA
  for(k in c(50, 100, 150, 200)) {
    topConfirmed = NULL
    for(i in 1:nmiR) {
      miRiTopk = bRank(temp, i, k, FALSE)
      miRiTopk <- as.matrix(miRiTopk)
      topConfirmedi = checkBRCAInteraction(miRiTopk, BRCAmiRs, BRCAmRNAs)
      if(!is.null(topConfirmedi)) {
        topConfirmed <- rbind(topConfirmed, topConfirmedi)
      }
    }
    write.table(topConfirmed, file = paste(directoryPath, outfile2, k, ".csv", sep = ""), sep=",", row.names = FALSE, col.names = FALSE)
  }
  ##======================
}

getList=function(directoryPath, fileName, hasRowName = TRUE, hasColName = TRUE) {
  
  if(hasRowName & hasColName) {
    f <- read.csv(paste(directoryPath, fileName, sep = ""), row.names = 1)
  } else if(hasRowName & !hasColName) {
    f <- read.csv(paste(directoryPath, fileName, sep = ""), row.names = 1, header = FALSE)
  } else if(!hasRowName & hasColName) {
    f <- read.csv(paste(directoryPath, fileName, sep = ""))
  } else {
    f <- read.csv(paste(directoryPath, fileName, sep = ""), header = FALSE)
  }
  
  len <- nrow(f)
  l <- ""
  for(i in 1:len) {
    l <- paste(l, f[i,1], f[i,2], sep="")
    if(i != len) {
      l <- paste(l, ",", sep = "")
    }
  }
  
  return(l)
}