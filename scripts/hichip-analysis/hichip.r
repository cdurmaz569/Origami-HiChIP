suppressMessages(library(rhdf5))
suppressMessages(library(bigmemory))


#### constants for the table structure

TABLE.BIN1 <- 1
TABLE.BIN2 <- 2
TABLE.IDX.DIFF <- 3
TABLE.PETCOUNT <- 4
TABLE.DISTANCE <- 5
TABLE.RE.CUT.NUM.BIN1 <- 6
TABLE.RE.CUT.NUM.BIN2 <- 7
TABLE.RE.CUT.TOTAL <- 8
TABLE.CHIP.BIN1 <- 9
TABLE.CHIP.BIN2 <- 10
TABLE.CHIP.TOTAL <- 11
TABLE.SELF <- 12
TABLE.WIDEREGION <- 13
TABLE.LOCALREGION <- 14
TABLE.SURROUNDING.NEIGHBORHOOD <- 15
TABLE.LOWERLEFT <- 16
TABLE.VERTICAL <- 17
TABLE.HORZIONTAL <- 18


DATATABLE.MATRIX <- c(TABLE.BIN1,TABLE.BIN2,TABLE.IDX.DIFF,TABLE.PETCOUNT,TABLE.DISTANCE)

GLM.MATRIX <- c(TABLE.CHIP.BIN1,TABLE.CHIP.BIN2,TABLE.RE.CUT.TOTAL,
                TABLE.SURROUNDING.NEIGHBORHOOD,TABLE.LOWERLEFT,TABLE.VERTICAL,TABLE.HORZIONTAL)



####

if( interactive()) 
  warning("This script is not meant to be run interactively")

execArgs <- commandArgs(T)

if(length(execArgs) < 18) {
  stop("Not enough arguments")
}

outputDir <- execArgs[1]
filePrefix <- execArgs[2]
peakFile <- execArgs[3]
intFile <- execArgs[4]
iterations <- as.integer(execArgs[5])
burnin <- as.integer(execArgs[6])
prune <- as.integer(execArgs[7])
saveMiniModel <- execArgs[8] == "yes"
verbose <- execArgs[9] == 'yes'
mergeRowsLessThan <- as.integer(execArgs[10])
mergeRowsAfter <- as.integer(execArgs[11])
minBinsWithPeak <- as.integer(execArgs[12])
glmEpsilon <- as.numeric(execArgs[13])
fileBacked <- as.logical(execArgs[14])
excludeNotWithinDist <- as.integer(execArgs[15])
onlyIntrachromosomal <- execArgs[16] == "yes"
parallel <- execArgs[17] == "yes"
parallelCores <- as.integer(execArgs[18])

if(verbose) {
  options(warn = 1) ### print all warnings
}

outputPrefix <- paste0(outputDir,"/",filePrefix)


modelFile <- paste(outputPrefix,"model.Rdata",sep='-')
outputFile <- paste(outputPrefix,"results.txt",sep='-')

### get execution path
fullArgs <- commandArgs()
scriptLocation <- sub("--file=","",fullArgs[grep("--file=",fullArgs)])
dbase <- dirname(scriptLocation)


source(paste(dbase,"hichip-model.r",sep='/'))


options <- list()

options[["outputDir"]] <- outputDir
options[["filePrefix"]] <- filePrefix
options[["iterations"]] <- iterations
options[["burnin"]] <- burnin
options[["pruning"]] <- prune
options[["showProgress"]] <- T
options[["mergeRowsLessThan"]] <- mergeRowsLessThan
options[["mergeRowsAfter"]] <- mergeRowsAfter
options[["miniModel"]] <- saveMiniModel
options[["parallel"]] <- parallel
options[["parallelCores"]] <- parallelCores

options[["minBinsWithPeak"]] <- minBinsWithPeak 

options[["glm.epsilon"]] <- glmEpsilon 

options[["filebacked"]] <- fileBacked


if(options$glm.epsilon<=0) {
  stop("The GLM epslion must be non-negative")
}

cat("Loading data files...\n")


dataOM <- read.big.matrix(intFile,
                         col.names=c("bin1","bin2","idxdiff","petcount","distance","recut1","recut2","retotal","chip1","chip2","chiptotal",
                                     "selfpets","wideregion","localregion","surrouding","lowerleft","vertical","horizontal"),
                         sep='\t',
                         type='integer',shared=F)


cat("Identifying ligation events...\n")

data <- list()

glm.bin.file <- paste0(filePrefix,"-glm.bin")
glm.desc.file <- paste0(filePrefix,"-glm.desc")


bin.file <- paste0(filePrefix,"-dataM.bin")
desc.file <- paste0(filePrefix,"-dataM.desc")

includeRows <- mwhich(dataOM,TABLE.IDX.DIFF,0,'gt')

pos <- rep(T,nrow(dataOM))

if( excludeNotWithinDist > 0 ) {
  z <- sapply(split(dataOM[,TABLE.IDX.DIFF],dataOM[,TABLE.BIN1],drop=F),function(v) any(v[v>0] <= excludeNotWithinDist))
  valid <- as.integer(names(z)[z])
  pos <- dataOM[,TABLE.BIN1] %in% valid
  rm(z,valid)
}

if( onlyIntrachromosomal ) {
  pos <- pos & !is.na(dataOM[,TABLE.DISTANCE])
}

pos <- which(pos[includeRows])

### remove self-bin hits

dataM <- data[["dataM"]] <- if(options$filebacked) as.big.matrix(cbind(dataOM[includeRows,DATATABLE.MATRIX][pos,],as.integer(1)),
                                                                 type='integer',
                                                                 backingpath=outputDir,
                                                                 backingfile=bin.file,
                                                                 descriptor=desc.file) else 
                                                                   as.big.matrix(cbind(dataOM[includeRows,DATATABLE.MATRIX][pos,],as.integer(1)),
                                               type='integer',
                                               shared=F)


temp <- cbind(1,dataOM[includeRows,GLM.MATRIX][pos,])
  
rm(dataOM)

class(temp) <- "double"
data[["dataGLM"]] <- if( options$filebacked ) as.big.matrix(temp,type='double',backingpath=outputDir,backingfile=glm.bin.file, descriptor=glm.desc.file) else 
                                              as.big.matrix(temp,type='double',shared=F)

rm(temp,pos,includeRows)


cat("Running mixture model...\n")



gbayes <- estimate.global.bayesian.mixture(data,
                                           outputPrefix,
                                           options)
gbayesp <- gbayes$prob
expected <- gbayes$expected


cat("Saving additional model output...\n")
save(gbayes,file=modelFile)


rm(gbayes)

cat("Loading bin locations...\n")

peaktable <- read.table(peakFile,sep='\t',colClasses = c("factor","integer","integer","NULL","NULL","NULL"))

cat("Estimating normalized values...\n")

obsp <- expp <- rep(0,nrow(peaktable))

nm <- cbind(dataM[,1],dataM[,2],dataM[,4],expected)

tmp <- sapply(split(rep(nm[,3],2),c(nm[,1],nm[,2])),sum)
obsp[as.integer(names(tmp))] <- as.vector(tmp)

tmp <- sapply(split(rep(nm[,4],2),c(nm[,1],nm[,2])),sum)
expp[as.integer(names(tmp))] <- as.vector(tmp)

normr <- (obsp+1e-6)/(expp+1e-6)
petr <- dataM[,4]/expected/normr[dataM[,1]]/normr[dataM[,2]]
#base <- normr[,dataM[,1]]/normr[,dataM[,2]]

cat("Writing out results table...\n")

leftanchor <- peaktable[dataM[,1],]
rightanchor <- peaktable[dataM[,2],]



m <- cbind(leftanchor,rightanchor,dataM[,c(1,2,4)],expected,petr,dataM[,5],gbayesp)

colnames(m) <- c("chromosome1","start1","end1","chromosome2","start2","end2","anchor idx1","anchor idx2","PET count","Expected PET Count","Normalized PET rate","Loop length","MM probability")

write.table(m,file=outputFile,row.names=F,quote=F,sep='\t')


