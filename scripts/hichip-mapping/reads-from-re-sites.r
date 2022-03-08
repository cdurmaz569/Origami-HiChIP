library(Rsamtools)
library(matrixStats)

args <- commandArgs(T)

if( length(args) < 3) {
  stop("Syntax: <RE file> <BAM file> <output directory>")
}

reFile <- args[1]
bamFile <- args[2]
outputDir <- args[3]

reDf <- read.table(reFile,sep='\t')
rePos <- split(reDf[,2],reDf[,1])
rm(reDf)


header <- scanBamHeader(bamFile)
targets <- header[[1]][[1]]

params <- ScanBamParam(what=c("flag","rname","pos","seq"))
rd <- scanBam(bamFile,param=params)[[1]]
w <- width(rd$seq)
wmax <- max(w)
wmb <- rowAlls(matrix(w==wmax,ncol=2,byrow=T))

mdist <- matrix(unlist(mapply(function(p,chr) {
  if(is.na(chr) || is.na(p)) return(NA)
  sort.int(abs(p-rePos[[chr]]),method="radix")[1] ## get distance to nearest RE site
},rd$pos,as.character(rd$rname),SIMPLIFY=F)),ncol=2,byrow=T)



if( !dir.exists(outputDir)) {
  dir.create(outputDir)
}

pdf(paste(outputDir,"distance-of-reads-from-RE-site.pdf",sep="/"),width=8,height=8)
layout(matrix(1:3,ncol=1))
d1 <- density(mdist,na.rm=T,from=0)
d2 <- density(mdist[!wmb,],na.rm=T,from=0)
d3 <- density(mdist[wmb,],na.rm=T,from=0)

m1 <- mean(mdist,na.rm=T)
m2 <- mean(mdist[!wmb,],na.rm=T)
m3 <- mean(mdist[wmb,],na.rm=T)

md1 <- median(mdist,na.rm=T)
md2 <- median(mdist[!wmb,],na.rm=T)
md3 <- median(mdist[wmb,],na.rm=T)


ymax <- max(c(d1$y,d2$y,d3$y))

plot(d1,ylim=c(0,ymax),xlim=c(0,1000),main="All PET ends",xlab="Distance to nearest RE site (bp)",frame.plot=F)
abline(v=0)
abline(v=m1,col='red')
abline(v=md1,col='blue')
legend("topright",legend=paste0(c("Mean (","Median ("),round(c(m1,md1),digits=2),c(")",")")),col=c("red","blue"),lty=1)

plot(d2,ylim=c(0,ymax),xlim=c(0,1000),main="PETs ends with junction detected",xlab="Distance to nearest RE site (bp)",frame.plot=F)
abline(v=0)
abline(v=m2,col='red')
abline(v=md2,col='blue')
legend("topright",legend=paste0(c("Mean (","Median ("),round(c(m2,md2),digits=2),c(")",")")),col=c("red","blue"),lty=1)

plot(d3,ylim=c(0,ymax),xlim=c(0,1000),main="PETs ends with no junction detected",xlab="Distance to nearest RE site (bp)",frame.plot=F)
abline(v=0)
abline(v=m3,col='red')
abline(v=md3,col='blue')
legend("topright",legend=paste0(c("Mean (","Median ("),round(c(m3,md3),digits=2),c(")",")")),col=c("red","blue"),lty=1)
layout(1)
dev.off()

