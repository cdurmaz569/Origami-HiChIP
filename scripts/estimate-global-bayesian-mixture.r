library(GenomicRanges,quietly = !interactive())
library(utils,quietly = !interactive())
library(matrixStats,quietly = !interactive())
library(parallel)
library(pbapply)

calc.zscore <- function(v) (v-mean(v))/sd(v)
bounded.prob <- function(v,low,high) pmin(pmax(v,low),high)

chain.process <- function(counts,lambda0,lambda1,lambdad0,lambdad1,pp,S,alphaparam,betaparam,
                          with.distance.weight,intdist,interchromosomal,minintdist,usedf,
                          suppress,mini.model,show.progress) {
  #lambda0 <- ret$lambda0[i]
  #lambda1 <- ret$lambda1[i]
  
  ret <- list()
  
  g1 <- pp*dpois(counts,lambda1 + lambdad1)
  g2 <- (1-pp) * dpois(counts,lambda0 + lambdad0)
  
  vp <- g1/rowSums(cbind(g1,g2)) 
  
  # Sometimes the counts value in dpois is so extreme the loss of precision causes the value to be 0, need to correct for this
  # if in this case the counts value is closer to lambda1 than lambda0, give it a probability of .999, otherwise .001
  if(any(is.na(vp))) {
    b<- is.na(vp)
    vp[b] <- ifelse(counts[b]>=lambda1,.999,.001)
  }
  
  ret[["z"]] <- vz <- rbinom(S,1,vp)
  ret[["pp"]] <- pp <- rbeta(S,alphaparam+vz,betaparam+(1-vz))
  
  lambda0prob <- bounded.prob(lambda0/(lambda0+lambdad0),.001,.999)
  lambda1prob <- bounded.prob(lambda1/(lambda1+lambdad1),.001,.999)
  
  groupcounts0 <- rbinom(length(counts),counts,lambda0prob)
  distcounts0 <- counts-groupcounts0
  
  groupcounts1 <- rbinom(length(counts),counts,lambda1prob)
  distcounts1 <- counts-groupcounts1
  
  isgroup0 <- vz == 0
  isgroup1 <- !isgroup0

  r0 <- sum(groupcounts0[isgroup0&!suppress])
  n0 <- sum(isgroup0&!suppress)
  l0 <- rgamma(1,r0,n0)
  
  l1 <- l0
  
  r1 <- sum(groupcounts1[isgroup1&!suppress])
  n1 <- sum(isgroup1&!suppress)
  l1x <- rgamma(1,r1,n1)
  l1 <- max(l1,l1x)
  
  ret[["lambda0"]] <- l0
  ret[["lambda1"]] <- l1
  
  if(with.distance.weight) {
    x <- log10(intdist[vz==1 & !interchromosomal]+1)
    
    s1 <- if( usedf > 0 ) smooth.spline(x,pmax(distcounts1[vz==1& !is.na(intdist)],0),df=usedf) else smooth.spline(x,pmax(distcounts1[vz==1& !is.na(intdist)],0))
    

    x <- log10(intdist[vz==0 & !interchromosomal]+1)
    
    s0 <- if( usedf > 0 ) smooth.spline(x,pmax(distcounts0[vz==0& !is.na(intdist)],0),df=usedf) else smooth.spline(x,pmax(distcounts0[vz==0& !is.na(intdist)],0))
    
    x <- log10(intdist+1)
    if(any(interchromosomal)) x[interchromosomal] <- log10(minintdist+1) ## set eact interchromsomal interaction to shortest distance (which should have the highest mean read count)
    
    lambdad1 <- pmax(predict(s1,x)$y,0) ### floor the value at 0
    lambdad0 <- pmax(predict(s0,x)$y,0) 
    
    ret[["lambdad1"]] <- lambdad1
    ret[["lambdad0"]] <- lambdad0
    
  }
  #if(show.progress) setTxtProgressBar(pb,i/N)
  
  ret
}

estimate.global.bayesian.mixture <- function(ints,depth,inttable,N=1000,burnin=100,pruning=NULL,
                                                        with.distance.weight=F,multiply=T,show.progress=F,
                                                        lambda0.init=1,lambda1.init=5,suppress.counts.higher.than=30,
                                                        mini.model=T,usedf=0,mode="long",
                                             parallel=F,parallel.cores=1) {
  if(!mini.model) stop("Only mini model storage is available right now")
  
  S <- nrow(ints)
  
  d <- GRanges(seqnames=as.character(depth$V1),ranges=IRanges(depth$V2,depth$V3),strand='*')
  g1 <- GRanges(seqnames=as.character(ints$V1),ranges=IRanges(ints$V2,ints$V3),strand='*')
  g2 <- GRanges(seqnames=as.character(ints$V4),ranges=IRanges(ints$V5,ints$V6),strand='*')
  
  m1 <- match(g1,d)
  m2 <- match(g2,d)
  counts <- ints[,7]
  d <- depth[,4]
  intdist <- distance(g1,g2)
  interchromosomal <- is.na(intdist)
  minintdist <- min(intdist[!interchromosomal])
  
  if(is.null(pruning) || pruning < 1) {
    pruning <- N+1 # will never prune
  } 
  
  if(!multiply) {
    sdepth <- rowSums(cbind(d[m1],d[m2]))
    msdepth <- median(sdepth)
  } else {
    sdepth <- rowProds(cbind(d[m1],d[m2]))
    msdepth <- median(sdepth)
  }
  
  f1 <- paste(ints$V1,ints$V2,ints$V3,sep='_')
  f2 <- paste(ints$V4,ints$V5,ints$V6,sep='_')
  
  m1 <- match(f1,rownames(inttable))
  m2 <- match(f2,rownames(inttable))
  
  depthscore <- floor(sdepth/msdepth)
  
  if(mode!="hichip") {
    countm <- t(mapply(function(x1,x2,n) {
      idx <- as.integer(colnames(inttable))<n
    
      z <- inttable[c(x1,x2),]

      a <- if(!any(idx)) 0 else sum(z[,idx])
      b <- sum(z)-a
      c(a,b)
    },as.list(m1),as.list(m2),as.list(counts),SIMPLIFY=T))
  
    alphaparam <- 1+countm[,1]
    betaparam <- 1+depthscore+countm[,2]
  } else {
    countm <- t(mapply(function(x1,x2,n) {
      idx <- as.integer(colnames(inttable))<n
      
      a <- sum(inttable[c(x1,x2),] > 0 & idx)
      b <- sum(inttable[c(x1,x2),] > 0 & !idx)
      
      #z <- inttable[c(x1,x2),]
      
      #a <- if(!any(idx)) 0 else sum(z[,idx])
      #b <- sum(z)-a
      c(a,b)
    },as.list(m1),as.list(m2),as.list(counts),SIMPLIFY=T))
    
    alphaparam <- 1+countm[,1]
    betaparam <- 1+countm[,2]

    rm(countm,depthscore,m1,m2,f1,f2,d,g1,g2)
    gc()
  }
  

  
  pp <- rep(.5,S)
  
  if (!mini.model) {
    ret <- list(
      z = vector("list",N),
      p1 = c(list(pp),vector("list",N)),
      mp = vector("list",N),
      lambda0 = c(lambda0.init,rep(NA_real_,N)),
      lambda1 = c(lambda1.init,rep(NA_real_,N)),
      lambdad1 = vector("list",N),
      lambdad0 = vector("list",N)
    )
  } else {
    burnin.l <- list(
      lambda0 = c(lambda0.init,rep(NA_real_,burnin)),
      lambda1 = c(lambda1.init,rep(NA_real_,burnin))
    )
    
    last.burnin <- list(lambda0=lambda0.init,lambda1=lambda1.init)
    
    #zm <- rep(0,length(counts))
  }
  
  totcounts <- sum(counts)
  
  lambdad1 <- lambdad0 <- rep(0,length(counts))
  
  cat("Burn-in...\n")
  
  if(show.progress) pb <- txtProgressBar()
  
  suppress <- counts > suppress.counts.higher.than
  
  for( i in 1:burnin ) {
    l <- chain.process(counts,burnin.l$lambda0[[i]],burnin.l$lambda1[[i]],
                       lambdad0,lambdad1,pp,S,alphaparam,betaparam,
                       with.distance.weight,intdist,interchromosomal,minintdist,usedf,
                       suppress,mini.model,show.progress)
    
    burnin.l$lambda0[[i+1]] <- l$lambda0
    burnin.l$lambda1[[i+1]] <- l$lambda1
    
    lambdad0 <- l$lambdad0
    lambdad1 <- l$lambdad1
    pp <- l$pp
    
    last.burnin <- l
    
    if( show.progress ) setTxtProgressBar(pb,i/burnin)
  }
  
  if(show.progress) close(pb)
  
  cat("MCMC...\n")
  
  
  if(parallel) {
    Niter <- floor(N/parallel.cores)
    
    il <- pblapply(rep(Niter,parallel.cores),function(v,lb,lambdad0,lambdad1) {
      last.iter <- lb
      
      ret <- list(
        lambda0 = c(rep(NA_real_,v)),
        lambda1 = c(rep(NA_real_,v))
      )
      
      zm <- rep(0,length(counts))
      
      for( i in 1:v ) {
        l <- chain.process(counts,last.iter$lambda0,last.iter$lambda1,
                           lambdad0,lambdad1,pp,S,alphaparam,betaparam,
                           with.distance.weight,intdist,interchromosomal,minintdist,usedf,
                           suppress,mini.model,show.progress)
        
        ret$lambda0[[i]] <- l$lambda0
        ret$lambda1[[i]] <- l$lambda1
        
        lambdad0 <- l$lambdad0
        lambdad1 <- l$lambdad1
        pp <- l$pp
        zm <- zm+l$z
        
        last.iter <- l
      }
      
      ret[["zm"]] <- zm
      ret[["lambdad0"]] <- lambdad0
      ret[["lambdad1"]] <- lambdad1
      
      ret
    },lb=last.burnin,lambdad0=lambdad0,lambdad1=lambdad1,cl=parallel.cores)
    #mc.preschedule = T,mc.cores=parallel.cores)
    
    if(any(e <- sapply(il,function(r) inherits(r, 'try-error')))){
      print(il[e]) 
  
      stop("Error encountered during paralleization")
    }
    
    zmat <- do.call(cbind,lapply(il,function(l) l$zm))
    
    lambdad0l <- lapply(il,function(l) l$lambdad0)
    lambdad1l <- lapply(il,function(l) l$lambdad1)
    
    
    ret <- list(lambda0 = do.call(c,lapply(il,function(l) l$lambda0)),
                lambda1 = do.call(c,lapply(il,function(l) l$lambda1)),
                zm=rowSums(zmat),
                lambdad0=lambdad0l[[length(lambdad0l)]],
                lambdad1=lambdad1l[[length(lambdad1l)]])
  } else {
    ret <- list(
      lambda0 = c(lambda0.init,rep(NA_real_,N)),
      lambda1 = c(lambda1.init,rep(NA_real_,N))
    )
    last.iter <- last.burnin
    zm <- rep(0,length(counts))
    
    if(show.progress) pb <- txtProgressBar()
    
    for( i in 1:N ) {
      l <- chain.process(counts,ret$lambda0[[i]],ret$lambda1[[i]],
                         lambdad0,lambdad1,pp,S,alphaparam,betaparam,
                         with.distance.weight,intdist,interchromosomal,minintdist,usedf,
                         suppress,mini.model,show.progress)
      
      ret$lambda0[[i+1]] <- l$lambda0
      ret$lambda1[[i+1]] <- l$lambda1
      
      lambdad0 <- l$lambdad0
      lambdad1 <- l$lambdad1
      pp <- l$pp
      zm <- zm+l$z
      
      last.iter <- l
      
      if(show.progress) setTxtProgressBar(pb,i/N)
    }
    
    ret[["zm"]] <- zm
    ret[["lambdad0"]] <- lambdad0
    ret[["lambdad1"]] <- lambdad1
    
    if(show.progress) close(pb)
  }
  
  
  dl <- list(burnin=burnin.l,lambda0=ret$lambda0,lambda1=ret$lambda1,zm=ret$zm,
             lambdad0=ret$lambdad0,lambdad1=ret$lambdad1)
  #ret <- list(s=l,l0=lambda0,l1=lambda1)


  dl <- c(dl,list(sdepth=sdepth,msdepth=msdepth,intdist=intdist))
  #if(mini.model) ret <- c(ret,list(zm=zm,lambdad1=lambdad1,lamdbad0=lamdbad0))
  
  dl
}

extract.global.bayesian.mixture.prob <- function(model) {
  if( "zm" %in% names(model)) {
    return(model$zm / length(model$lambda0))
  } else {
    m <- do.call(rbind,model$z)
    N <- nrow(m)
    return(colSums(m)/N)
  }
}
