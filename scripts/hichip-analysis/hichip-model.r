suppressMessages(library(GenomicRanges))
suppressMessages(library(utils))
suppressMessages(library(matrixStats))
suppressMessages(library(parallel))
suppressMessages(library(pbapply))
suppressMessages(library(rhdf5))
suppressMessages(library(bigmemory))
suppressMessages(library(bigalgebra))
suppressMessages(library(biganalytics))
suppressMessages(library(VGAM))
suppressMessages(library(Rcpp))

##### C++


cppFunction('

List updateGroupStructure(S4 dMat, 
                                   NumericVector l0, 
                                   NumericVector l1, 
                                   NumericVector mixturep, // Will be modified
                                   IntegerVector zv, // Will be modified
                                   IntegerVector gid, // will be modified
                                   IntegerVector bip,
                                   NumericVector alpha,
                                   NumericVector beta) {
  std::vector<int> wisgroup0, wisgroup1;
  XPtr<BigMatrix> pMat(Rcpp::as<SEXP>(dMat.slot("address")));
  MatrixAccessor<int> iAcc(*pMat);
  int PETCOLUMN=3;

  for( int i = 0; i < l0.size(); i++ ) {
    int p = iAcc[PETCOLUMN][i];
    double lambda0 = l0[i];
    double lambda1 = l1[i];
    double mp = mixturep[i];

    //Use a zero-truncated Poisson
    double g0 = (1-mp) * R::dpois(p,lambda0,false)/(1-R::dpois(0,lambda0,false));
    double g1 = mp * R::dpois(p,lambda1,false)/(1-R::dpois(0,lambda1,false));

    double l = g1/(g0+g1);
    if(std::isnan(l)) {
      if(p>lambda1) l = .999;
      else l = .001;
    }

    // Bound the probabilities
    if(l > .999 ) l = .999;
    if(l < .001 ) l = .001;

    int ga = 0;
    if(bip[i] == 1) ga = R::rbinom(1,l);

    gid[i] = ga;
    zv[i] += ga;

    mixturep[i] = R::rbeta(alpha[i]+ga,beta[i]+(1-ga));

    //This will return index lists of which contact is in which group (gor this iteration)
    //Add one to each because R using 1-based coordinates
    if(ga==0) wisgroup0.push_back(i+1);
    else wisgroup1.push_back(i+1);
  }
  
  return Rcpp::List::create(Rcpp::Named("wisgroup0") = Rcpp::IntegerVector::import(wisgroup0.begin(),wisgroup0.end()),
                            Rcpp::Named("wisgroup1") = Rcpp::IntegerVector::import(wisgroup1.begin(),wisgroup1.end()));
}',plugins=c("cpp11"),depends=c("BH","bigmemory"),
            includes=c('#include "bigmemory/BigMatrix.h"','#include "bigmemory/MatrixAccessor.hpp"')  )


cppFunction('
SEXP updatePeakDist(IntegerVector bothinpeak,
                    IntegerVector binpc,
                    S4 dMat,
                    IntegerVector p1lambdav,
                    IntegerVector p2lambdav,
                    IntegerVector minbinsv) 
{
  double p1lambda = p1lambdav[0];
  double p2lambda = p2lambdav[0];
  int minbins = minbinsv[0];
  int CHIP1=1, CHIP2=2;
  XPtr<BigMatrix> pMat(Rcpp::as<SEXP>(dMat.slot("address")));
  MatrixAccessor<double> iAcc(*pMat);


  for( int i = 0; i < bothinpeak.size(); i++ ) {
    int chip1 = iAcc[CHIP1][i], chip2 = iAcc[CHIP2][i];

    double p1p = R::ppois(chip1,p1lambda,1,0);
    double p2p = R::ppois(chip2,p2lambda,1,0);

    if(std::isnan(p1p)) {
      if(chip1>p1lambda) p1p = .999;
      else p1p = .001;
    }

    if(std::isnan(p2p)) {
      if(chip2>p2lambda) p2p = .999;
      else p2p = .001;
    }

    int b1 = R::rbinom(1,p1p);
    int b2 = R::rbinom(1,p2p);

    int b = b1+b2>=minbins ? 1 : 0;

    bothinpeak[i] = b;
    binpc[i] += b;
  }

  return R_NilValue;

}
',plugins=c("cpp11"),depends=c("BH","bigmemory"),
            includes=c('#include "bigmemory/BigMatrix.h"','#include "bigmemory/MatrixAccessor.hpp"') )

cppFunction('SEXP swapLambdas(NumericVector lambda0,NumericVector lambda1) {
  for( int i = 0; i < lambda0.size(); i++ ) {
    if(lambda0[i] < lambda1[i]) continue;
    double tmp = lambda0[i];
    lambda0[i] = lambda1[i];
    lambda1[i] = tmp;
  }

  return R_NilValue;
}',plugins=c("cpp11"))

cppFunction('SEXP updateExpected(NumericVector lambda0, NumericVector lambda1,NumericVector mixturep, NumericVector expected) {
  for( int i = 0; i < lambda0.size(); i++ ) {
    double mp = mixturep[i];

    if(lambda0[i] > lambda1[i]) {
      double tmp = lambda0[i];
      lambda0[i] = lambda1[i];
      lambda1[i] = tmp;
    }

    expected[i] += mp * lambda1[i] + (1-mp) * lambda0[i];
  }

  return R_NilValue;
}',plugins=c("cpp11"))

cppFunction('
SEXP updateLambdas(NumericVector lambda0, //modified
                   NumericVector lambda1, //modified
                             S4 glmMat, //not modified
                             S4 dataMat, //not modified
                             NumericMatrix coefs0, //not modified
                             NumericMatrix coefs1) //not modified
{
  XPtr<BigMatrix> pGlm(Rcpp::as<SEXP>(glmMat.slot("address")));
  XPtr<BigMatrix> pDmat(Rcpp::as<SEXP>(dataMat.slot("address")));
  MatrixAccessor<double> gAcc(*pGlm);
  MatrixAccessor<int> dAcc(*pDmat);

  int nLevels = coefs0.nrow();
  int nCoefs = coefs0.ncol();
  int nContacts = lambda0.size();

  //lambda0.fill(0);
  //lambda1.fill(0);

  double *l0 = lambda0.begin(), *l1 = lambda1.begin();
  double *c0 = coefs0.begin(), *c1 = coefs1.begin();
  int *intIdx = dAcc[5];

  double *v0 = gAcc[0];
  double *v1 = gAcc[1];
  double *v2 = gAcc[2];
  double *v3 = gAcc[3];
  double *v4 = gAcc[4];
  double *v5 = gAcc[5];
  double *v6 = gAcc[6];
  double *v7 = gAcc[7];


  for( int i = 0; i < nContacts; i++ ) {
    int iIdx = intIdx[i]-1;
    l0[i] = std::exp(v0[i]*c0[iIdx] + 
                     v1[i]*c0[nLevels+iIdx] +
                     v2[i]*c0[2*nLevels+iIdx] +
                     v3[i]*c0[3*nLevels+iIdx] +
                     v4[i]*c0[4*nLevels+iIdx] +
                     v5[i]*c0[5*nLevels+iIdx] +
                     v6[i]*c0[6*nLevels+iIdx] +
                     v7[i]*c0[7*nLevels+iIdx]);

    l1[i] = std::exp(v0[i]*c1[iIdx] + 
                     v1[i]*c1[nLevels+iIdx] +
                     v2[i]*c1[2*nLevels+iIdx] +
                     v3[i]*c1[3*nLevels+iIdx] +
                     v4[i]*c1[4*nLevels+iIdx] +
                     v5[i]*c1[5*nLevels+iIdx] +
                     v6[i]*c1[6*nLevels+iIdx] +
                     v7[i]*c1[7*nLevels+iIdx]);
  
  }

  return R_NilValue;

}',plugins=c("cpp11"),
   depends=c("BH","bigmemory"),
   includes=c('#include "bigmemory/BigMatrix.h"','#include "bigmemory/MatrixAccessor.hpp"')  )



#####

GLM.INTERCEPT <- 1
GLM.CHIP.BIN1 <- 2
GLM.CHIP.BIN2 <- 3

DTABLE.BIN1 <- 1
DTABLE.BIN2 <- 2
DTABLE.IDX.DIFF <- 3
DTABLE.PETCOUNT <- 4
DTABLE.DISTANCE <- 5
DTABLE.ADJIDX <- 6

#####


chain.process <- function(h5datafile,data,model.init,options,iterations,prune,pb=NULL,full.return=F) {
  trianglesets <- data$trianglesets
  nsites <- data$nsites
  tmidx <- data$trianglemaxidx
  
  if(!is.big.matrix(model.init$paramd)) {
    paramd <- attach.big.matrix(model.init$paramd,path=options$outputDir)
  } else {
    paramd <- model.init$paramd
  }
  
  if(!is.big.matrix(model.init$parami)) {
    parami <- attach.big.matrix(model.init$parami,path=options$outputDir)
  } else {
    parami <- model.init$parami
  }


  alphaparam <- paramd[,4]
  betaparam <- paramd[,5]
  pp <- paramd[,3]
  lambda0 <- paramd[,1]
  lambda1 <- paramd[,2]
  bothinpeak <- parami[,1]
  p1lambda <- model.init$p1lambda
  p2lambda <- model.init$p2lambda
  
  rm(model.init)
  
  chip1.signal <- as.integer(floor(colsum(data$dataGLM,GLM.CHIP.BIN1)))
  chip2.signal <- as.integer(floor(colsum(data$dataGLM,GLM.CHIP.BIN2)))
  
  zm <- as.integer(rep(0,nsites))
  binpc <- as.integer(rep(0,nsites))
  vz <- as.integer(rep(0,nsites))
  expected <- as.numeric(rep(0,nsites))
  
  coef0 <- coef1 <- matrix(0,tmidx,ncol(data$dataGLM))
  poissonfam <- poisson() ### only need one copy of this
  
  h5write(alphaparam,h5datafile,"alphaparam")
  h5write(betaparam,h5datafile,"betaparam")

  for( i in 1:iterations ) { 
    ### WARNING: many of the functions in this loop operate via pass by reference semantics, 
    ### so the values may be different post function call
    lIdx <- updateGroupStructure(data$dataM,
                              lambda0,
                              lambda1,
                              pp, #Will be modified
                              zm, #Will be modified
                              vz, #will be modified
                              bothinpeak,
                              alphaparam,
                              betaparam)
    
    ng0 <- length(lIdx$wisgroup0)
    ng1 <- length(lIdx$wisgroup1)
    

    if(!options$miniModel) h5write(pp,h5datafile,paste0("pp/ELT",i))

    
    updatePeakDist(bothinpeak, # will be modified
                   binpc, # will be modified
                   data$dataGLM,
                   p1lambda,
                   p2lambda,
                   options$minBinsWithPeak) 

    p1lambda <- rgamma(1,1+chip1.signal,1+data$nsites)
    p2lambda <- rgamma(1,1+chip2.signal,1+data$nsites)
    
    h5write(p1lambda,h5datafile,paste0("p1lambda/ELT",i))
    h5write(p2lambda,h5datafile,paste0("p2lambda/ELT",i))

    
    #### update GLM
    
    idx <- 0

    for( s in trianglesets ) {
      if(length(s)==0) {
        idx <- idx+1
        next
      }
     # fullm <- data$dataGLM[s,]
      
      w <- intersect(lIdx$wisgroup0,s)
      if( length(w)>1) {
        #cf0 <- glm.fit(data$dataGLM[w,],data$dataM[w,DTABLE.PETCOUNT],start=coef0[[idx+1]],family=poissonfam)
        cf0 <- glm.fit(data$dataGLM[w,],data$dataM[w,DTABLE.PETCOUNT],family=poissonfam,control=list(epsilon=options$glm.epsilon)) #start=coef0[idx+1,],
        
        cv <- coef(cf0)
        if(any(n <- is.na(cv))) cv[n] <- 0
        #lambda0[s] <- exp((fullm %*% cbind(cv))[,1])
        #coef0[[idx+1]] <- cv
        coef0[idx+1,] <- cv
      } else {
        ### use the previous row if there are too few data points, set to max value if this happens on the first row
        #coef0[[idx+1]] <- if(idx>0) coef0[[idx]] else coef0[[1]]
        coef0[idx+1,] <- if(idx>0) coef0[idx,] else coef0[1,]
      }
        
      w <- intersect(lIdx$wisgroup1,s)
      if( length(w)>1) {
        #cf1 <- glm.fit(data$dataGLM[w,],data$dataM[w,DTABLE.PETCOUNT],start=coef1[[idx+1]],family=poissonfam)
        cf1 <- glm.fit(data$dataGLM[w,],data$dataM[w,DTABLE.PETCOUNT],family=poissonfam,control=list(epsilon=options$glm.epsilon)) #start=coef1[idx+1,],
        cv <- as.vector(coef(cf1))
        if(any(n <- is.na(cv))) cv[n] <- 0
        #lambda1[s] <- exp((fullm %*% cbind(cv))[,1])
        #coef1[[idx+1]] <- cv
        coef1[idx+1,] <- cv
        
      } else {
        #coef1[[idx+1]] <- if(idx>0) coef1[[idx]] else coef1[[1]]
        coef1[idx+1,] <- if(idx>0) coef1[idx,] else coef1[1,]
        
      }
      
      #rm(fullm)
      idx <- idx + 1
    }
    
    updateLambdas(lambda0,lambda1,data$dataGLM,data$dataM,coef0,coef1)

    h5write(coef0,h5datafile,paste0("lambda0coef/ELT",i))
    h5write(coef1,h5datafile,paste0("lambda1coef/ELT",i))
    

    updateExpected(lambda0,lambda1,pp,expected)

    if(!options$miniModel)
    {
      h5write(lambda0,h5datafile,paste0("lambda0/ELT",i))
      h5write(lambda1,h5datafile,paste0("lambda1/ELT",i))
    }
    
    if(!is.null(pb)) setTxtProgressBar(pb,i/iterations)
  }
  
  h5write(zm,h5datafile,"zm")
  h5write(binpc,h5datafile,"bothinpeakcount")
  h5write(expected,h5datafile,"expectedsum")
  
  if(full.return) {
    paramd[,1] <- lambda0
    paramd[,2] <- lambda1
    paramd[,3] <- pp
    paramd[,4] <- alphaparam
    paramd[,5] <- betaparam
    
    if( is.filebacked(paramd)) flush(paramd)
    
    parami[,1] <- bothinpeak
    
    if( is.filebacked(parami)) flush(parami)
    
    list(#zm=zm,
         paramd=paramd,
         parami=parami,
         p1lambda=p1lambda,
         p2lambda=p2lambda)
  } else {
    list(zm=zm,expected=expected)
  }
}

setup.h5.model.output.file <- function(fname) {
  if(file.exists(fname)) stop(paste("H5 file",fname,"already exists"))
  
  h5createFile(fname)
  
  h5createGroup(fname,"lambda0coef")
  h5createGroup(fname,"lambda1coef")
  h5createGroup(fname,"lambda0")
  h5createGroup(fname,"lambda1")
  h5createGroup(fname,"p1lambda")
  h5createGroup(fname,"p2lambda")
  h5createGroup(fname,"pp")

}

estimate.global.bayesian.mixture <- function(data,
                                             hdf5.file.prefix,
                                             options) {
  
  paramd.bin.file <- paste0(options$filePrefix,"-double-param.bin")
  paramd.desc.file <- paste0(options$filePrefix,"-double-param.desc")
  
  parami.bin.file <- paste0(options$filePrefix,"-integer-param.bin")
  parami.desc.file <- paste0(options$filePrefix,"-integer-param.desc")
  

  bm.file <- paste0(hdf5.file.prefix,"-bm.txt")
  data.file <- paste0(hdf5.file.prefix,"-data.h5")
  burnin.file <- paste0(hdf5.file.prefix,"-burnin.h5")
  full.file <- paste0(hdf5.file.prefix,"-full.h5")
  parallel.files <- sapply(1:options$parallelCores,function(i) paste0(hdf5.file.prefix,paste0("-thread",i,".h5")))
  
  setup.h5.model.output.file(burnin.file)
  
  model <- list()
  dataM <- data$dataM

  interchromosomal <- is.na(dataM[,DTABLE.DISTANCE])
  intfreq <- table(c(dataM[,DTABLE.BIN1],dataM[,DTABLE.BIN2]),rep(dataM[,DTABLE.PETCOUNT],2))
  
  
  intn <- as.integer(colnames(intfreq))
  intrn <- as.integer(rownames(intfreq))
  idx1 <- match(dataM[,DTABLE.BIN1],intrn)
  idx2 <- match(dataM[,DTABLE.BIN2],intrn)
  idxl <- lapply(1:length(idx1),function(i) c(idx1[i],idx2[i]))
  idxd <- as.vector(dataM[,DTABLE.IDX.DIFF])
  
  if(options$mergeRowsLessThan > 0 ) {
    cat(paste("Merging rows with less than",options$mergeRowsLessThan,"data points\n"))
    midx <- max(idxd[!interchromosomal])
    
    nTable <- table(idxd)
    
    for( i in midx:1 ) {
      n <- nTable[as.character(i)]
      if( is.na(n) || n == 0 ) next;
      if( n >= options$mergeRowsLessThan ) break
    }

    if(i==1)
      warning("After merging, all the rows had less than the given threshold of points, so everything was collapsed into a single bin",immediate. = T)
        
    w <- which(!interchromosomal & idxd>i)
    
    idxd[w] <- i
    
    cat(paste("There are now",i,"out of",midx,"rows left\n"))
  }
  
  if(options$mergeRowsAfter > 0 ) {
    w <- which(!interchromosomal & idxd>options$mergeRowsAfter)
    
    idxd[w] <- options$mergeRowsAfter
  }
  
  if(any(interchromosomal)) idxd[interchromosomal] <- max(idxd[!interchromosomal])+1 ### make an interchromosomal group
  
  cat("Building putative interactions...\n")
  
  

  data[["nsites"]] <- as.numeric(nrow(dataM))
  data[["totcounts"]] <- colsum(dataM,DTABLE.PETCOUNT)
  data[["trianglemaxidx"]] <- max(idxd)
  data[["trianglesets"]] <- lapply(1:data$trianglemaxidx,function(i) which(idxd==i))
  
  data$dataM[,DTABLE.ADJIDX] <- as.integer(idxd)
  
  
  options[["minintdist"]] <- min(data$dataM[,DTABLE.DISTANCE][!interchromosomal])
  
  options[["origPruning"]] <- options$pruning
  if(is.null(options$pruning) || options$pruning < 1) {
    options$pruning <- options$iterations+1 # will never prune
  } 
  


  cat("Building mixture probability priors...\n")
  
  
  
  u.counts <- unique(data$dataM[,DTABLE.PETCOUNT]) ### pet counts
  names(u.counts) <- u.counts
  u.counts.b <- lapply(u.counts,function(p) intn < p)
  intfreqb <- intfreq>0

  countm <- t(mapply(function(cidx,iv,id){
    
    if(id<2) return(c(0,8)) ## downweigh those right on the diagonal
    
    cb <- u.counts.b[[cidx]]
    intf <- intfreqb[iv,]
    
    a <- sum2(as.integer(intf & cb),mode='integer')
    b <- sum2(as.integer(intf & !cb),mode='integer')
    
    c(a,b)
  },as.list(match(data$dataM[,DTABLE.PETCOUNT],as.integer(names(u.counts.b)))),as.list(idxl),idxd,SIMPLIFY=T))


  
  alphaparam <- 1+countm[,1]
  betaparam <- 1+countm[,2]


  cat("Randomly generating initial parameter values...\n")
  
  pp <- rbeta(data$nsites,alphaparam,betaparam)
  
  ### set initial lambdas by sampling the PET counts
  cl0 <- sample(data$dataM[,DTABLE.PETCOUNT],data$nsites,replace=T)
  cl1 <- sample(data$dataM[,DTABLE.PETCOUNT],data$nsites,replace=T)

  swapLambdas(cl0,cl1)
  

  
  model[["p1lambda"]] <- rgamma(1,1,1)
  model[["p2lambda"]] <- rgamma(1,1,1)
  
  bothinpeak <- as.integer(rbinom(data$nsites,1,.5))
  
  model[["paramd"]] <- if( options$filebacked ) as.big.matrix(cbind(cl0,cl1,pp,alphaparam,betaparam),
                          type='double',backingpath=options$outputDir,backingfile=paramd.bin.file,descriptor=paramd.desc.file)
                       else as.big.matrix(cbind(cl0,cl1,pp,alphaparam,betaparam), type='double',shared=F)
  model[["parami"]] <- if(options$filebacked) as.big.matrix(cbind(bothinpeak),
                          type='integer',backingpath=options$outputDir,backingfile=parami.bin.file,descriptor=parami.desc.file)
                      else as.big.matrix(cbind(bothinpeak), type='integer',shared=F)

  
  rm(cl0,cl1,bothinpeak,pp,countm,u.counts,u.counts.b,intfreqb,intn,idx1,idx2,idxl,intfreq,alphaparam,betaparam,interchromosomal)
  gc()
  
  cat("Burn-in...\n")
  
  if(options$showProgress) pb <- txtProgressBar()
  
  burnin.l <- chain.process(burnin.file,data,model,options,options$burnin,options$pruning,pb,full.return = T)

  if(options$showProgress) close(pb)
  
  rm(model)
  gc()

  cat("MCMC...\n")
  
  
  if(options$parallel) {
    cat("(Note: there is no progress bar displayed when running in parallel mode\n")
    Niter <- floor(options$iterations/options$parallelCores)
    
    lapply(parallel.files,setup.h5.model.output.file)
    
    il <- mclapply(parallel.files,
             chain.process,
             data=data,
             model=burnin.l,
             options=options,
             iterations=Niter,
             prune=options$pruning,
             full.return=F,
             mc.preschedule=F,
             mc.cores=options$parallelCores)
    
    if(any(e <- sapply(il,function(r) !any(names(r) == "zm")))) {
      print(il)
      stop("There was an error during parallelization")
    }
    
    zmat <- do.call(cbind,lapply(il,function(l) l$zm))
    emat <- rowSums(do.call(cbind,lapply(il,function(l) l$expected)))/options$iterations
    
    ret <- list(zm=rowSums(zmat),expected=emat)
  } else {
    if(options$showProgress) pb <- txtProgressBar()
    
    setup.h5.model.output.file(full.file)

    ret <- chain.process(full.file,data,burnin.l,options,options$iterations,options$pruning,pb,full.return = F)
    ret$expected <- ret$expected/options$iterations
    
    if(options$showProgress) close(pb)
  }
  
  
  dl <- list(prob=ret$zm/options$iterations,expected=ret$expected)

  dl
}


