library(dnar)
library(parallel)
source('functions.R')

fastqs<-list.files('data','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined',fastqs)]
primers<-sub('.*(matK|rbcL).*_R([0-9]+)_.*','\\1\\2',basename(fastqs))


for(ii in unique(primers)){
  message(ii)
  thisPrimer<-primerSeqs[[sub('[12]$','',tolower(ii))]][as.numeric(substring(ii,nchar(ii)))]
  thisFiles<-fastqs[primers==ii]
  reads<-mclapply(thisFiles,function(xx){library(dnar);cat('.');read.fastq(xx)$seq},mc.cores=6,mc.preschedule=FALSE)
  #WARNING cutting reads 
  trimReads<-lapply(reads,function(xx)substring(xx[!grepl('[^ACTG]',xx)],nchar(thisPrimer)+1,200))
  samples<-rep(basename(thisFiles),sapply(trimReads,length))
  otus<-runSwarm(unlist(trimReads),'~/installs/swarm/swarm',swarmArgs='-f -t 30')
  browser()
  outDir<-sprintf('work/qiime_%s',ii)
  if(!dir.exists(outDir)){
    out<-runQiime(unlist(reads),storeDir=outDir)
  }
  #write.csv(cbind(samples,'otu'=otus[['otus']]),sprintf('work/swarm_%s.csv',ii))
  #write.csv(cbind(samples,otus[['otus']]),sprintf('work/swarm_%s.fa'))
}
