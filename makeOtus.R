library(dnar)
library(parallel)

fastqs<-list.files('data','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined',fastqs)]
primers<-sub('.*(matK|rbcL).*_R([0-9]+)_.*','\\1\\2',basename(fastqs))

runSwarm<-function(seqs,swarmBin='swarm',swarmArgs='-f'){
  if(any(is.na(seqs)))stop(simpleError('NAs in seqs'))
  seqIds<-as.numeric(as.factor(seqs))
  seqCounts<-ave(seqIds,seqIds,FUN=length)
  seqNames<-sprintf('%08d_%d',seqIds,seqCounts)
  readFile<-tempfile()
  outFile<-tempfile()
  seqFile<-tempfile()
  uniqSelector<-!duplicated(seqs)
  write.fa(seqNames[uniqSelector],seqs[uniqSelector],readFile)
  cmd<-sprintf('%s %s %s -o %s -w %s',swarmBin,swarmArgs,readFile,outFile,seqFile)
  system(cmd)
  swarm<-readLines(outFile)
  otuSplit<-strsplit(swarm,' ')
  otus<-rep(1:length(otuSplit),sapply(otuSplit,length))
  names(otus)<-unlist(otuSplit)  
  out<-otus[seqNames]
  seedSeqs<-read.fa(seqFile)
  return(list('otus'=out,'seqs'=seedSeqs))
}

for(ii in unique(primers)){
  message(ii)
  thisFiles<-fastqs[primers==ii]
  reads<-cleanMclapply(thisFiles,function(xx){library(dnar);cat('.');read.fastq(xx)$seq},mc.cores=5)
  reads<-lapply(reads,function(xx)xx[!grepl('[^ACTG]',xx)])
  samples<-rep(basename(thisFiles),sapply(reads,length))
  otus<-runSwarm(unlist(reads),'~/installs/swarm/swarm')
  browser()
  #write.csv(cbind(samples,'otu'=otus[['otus']]),sprintf('work/swarm_%s.csv',ii))
  #write.csv(cbind(samples,otus[['otus']]),sprintf('work/swarm_%s.fa'))
}
