library(dnar)
library(parallel)
source('functions.R')
source('readSamples.R')

fastqs<-list.files('data/','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined',fastqs)]
primers<-sub('.*(matK|rbcL).*_R([0-9]+)_.*','\\1\\2',basename(fastqs))

swarmData<-lapply(unique(primers),function(ii){
  message(ii)
  outMat<-sprintf('work/swarm/%s.Rdat',ii)
  outFa<-sprintf('work/swarm/%s.fa.gz',ii)
  outTaxa<-sprintf('work/swarm/%s_taxa.csv',ii)
  outHits<-sprintf('work/swarm/%s_allHits.csv',ii)
  if(!file.exists(outMat)|!file.exists(outFa))source('makeOtus.R')
  if(!file.exists(outTaxa)|!file.exists(outHits))source('parseBlast.R')
  tmp<-environment()
  load(outMat,tmp)
  swarmOtus<-get('swarmOtus',tmp)
  rm(tmp)
  swarmSeqs<-read.fa(outFa)
  rownames(swarmSeqs)<-swarmSeqs$name
  swarmTaxa<-read.csv(outTaxa,row.names=1)
  swarmTaxa$seq<-swarmSeqs[rownames(swarmTaxa),'seq']
  return(list('otus'=swarmOtus,'taxa'=swarmTaxa))
})
names(swarmData)<-unique(primers)

nReads<-lapply(swarmData,function(xx)apply(xx[['otus']],1,sum))
isEnough<-sapply(nReads,function(xx)xx>5000)

