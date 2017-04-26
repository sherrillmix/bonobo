library(dnar)
library(parallel)
source('functions.R')
source('16s/functions.R')
source('readSamples.R')

fastqs<-list.files('data/','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined',fastqs)]
primers<-sub('.*(matK|rbcL).*_R([0-9]+)_.*','\\1\\2',basename(fastqs))

nReadCut<-5000
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
  swarmTaxa<-read.csv(outTaxa,row.names=1,stringsAsFactors=FALSE)
  swarmTaxa$seq<-swarmSeqs[rownames(swarmTaxa),'seq']
  swarmTaxa$bestId<-ave(swarmTaxa$best,naReplace(swarmTaxa$best,'__NAFILLER__'),FUN=function(x){sprintf('%s #%d',ifelse(is.na(x),'Unknown',x),1:length(x))})
  #trashing singletons
  swarmOtus<-swarmOtus[,apply(swarmOtus,2,sum)>1]
  props<-t(apply(swarmOtus,1,function(xx)xx/sum(xx)))
  nReads<-apply(swarmOtus,1,sum)
  isEnough<-nReads>nReadCut
  return(list('otus'=swarmOtus,'props'=props,'taxa'=swarmTaxa,'nReads'=nReads,'isEnough'=isEnough))
})
names(swarmData)<-unique(primers)


