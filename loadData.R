library(dnar)
library(parallel)
source('functions.R')
source('16s/functions.R')
source('readSamples.R')

#fastqs<-list.files('data/','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
#fastqs<-fastqs[!grepl('Undetermined',fastqs)]
#primers<-sub('.*(matK|rbcL).*_R([0-9]+)_.*','\\1\\2',basename(fastqs))
primers<-c('matK','rbcL')

nReadCut<-5000
swarmData<-lapply(unique(primers),function(ii){
  message(ii)
  outMat<-sprintf('work/swarmPair/%s.Rdat',ii)
  outFa<-sprintf('work/swarmPair/%s.fa.gz',ii)
  outTaxa<-sprintf('work/swarmPair/%s_taxa.csv',ii)
  outHits<-sprintf('work/swarmPair/%s_allHits.csv',ii)
  outAlign<-sprintf('work/swarmPair/%s_align.fa.gz',ii)
  if(!file.exists(outMat)|!file.exists(outFa))source('makeOtus.R')
  if(!file.exists(outTaxa)|!file.exists(outHits))source('parseBlast.R')
  tmp<-environment()
  load(outMat,tmp)
  swarmOtus<-get('swarmOtus',tmp)
  rownames(swarmOtus)<-sub('^([A-Z]+[0-9]+).*','\\1',rownames(swarmOtus))
  rm(tmp)
  swarmSeqs<-read.fa(outFa)
  rownames(swarmSeqs)<-swarmSeqs$name
  swarmTaxa<-read.csv(outTaxa,row.names=1,stringsAsFactors=FALSE)
  swarmTaxa$seq<-swarmSeqs[rownames(swarmTaxa),'seq']
  swarmTaxa$bestAnnot<-apply(swarmTaxa[,c('phylum','class','order','family','genus','species')],1,function(x)ifelse(is.na(x),x,sprintf('%s_%s',substring(names(x),1,1),x))[max(c(1,which(!is.na(x))))])
  swarmTaxa$bestId<-ave(swarmTaxa$bestAnnot,naReplace(swarmTaxa$best,'__NAFILLER__'),FUN=function(x){sprintf('%s #%d',ifelse(is.na(x),'Unknown',x),1:length(x))})
  #trashing singletons
  swarmOtus<-swarmOtus[,apply(swarmOtus,2,sum)>1]
  props<-t(apply(swarmOtus,1,function(xx)xx/sum(xx)))
  nReads<-apply(swarmOtus,1,sum)
  isEnough<-nReads>nReadCut
  subsampledOtus<-cacheOperation(sprintf('work/%s_rarefyOtus.Rdat',ii),apply,swarmOtus,2,function(xx)if(sum(xx)<nReadCut)rep(NA,length(xx)) else rarefyCounts(xx,nReadCut))
  return(list('otus'=swarmOtus,'props'=props,'rare'=subsampledOtus,'taxa'=swarmTaxa,'nReads'=nReads,'isEnough'=isEnough))
})
names(swarmData)<-unique(primers)
