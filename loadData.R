library(dnar)
library(parallel)
source('functions.R')
source('readSamples.R')

fastqs<-list.files('data/','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined',fastqs)]
primers<-sub('.*(matK|rbcL).*_R([0-9]+)_.*','\\1\\2',basename(fastqs))

for(ii in unique(primers)){
  message(ii)
  outMat<-sprintf('work/swarm/%s.Rdat',ii)
  outFa<-sprintf('work/swarm/%s.fa.gz',ii)
  outTaxa<-sprintf('work/swarm/%s_taxa.csv',ii)
  outHits<-sprintf('work/swarm/%s_allHits.csv',ii)
  if(!file.exists(outMat)|!file.exists(outFa))source('makeOtus.R')
  if(!file.exists(outTaxa)|!file.exists(outHits))source('parseBlast.R')
}
