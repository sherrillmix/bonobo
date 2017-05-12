library(dnar)
library(parallel)
library(dnaplotr)
source('functions.R')
source('readSamples.R')

fastqs<-list.files('data/','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined',fastqs)]
primers<-sub('.*(matK|rbcL).*_R([0-9]+)_.*','\\1\\2',basename(fastqs))

for(ii in unique(primers)){
  message(ii)
  outMat<-sprintf('work/swarm/%s.Rdat',ii)
  outFa<-sprintf('work/swarm/%s.fa.gz',ii)
  if(all(file.exists(outMat,outFa))){
    message('Already processed. Skipping')
    next()
  }
  thisPrimer<-primerSeqs[[sub('[12]$','',tolower(ii))]][as.numeric(substring(ii,nchar(ii)))]
  thisFiles<-fastqs[primers==ii]
  reads<-mclapply(thisFiles,function(xx){library(dnar);cat('.');read.fastq(xx)$seq},mc.cores=6,mc.preschedule=FALSE)
  trimReads<-lapply(reads,function(xx)substring(xx[!grepl('[^ACTG]',xx)],nchar(thisPrimer)+1))
  samples<-rep(basename(thisFiles),sapply(trimReads,length))
  otus<-runSwarm(unlist(trimReads),'~/installs/swarm/swarm',swarmArgs='-f -t 30')
  swarmOtus<-as.data.frame.matrix(table(samples,otus[['otus']]))
  write.fa(otus[['seqs']]$name,otus[['seqs']]$seq,outFa)
  save(swarmOtus,file=outMat)
}

primerBases<-sub('[0-9]$','',primers)
for(primerBase in unique(primerBases)){
  message(primerBase)
  outMat<-sprintf('work/swarmPair/%s.Rdat',primerBase)
  outFa<-sprintf('work/swarmPair/%s.fa.gz',primerBase)
  outAlign<-sprintf('work/swarmPair/%s_align.fa.gz',primerBase)
  if(all(file.exists(outMat,outFa,outAlign))){
    message('Already processed. Skipping')
    next()
  }
  trimReads<-lapply(sprintf('%s%d',primerBase,1:2),function(ii){
    thisPrimer<-primerSeqs[[sub('[12]$','',tolower(ii))]][as.numeric(substring(ii,nchar(ii)))]
    thisFiles<-fastqs[primers==ii]
    reads<-mclapply(thisFiles,function(xx){library(dnar);cat('.');read.fastq(xx)},mc.cores=10,mc.preschedule=FALSE)
    if(mean(unlist(lapply(reads,function(xx)substring(xx$seq,1,nchar(thisPrimer))))%in% expandAmbiguous(thisPrimer)[[1]])<.75)stop(simpleError('Expected primer does not match read start'))
    trimReads<-lapply(reads,function(xx){xx$seq<-substring(xx$seq,nchar(thisPrimer)+1);return(xx)})
    names(trimReads)<-thisFiles
    return(trimReads)
  })
  if(any(names(trimReads[[1]])!=sub('_R2_','_R1_',names(trimReads[[2]]))))stop('Read 1 vs reads 2 file name mismatch')
  readCounts<-as.data.frame(sapply(trimReads,function(xx)sapply(xx,nrow)))
  colnames(readCounts)<-c('raw1','raw2')
  #concatenate high quality reads
  trimReads<-mcmapply(function(left,right,...){
    cat('.')
    if(any(sub(' .*$','',left$name)!=sub(' .*$','',right$name)))stop('Read 1 vs read 2 name mismatch')
    #last base always low qual so ignore
    q1<-sapply(qualToInts(substring(left$qual,1,nchar(left$qual)-1)),function(xx)sum(10^(-xx/10)))
    q2<-sapply(qualToInts(substring(right$qual,1,nchar(right$qual)-1)),function(xx)sum(10^(-xx/10)))
    #less than 1 expected error in both reads
    selector<-q1<1&q2<1 & !grepl('[^ACTG]',left$seq)&!grepl('[^ACTG]',right$seq)
    seqs<-paste(left[selector,'seq'],revComp(right[selector,'seq']),sep='')
    return(seqs)
  },trimReads[[1]],trimReads[[2]],mc.cores=5,SIMPLIFY=FALSE)
  readCounts$filter<-sapply(trimReads,length)
  write.csv(readCounts,'work/swarmPair/%s_counts.csv')
  samples<-rep(basename(names(trimReads)),sapply(trimReads,length))
  otus<-runSwarm(unlist(trimReads),'~/installs/swarm/swarm',swarmArgs='-f -t 40')
  swarmOtus<-as.data.frame.matrix(table(samples,otus[['otus']]))
  write.fa(otus[['seqs']]$name,otus[['seqs']]$seq,outFa)
  save(swarmOtus,file=outMat)
  load(outMat)
  tmpFile<-tempfile()
  otus<-list('seqs'=read.fa(outFa))
  #throwing out singletons for aligning
  swarmOtus<-swarmOtus[,apply(swarmOtus,2,sum)>1]
  seqs<-otus[['seqs']]
  rownames(seqs)<-seqs$name
  write.fa(colnames(swarmOtus),seqs[colnames(swarmOtus),'seq'],tmpFile)
  cmd<-sprintf('~/installs/mafft/bin/mafft --thread 50 %s|gzip>%s',tmpFile,outAlign)
  message(cmd)
  system(cmd)
  align<-read.fa(outAlign)
  png(sprintf('out/align_%s.png',primerBase),width=4000,height=4000,res=250);plotDNA(align$seq);dev.off()
}
