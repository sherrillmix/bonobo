library(parallel)
library(dnar)

source('../readSamples.R',chdir=TRUE)

runQiime<-function(seqs,storeDir=NULL){
  if(any(is.na(seqs)))stop(simpleError('NAs in seqs'))
  seqIds<-1:length(seqs)
  readDir<-tempfile()
  dir.create(readDir)
  readFile<-file.path(readDir,'XXX.fa')
  outDir<-tempfile()
  seqNames<-sprintf('XXX_%d',seqIds)
  write.fa(seqNames,seqs,readFile)
  #miniconda doesn't like sh so need to use bash
  cmd<-sprintf('echo "source activate qiime1; pick_de_novo_otus.py --input %s --output %s --parallel --jobs_to_start 32 --force"|bash',readFile,outDir)
  message(cmd)
  exit<-system(cmd)
  if(exit!=0)stop(simpleError('Problem running qiime'))
  #get otu assignments
  assigns<-strsplit(readLines(file.path(outDir,'uclust_picked_otus/XXX_otus.txt')),'\t')
  names(assigns)<-sapply(assigns,'[[',1)
  assigns<-lapply(assigns,'[',-1)
  otus<-rep(names(assigns),sapply(assigns,length))
  names(otus)<-unlist(assigns)
  out<-otus[seqNames]
  #get taxa assignments
  taxa<-strsplit(readLines(file.path(outDir,'uclust_assigned_taxonomy/XXX_rep_set_tax_assignments.txt')),'\t')
  names(taxa)<-sapply(taxa,'[[',1)
  taxa<-sapply(taxa,'[[',2)
  #get sequences
  seqs<-read.fa(file.path(outDir,'pynast_aligned_seqs/XXX_rep_set_aligned_pfiltered.fasta'))
  seqs<-structure(seqs$seq,.Names=seqs$name)
  if(!is.null(storeDir)){
    #avoid problems with cross filesystem copies
    tmpDir<-tempfile(tmpdir='.')
    dir.create(tmpDir)
    file.copy(outDir,tmpDir)
    file.rename(file.path(tmpDir,basename(outDir)),storeDir)
  }
  return(list('otus'=out,'seqs'=seqs,'taxa'=taxa))
}

fastqs<-list.files('data/joined','.fastq.gz',full.name=TRUE)
#just get sequence to reduce memory
allSeq<-lapply(fastqs,function(xx)read.fastq(xx)$seq)
if(!file.exists('work/qiime')){
  out<-runQiime(unlist(allSeq),storeDir='work/qiime')
  withAs(outFile=gzfile('work/qiimeOtuIds.csv.gz'),write.csv(data.frame('file'=rep(sub('_Feces_1_1.fastq.gz','',basename(fastqs)),sapply(allSeq,length)),'otu'=out[['otus']],stringsAsFactors=FALSE),outFile,row.names=FALSE))
  write.fa(names(out[['seqs']]),out[['seqs']],'work/qiimeOtus.fa.gz')
  write.csv(data.frame('name'=names(out[['taxa']]),'taxa'=out[['taxa']],stringsAsFactors=FALSE),'work/qiimeOtus.taxa',row.names=FALSE)
}
otus<-read.csv('work/qiimeOtuIds.csv.gz',stringsAsFactors=FALSE)
otuTab<-as.data.frame.matrix(table(otus$otu,otus$file))
otuTab<-otuTab[apply(otuTab,1,sum)>3,]

samples$name<-sapply(samples$Code,function(x)colnames(otuTab)[grep(sprintf('_%s$',x),colnames(otuTab))])
rownames(samples)<-samples$name
#make sure same order as otu table
if(any(!samples$name %in% colnames(otuTab)))stop('Sample missing from OTUs')
samples<-samples[colnames(otuTab)[colnames(otuTab) %in% samples$name],]

