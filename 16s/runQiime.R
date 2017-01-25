library(parallel)
library(dnar)


runQiime<-function(seqs){
  if(any(is.na(seqs)))stop(simpleError('NAs in seqs'))
  seqIds<-1:length(seqs)
  readDir<-tempfile()
  dir.create(readDir)
  readFile<-file.path(readDir,'XXX.fa')
  outDir<-tempfile()
  seqNames<-sprintf('XXX_%d',seqIds)
  write.fa(seqNames,seqs,readFile)
  #miniconda doesn't like sh so need to use bash
  cmd<-sprintf('echo "source activate qiime1; pick_de_novo_otus.py --input %s --output %s --parallel --jobs_to_start 16 --force"|bash',readFile,outDir)
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
  return(list('otus'=out,'seqs'=seqs,'taxa'=taxa))
}

fastqs<-list.files('data/joined','.fastq.gz',full.name=TRUE)
for(fastq in fastqs){
  message(fastq)
  outFiles<-sprintf(c('%s.otu','%s.otuSeq','%s.taxa'),sprintf('work/%s',sub('.fastq.gz$','',basename(fastq))))
  seqs<-read.fastq(fastq)
  out<-runQiime(seqs$seq)
  writeLines(out[['otus']],outFiles[1])
  write.fa(names(out[['seqs']]),out[[2]],outFiles[2])
  write.csv(out[['taxa']],outFiles[3])
}
