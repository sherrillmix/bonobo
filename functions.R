
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
parseQiimeTaxa<-function(taxas,desiredTaxa=c('k','p','c','o','f','g','s'),concatLastTwo=TRUE){
  taxa<-strsplit(taxaRaw$taxa,'[;] ?')
  taxa<-do.call(rbind,lapply(taxa,function(xx){
    y<-xx[grep('^[a-z]__',xx)]
    names(y)<-sub('^([a-z])__.*','\\1',y)
    y<-sub('^[a-z]__','',y)
    y[y=='']<-NA
    out<-y[desiredTaxa]
    if(!is.na(out[length(desiredTaxa)]))out[length(desiredTaxa)]<-paste(out[length(desiredTaxa)-1:0],collapse=' ')
    return(out)
  }))
  return(as.data.frame(taxa,stringsAsFactors=FALSE))
}


