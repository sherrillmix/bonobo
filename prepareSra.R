source('readSamples.R')

fastqs<-list.files('data','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined',fastqs)]
primers<-sub('.*(matK|rbcL).*_R([0-9]+)_.*','\\1',basename(fastqs))
ends<-sub('.*(matK|rbcL).*_R([0-9]+)_.*','\\2',basename(fastqs))
sampleNames<-sub('^(.*)(matK|rbcL).*_R([0-9]+)_.*','\\1',basename(fastqs))
out<-data.frame('file'='__NO__','sample'='__NO__','primer'='__NO__','end'=1,stringsAsFactors=FALSE)[0,]
for(ii in unique(sampleNames)){
  if(!ii %in% rownames(samples))next()
  for(jj in unique(primers)){
    message(ii,jj)
    seqs<-lapply(c('1','2'),function(kk){
      thisFiles<-fastqs[ends==kk&primers==jj&sampleNames==ii]
      do.call(rbind,lapply(thisFiles,dnar::read.fastq))
    })
    if(nrow(seqs[[1]])!=nrow(seqs[[2]]))stop('Left right mismatch')
    for(kk in 1:2){
      outFastq<-file.path('work/sra/',sprintf('%s_%s_%d.fastq.gz',ii,jj,kk))
      message(outFastq)
      write.fastq(seqs[[kk]]$name,seqs[[kk]]$seq,seqs[[kk]]$qual,outFastq)
      out<-rbind(out,data.frame('file'=outFastq,'sample'=ii,'primer'=jj,'end'=kk,stringsAsFactors=FALSE))
    }
  }
}

fastqs<-list.files('16s/data','_[12].fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined|joined',fastqs)]
ends<-sub('^.*_([12])\\.fastq.gz','\\1',basename(fastqs))
sampleNames<-sub('^.*_([^_]+)_Feces_.*\\.fastq.gz','\\1',basename(fastqs))
for(ii in unique(sampleNames)){
  if(!ii %in% rownames(samples))next()
  message(ii)
  seqs<-lapply(c('1','2'),function(kk){
    thisFiles<-fastqs[ends==kk&sampleNames==ii]
    do.call(rbind,lapply(thisFiles,dnar::read.fastq))
  })
  if(nrow(seqs[[1]])!=nrow(seqs[[2]]))stop('Left right mismatch')
  for(kk in 1:2){
    outFastq<-file.path('work/sra/',sprintf('%s_%s_%d.fastq.gz',ii,'16s',kk))
    message(outFastq)
    write.fastq(seqs[[kk]]$name,seqs[[kk]]$seq,seqs[[kk]]$qual,outFastq)
    out<-rbind(out,data.frame('file'=outFastq,'sample'=ii,'primer'='16s','end'=kk,stringsAsFactors=FALSE))
  }
}


sraOut<-read.table('Metagenome.environmental.1.0.tsv',header=TRUE,check.names=FALSE)
mandatory<-grepl('\\*',colnames(sraOut))
colnames(sraOut)<-sub('\\*','',colnames(sraOut))
out$organism<-'feces metagenome'
out$host<-sub('P\\.t\\.','Pan troglodytes',samples[out$sample,'Species'])
out$lat_lon<-paste(samples[out$sample,'cleanLat'],samples[out$sample,'cleanLon'])
out$collection_date<-as.character(samples[out$sample,'rDate'])
out$sample_name<-out$sample
out$geo_loc_name<-sprintf('Democratic Republic of the Congo: %s',samples[out$sample,'area2'])
out$isolation_source<-'Feces'
out$bioproject_accession<-'PRJNA389566'
out$samp_mat_process<-sprintf('%s amplification',out$primer)
out$label<-out$sample
out[out$lat_lon==' ','lat_lon']<-'not collected'
write.table(out[out$end==1&out$primer=='16s',c('sample_name','bioproject_accession','organism','isolation_source','collection_date','geo_loc_name','lat_lon','label','host')],'work/sra/sraSampleDeposit.tsv',sep='\t',row.names=FALSE)


#downloaded from SRA after loading samples
sraDeposit<-read.table('work/sra/sraSampleIds.tsv',sep='\t',header=TRUE,stringsAsFactors=FALSE)
rownames(sraDeposit)<-sraDeposit$sample_name
out$biosample_accession<-sraDeposit[out$sample,'accession']
out$library_ID<-paste(out$sample,out$primer,sep='_')
out$title<-sprintf('%s amplification of %s feces: %s',out$primer,out$host,out$sample)
out$platform<-'Illumina'
out$instrument_model<-'Illumina MiSeq'
out$library_strategy<-'AMPLICON'
out$library_source<-'METAGENOMIC'
out$library_selection<-'PCR'
out$library_layout<-'paired'
out$design_description<-sprintf('DNA was extracted from a %s fecal sample and PCR amplified using %s specific primers',out$host,out$primer)
out$filetype<-'fastq'
out$filename2<-basename(ave(out$file,out$title,FUN=function(x)x[2]))
out$filename<-basename(ave(out$file,out$title,FUN=function(x)x[1]))
write.table(out[out$end==1,c('bioproject_accession','biosample_accession','library_ID','title','library_strategy','library_source','library_selection','library_layout','platform','instrument_model','design_description','filetype','filename','filename2')],'work/sra/sraReadDeposit.tsv',row.names=FALSE,sep='\t')
