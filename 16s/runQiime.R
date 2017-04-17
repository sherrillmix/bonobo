library(parallel)
library(dnar)

source('../readSamples.R',chdir=TRUE)
source('../functions.R')
source('functions.R')
nRequiredReads<-15000

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
taxaRaw<-read.csv('work/qiimeOtus.taxa',stringsAsFactors=FALSE)
taxa<-parseQiimeTaxa(taxaRaw$taxa)
rownames(taxa)<-taxaRaw$name
taxa$best<-apply(taxa,1,function(x)ifelse(is.na(x),x,sprintf('%s_%s',names(x),x))[max(c(1,which(!is.na(x))))])
taxa$bestId<-ave(taxa$best,naReplace(taxa$best,'__NAFILLER__'),FUN=function(x){sprintf('%s #%d',ifelse(is.na(x),'Unknown',x),1:length(x))})
#filter chloroplast reads
isChloro<-taxa[rownames(otuTab),'c']=='Chloroplast'&!is.na(taxa[rownames(otuTab),'c'])
otuTab<-otuTab[!isChloro,]
otuSeq<-read.fa('work/qiimeOtus.fa.gz')
rownames(otuSeq)<-otuSeq$name
taxa$seq<-otuSeq[rownames(taxa),'seq']

samples$name<-sapply(samples$Code,function(x)colnames(otuTab)[grep(sprintf('_%s$',x),colnames(otuTab))])
rownames(samples)<-samples$name
#make sure same order as otu table
if(any(!samples$name %in% colnames(otuTab)))stop('Sample missing from OTUs')
samples<-samples[colnames(otuTab)[colnames(otuTab) %in% samples$name],]
samples$isEnough<-apply(otuTab[,samples$name],2,sum)>nRequiredReads
head(sort(apply(otuTab[,samples$name],2,sum)))

otuProp<-apply(otuTab,2,function(x)x/sum(x))

readCounts<-apply(otuTab,2,sum)
mean(readCounts[samples[samples$isEnough,'name']])
