#library(dnar)
#fastqs<-list.files('data','_R[12]_.*fastq.gz$',recursive=TRUE,full.names=TRUE)
#nReads<-unlist(cleanMclapply(fastqs,function(xx){library(dnar);nrow(read.fastq(xx))},mc.cores=30))
#names(nReads)<-sub('_S[0-9]+_L[0-9]+_R([0-9]).*$','',basename(fastqs))

if(!exists('swarmData'))source('loadData.R')

readCounts<-sapply(swarmData,function(xx){
  readCounts<-apply(xx[['otus']],1,sum)
})
rawReadCounts<-do.call(rbind,cacheOperation('work/rawReadCounts.Rdat',mclapply,rownames(readCounts),function(xx){
  rbcl<-list.files('data',sprintf('%srbcL.*_R1_.*\\.fastq.gz',xx),recursive=TRUE,full.names=TRUE)
  matk<-list.files('data',sprintf('%smatK.*_R1_.*\\.fastq.gz',xx),recursive=TRUE,full.names=TRUE)
  rbclN<-sum(sapply(rbcl,function(xx)as.numeric(system(sprintf('zcat %s|wc -l',xx),intern=TRUE))/4))
  matkN<-sum(sapply(matk,function(xx)as.numeric(system(sprintf('zcat %s|wc -l',xx),intern=TRUE))/4))
  return(c('rawMatk'=matkN,'rawRbcl'=rbclN))
},mc.cores=5))
readCounts<-cbind(readCounts,rawReadCounts)
rbcl<-read.csv('work/swarmPair/rbcL_counts.csv',row.names=1)
matk<-read.csv('work/swarmPair/matK_counts.csv',row.names=1)
if(any(rbcl$raw1!=rbcl$raw2))stop('Mismatch between left and right')
if(any(matk$raw1!=matk$raw2))stop('Mismatch between left and right')
rbcl$base<-sub('rbcL.*$','',basename(rownames(rbcl)))
matk$base<-sub('matK.*$','',basename(rownames(matk)))
swarmN<-data.frame(
  'rawRbcl'=tapply(rbcl$raw1,rbcl$base,sum),
  'rawMatk'=tapply(matk$raw1,matk$base,sum),
  'filterRbcl'=tapply(rbcl$filter,rbcl$base,sum),
  'filterMatk'=tapply(matk$filter,matk$base,sum)
)
if(!any(rownames(readCounts) %in% rownames(swarmN)))stop('Missing names in swarm')
swarmN<-swarmN[rownames(readCounts),]
if(any(swarmN[,'rawRbcl']!=readCounts[,'rawRbcl']))stop('Mismatch in raws')
if(any(swarmN[,'rawMatk']!=readCounts[,'rawMatk']))stop('Mismatch in raws')
readCounts<-cbind(readCounts,swarmN[rownames(readCounts),c('filterRbcl','filterMatk')])

tableOrder<-readLines('tableOrder.csv')
if(any(!rownames(readCounts)[rownames(readCounts) %in% rownames(samples)] %in% tableOrder) || any(!tableOrder %in% rownames(readCounts)))stop('Table order problem')
write.csv(readCounts[tableOrder,c('rawMatk','filterMatk','matK','rawRbcl','filterRbcl','rbcL')],'out/readCounts.csv',)

bact<-read.csv('16s/out/16s_readCounts.csv',row.names=1)[tableOrder,,drop=FALSE]
readCounts<-readCounts[tableOrder,]
if(any(rownames(bact)!=rownames(readCounts)))stop('Name mismatch in counts')
readCounts$"16S"<-bact[,"X16s"]
write.csv(readCounts[tableOrder,c('matK','rbcL','16S')],'out/finalReadCounts.csv',)
