library(dnar)
library(dnaplotr)
library(levenR)
library(parallel)


fastqs<-list.files('data/','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined',fastqs)]
allQs<-mclapply(fastqs,function(xx){
  library(dnar)
  seqs<-read.fastq(xx)
  if(nrow(seqs)<500)return(NULL)
  quals<-do.call(rbind,qualToInts(seqs$qual))
  return(apply(quals,2,mean))
},mc.cores=5)
names(allQs)<-fastqs

fastqs16s<-list.files('16s/data/','_[12]\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
allQs16s<-cleanMclapply(fastqs16s,function(xx){
  library(dnar)
  seqs<-read.fastq(xx)
  if(nrow(seqs)<500)return(NULL)
  quals<-do.call(rbind,qualToInts(seqs$qual))
  return(apply(quals,2,mean))
},mc.cores=10)
names(allQs16s)<-fastqs16s

r1s<-allQs[grepl('_R1_',names(allQs))&!grepl('Weimin_Plants_3_6_17',names(allQs))]
r2s<-allQs[grepl('_R2_',names(allQs))&!grepl('Weimin_Plants_3_6_17',names(allQs))]
r1s16s<-allQs16s[!grepl('/joined/',names(allQs16s))&grepl('_1\\.fastq',names(allQs16s))]
r2s16s<-allQs16s[!grepl('/joined/',names(allQs16s))&grepl('_2\\.fastq',names(allQs16s))]
newR1s<-allQs[grepl('_R1_',names(allQs))&grepl('Weimin_Plants_3_6_17',names(allQs))]
newR2s<-allQs[grepl('_R2_',names(allQs))&grepl('Weimin_Plants_3_6_17',names(allQs))]

pdf('out/quals.pdf')
  plot(1,1,type='n',ylim=c(0,40),xlim=c(1,length(allQs[[1]])),main='All',xlab='Position',ylab='Average quality')
  lines(apply(do.call(rbind,r1s),2,mean),col='red',lwd=2)
  lines(apply(do.call(rbind,r2s),2,mean),col='orange',lwd=2)
  lines(apply(do.call(rbind,r1s16s),2,mean),col='purple',lwd=2)
  lines(apply(do.call(rbind,r2s16s),2,mean),col='blue',lwd=2)
  lines(apply(do.call(rbind,newR1s),2,mean),col='green',lwd=2)
  lines(apply(do.call(rbind,newR2s),2,mean),col='cyan',lwd=2)
  abline(h=30,lty=2)
  legend('bottomleft',c('Chloro 1','Chloro 2','16S 1','16S 2','New chloro 1','New chloro 2'),col=c('red','orange','purple','blue','green','cyan'),lty=1,inset=.01,lwd=2)
  plot(1,1,type='n',ylim=c(0,40),xlim=c(1,length(allQs[[1]])),main='Chloro Read 1',xlab='Position',ylab='Average quality per sample')
  lapply(r1s,function(xx)lines(1:length(xx),xx,col='#00000044'))
  plot(1,1,type='n',ylim=c(0,40),xlim=c(1,length(allQs[[1]])),main='Chloro Read 2',xlab='Position',ylab='Average quality per sample')
  lapply(r2s,function(xx)lines(1:length(xx),xx,col='#00000044'))
  plot(1,1,type='n',ylim=c(0,40),xlim=c(1,length(allQs[[1]])),main='16S Read 1',xlab='Position',ylab='Average quality per sample')
  lapply(r1s16s,function(xx)lines(1:length(xx),xx,col='#00000044'))
  plot(1,1,type='n',ylim=c(0,40),xlim=c(1,length(allQs[[1]])),main='16S Read 2',xlab='Position',ylab='Average quality per sample')
  lapply(r2s16s,function(xx)lines(1:length(xx),xx,col='#00000044'))
  plot(1,1,type='n',ylim=c(0,40),xlim=c(1,length(allQs[[1]])),main='New chloro Read 1',xlab='Position',ylab='Average quality per sample')
  lapply(newR1s,function(xx)lines(1:length(xx),xx,col='#00000044'))
  plot(1,1,type='n',ylim=c(0,40),xlim=c(1,length(allQs[[1]])),main='New chloro Read 2',xlab='Position',ylab='Average quality per sample')
  lapply(newR2s,function(xx)lines(1:length(xx),xx,col='#00000044'))
dev.off()

r1<-read.fastq('data/Weimin_10_13_16_2/PA1044matK-bd_S51_L001_R1_001.fastq.gz')
r2<-read.fastq('data/Weimin_10_13_16_2/PA1044matK-bd_S51_L001_R2_001.fastq.gz')
seqs1<-seqSplit(r1$seq)
seqs2<-seqSplit(r2$seq)
quals1<-do.call(rbind,qualToInts(r1$qual))
quals2<-do.call(rbind,qualToInts(r2$qual))
pdf('test.pdf')
plot(apply(quals1,2,mean),type='l',col='red')
lines(apply(quals2,2,mean),col='blue')
dev.off()

align1<-levenAlign(r1$seq,mostAbundant(r1$seq),nThreads=40)
align2<-levenAlign(r2$seq,mostAbundant(r2$seq),nThreads=40)
png('r1_align.png',width=6000,height=6000,res=500)
plotDNA(removeGapCols(align1[[2]][order(leven(r1$seq,mostAbundant(r1$seq),nThreads=40),r1$seq)]))
dev.off()
png('r2_align.png',width=6000,height=6000,res=500)
plotDNA(removeGapCols(align2[[2]][order(leven(r2$seq,mostAbundant(r2$seq),nThreads=40),r2$seq)]))
dev.off()
png('r1.png',width=6000,height=6000,res=500)
plotDNA(sort(r1$seq))
dev.off()
png('r2.png',width=6000,height=6000,res=500)
plotDNA(sort(r2$seq))
dev.off()
png('r1_rev.png',width=6000,height=6000,res=500)
plotDNA(sort(reverseString(r1$seq)))
dev.off()
png('r2_rev.png',width=6000,height=6000,res=500)
plotDNA(sort(reverseString(r2$seq)))
dev.off()
