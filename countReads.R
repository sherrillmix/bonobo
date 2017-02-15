library(dnar)
fastqs<-list.files('data','_R[12]_.*fastq.gz$',recursive=TRUE,full.names=TRUE)
nReads<-unlist(cleanMclapply(fastqs,function(xx){library(dnar);nrow(read.fastq(xx))},mc.cores=30))
names(nReads)<-sub('_S[0-9]+_L[0-9]+_R([0-9]).*$','',basename(fastqs))

