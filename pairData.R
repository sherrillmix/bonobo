
library(parallel)

fastqs<-list.files('data','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined',fastqs)]
bases<-file.path(basename(dirname(fastqs)),sub('_L[0-9]+_R[12]_.*$','',basename(fastqs)))

if(!dir.exists('work/pair'))dir.create('work/pair')
mclapply(bases,function(base){
  thisFiles<-sort(fastqs[bases==base])
  if(!dir.exists(file.path('work/pair',dirname(base))))dir.create(file.path('work/pair',dirname(base)))
  outFiles<-file.path('work/pair',c(sprintf(c('%s.match.fastq','%s.unmatch1.fastq','%s.unmatch2.fastq','%s.hist'),base)))
  if(length(thisFiles)!=2)stop('Number of files not equal to 2')
  cmd<-sprintf('~/installs/bbmap/bbmerge.sh in1=%s in2=%s out=%s outu1=%s outu2=%s t=3 ihist=%s',thisFiles[1],thisFiles[2],outFiles[1],outFiles[2],outFiles[3],outFiles[4])
  message(cmd)
  system(cmd)
  sapply(outFiles,function(xx)system(sprintf('gzip %s',xx)))
},mc.cores=3)
