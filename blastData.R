library(parallel)

fastqs<-list.files('data/','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined',fastqs)]

if(!file.exists('work/matk.nin'))system('makeblastdb -in ../chlorophyll/work/trimMatk.fa -dbtype nucl -out work/matk')
if(!file.exists('work/rbcl.nin'))system('makeblastdb -in ../chlorophyll/work/trimRbcl.fa -dbtype nucl -out work/rbcl')


if(!dir.exists('work/blast'))dir.create('work/blast')
check<-mclapply(fastqs,function(fastq){
  outFile<-sprintf('work/blast/%s',sub('\\.fastq.gz$','.blast.gz',basename(fastq)))
  if(grepl('rbcL',fastq))gene<-'rbcl'
  else if(grepl('matK',fastq))gene<-'matk'
  else stop('Unknown gene')
  cmd<-sprintf("zcat %s|awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,\">\");print}; if(P==4)P=0; P++}'|blastn -db work/%s -num_threads 5 -culling_limit 10 -outfmt 6|gzip > %s",fastq,gene,outFile)
  if(!file.exists(outFile)){
    message(cmd)
    system(cmd)
  }else{
    message("ALREADY DONE: ",cmd)
  }
},mc.cores=10)

