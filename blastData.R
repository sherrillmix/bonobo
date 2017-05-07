library(parallel)

fastqs<-list.files('data/','_R[12]_.*\\.fastq\\.gz$',recursive=TRUE,full.names=TRUE)
fastqs<-fastqs[!grepl('Undetermined',fastqs)]

#if(!file.exists('work/matk.nin'))system('makeblastdb -in ../chlorophyll/work/trimMatk.fa -dbtype nucl -out work/matk')
#if(!file.exists('work/rbcl.nin'))system('makeblastdb -in ../chlorophyll/work/trimRbcl.fa -dbtype nucl -out work/rbcl')

if(!file.exists('work/matk.nin')){system('zcat ../chlorophyll/matk/matk.fa.gz | makeblastdb -in - -title matk -dbtype nucl -out work/matk')}
if(!file.exists('work/rbcl.nin')){system('zcat ../chlorophyll/rbcl/rbcl.fa.gz | makeblastdb -in - -title rbcl -dbtype nucl -out work/rbcl')}

if(!dir.exists('work/blast'))dir.create('work/blast')
check<-mclapply(fastqs,function(fastq){
  outFile<-sprintf('work/blast/%s',sub('\\.fastq.gz$','.blast.gz',basename(fastq)))
  if(grepl('rbcL',fastq))gene<-'rbcl'
  else if(grepl('matK',fastq))gene<-'matk'
  else stop('Unknown gene')
  cmd<-sprintf("zcat %s|awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,\">\");print}; if(P==4)P=0; P++}'|blastn -db work/%s -num_threads 4 -culling_limit 10 -outfmt 6|gzip > %s",fastq,gene,outFile)
  if(!file.exists(outFile)){
    message(cmd)
    system(cmd)
  }else{
    message("ALREADY DONE: ",cmd)
  }
},mc.cores=14)

contigFa<-list.files('work/swarm','.fa.gz$',full.names=TRUE)
for(fasta in contigFa){
  outFile<-sub('fa.gz$','blast.gz',fasta)
  if(grepl('rbcL',fasta))gene<-'rbcl'
  else if(grepl('matK',fasta))gene<-'matk'
  else stop('Unknown gene')
  cmd<-sprintf("zcat %s|blastn -db work/%s -num_threads 30 -culling_limit 10 -outfmt 6|gzip > %s",fasta,gene,outFile)
  if(!file.exists(outFile)){
    message(cmd)
    system(cmd)
  }else{
    message("ALREADY DONE: ",cmd)
  }
}

contigFa<-list.files('work/swarmPair','.fa.gz$',full.names=TRUE)
for(fasta in contigFa){
  outFile<-sub('fa.gz$','blast.gz',fasta)
  if(grepl('rbcL',fasta))gene<-'rbcl'
  else if(grepl('matK',fasta))gene<-'matk'
  else stop('Unknown gene')
  cmd<-sprintf("zcat %s|parallel --block 1m --recstart '>' -L2 -j 12 --pipe blastn -db work/%s -num_threads 5 -culling_limit 30 -outfmt 6 -query -|gzip > %s",fasta,gene,outFile)
  if(!file.exists(outFile)){
    message(cmd)
    system(cmd)
  }else{
    message("ALREADY DONE: ",cmd)
  }
}
