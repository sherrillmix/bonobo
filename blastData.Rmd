library(parallel)

#TODO
if(!file.exists('work/matk.nin')){system('zcat ../chlorophyll/matk/matk.fa.gz | makeblastdb -in - -title matk -dbtype nucl -out work/matk')}
if(!file.exists('work/rbcl.nin')){system('zcat ../chlorophyll/rbcl/rbcl.fa.gz | makeblastdb -in - -title rbcl -dbtype nucl -out work/rbcl')}

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
