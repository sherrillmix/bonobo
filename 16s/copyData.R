source('../readSamples.R',chdir=TRUE)
#TODO
samples<-readLines('primates.txt')
targetDir<-'~/projects/animalPoop/bigRun/split/'
for(ii in samples){
  targets<-sort(list.files(targetDir,sprintf('^%s',ii),full.name=TRUE))
  if(length(targets)!=2)stop("Incorrect files in ",ii)
  links<-file.path('data',basename(targets))
  needLinked<-!file.exists(links)
  if(any(needLinked))file.symlink(targets[needLinked],links[needLinked])
  out<-sprintf('data/joined/%s.fastq',ii)

  cmd<-sprintf('~/installs/bbmap/bbmerge.sh in1=%s in2=%s out=%s t=10 2>%s',links[1],links[2],out,sub('fastq$','out',out))
  message(cmd)
  exit<-system(cmd)
  if(exit!=0)stop("Problem running pairing in",ii)
  system(sprintf('gzip %s',out))
}


