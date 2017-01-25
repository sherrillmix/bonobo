source('../readSamples.R',chdir=TRUE)
samples<-readLines('primates.txt')
targetDir<-'~/projects/animalPoop/bigRun/split/'
for(ii in samples){
  targets<-list.files(targetDir,sprintf('^%s',ii),full.name=TRUE)
  if(length(targets)!=2)stop("Incorrect files in ",ii)
  file.symlink(targets,file.path('data',basename(targets)))
}
