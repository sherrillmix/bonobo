library(parallel)
library(dnar)

source('../readSamples.R',chdir=TRUE)

runQiime<-function(seqs,storeDir=NULL){
  if(any(is.na(seqs)))stop(simpleError('NAs in seqs'))
  seqIds<-1:length(seqs)
  readDir<-tempfile()
  dir.create(readDir)
  readFile<-file.path(readDir,'XXX.fa')
  outDir<-tempfile()
  seqNames<-sprintf('XXX_%d',seqIds)
  write.fa(seqNames,seqs,readFile)
  #miniconda doesn't like sh so need to use bash
  cmd<-sprintf('echo "source activate qiime1; pick_de_novo_otus.py --input %s --output %s --parallel --jobs_to_start 32 --force"|bash',readFile,outDir)
  message(cmd)
  exit<-system(cmd)
  if(exit!=0)stop(simpleError('Problem running qiime'))
  #get otu assignments
  assigns<-strsplit(readLines(file.path(outDir,'uclust_picked_otus/XXX_otus.txt')),'\t')
  names(assigns)<-sapply(assigns,'[[',1)
  assigns<-lapply(assigns,'[',-1)
  otus<-rep(names(assigns),sapply(assigns,length))
  names(otus)<-unlist(assigns)
  out<-otus[seqNames]
  #get taxa assignments
  taxa<-strsplit(readLines(file.path(outDir,'uclust_assigned_taxonomy/XXX_rep_set_tax_assignments.txt')),'\t')
  names(taxa)<-sapply(taxa,'[[',1)
  taxa<-sapply(taxa,'[[',2)
  #get sequences
  seqs<-read.fa(file.path(outDir,'pynast_aligned_seqs/XXX_rep_set_aligned_pfiltered.fasta'))
  seqs<-structure(seqs$seq,.Names=seqs$name)
  if(!is.null(storeDir)){
    #avoid problems with cross filesystem copies
    tmpDir<-tempfile(tmpdir='.')
    dir.create(tmpDir)
    file.copy(outDir,tmpDir)
    file.rename(file.path(tmpDir,basename(outDir)),storeDir)
  }
  return(list('otus'=out,'seqs'=seqs,'taxa'=taxa))
}

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

samples$name<-sapply(samples$Code,function(x)colnames(otuTab)[grep(sprintf('_%s$',x),colnames(otuTab))])
rownames(samples)<-samples$name
#make sure same order as otu table
if(any(!samples$name %in% colnames(otuTab)))stop('Sample missing from OTUs')
samples<-samples[colnames(otuTab)[colnames(otuTab) %in% samples$name],]

library(phyloseq)
tree<-multi2di(read_tree('work/qiime/rep_set.tre'))
phyOtu<-otu_table(otuTab[rownames(otuTab) %in% tree$tip.label,samples$name],taxa_are_rows=TRUE)
qiimeData<-phyloseq(otu_table=phyOtu,phy_tree=tree)
#make sure tree is bifurcating or breaks UniFrac without error
uniDist<-UniFrac(qiimeData,weighted=TRUE)
library(ape)
uniPca<-pcoa(uniDist)
predictors<-model.matrix(~0+Species+malaria+SIV+area,samples)
colnames(predictors)<-sub('^Species','',colnames(predictors))
colnames(predictors)[colnames(predictors)=='malariaTRUE']<-'malariaPos'


#areaCols<-rainbow.lab(length(unique(samples$area)),alpha=.8)
source('myBiplot.R')
colorBrew<-c('#e41a1cBB','#377eb8BB','#4daf4aBB','#984ea3BB','#ff7f00BB','#ffff33BB','#a65628BB','#f781bfBB','#999999BB')
nArea<-length(unique(samples$area))
if(nArea>length(colorBrew))stop('Need to adjust colors for more areas')
areaCols<-colorBrew[1:nArea]
names(areaCols)<-unique(samples$area)
areaPch<-sapply(names(areaCols),function(x)mostAbundant(samples$Species[samples$area==x]))
malariaCols<-c('#00000022','#000000CC')
malariaCols2<-rainbow.lab(2,alpha=.9,lightMultiple=.7)
names(malariaCols2)<-names(malariaCols)<-c('PlasmoNeg','PlasmoPos')
speciesPch<-20+1:length(unique(samples$Species))
speciesCols<-rainbow.lab(length(unique(samples$Species)),start=-2,end=2,alpha=.9,lightMultiple=.8)
names(speciesCols)<-names(speciesPch)<-unique(samples$Species)
#split out TL2-E and -W, bonobos and KR and chimps
splitAreas<-ifelse(samples$area %in% c('TL-E','TL-W','KR'),samples$area,samples$Species)
splitAreaCols<-rainbow.lab(length(unique(splitAreas)))
names(splitAreaCols)<-unique(splitAreas)[order(unique(splitAreas) %in% c('P.t. schweinfurthii'))]
pdf('out/pcoa.pdf',height=9,width=9)
  sapply(list(1:2,3:4,5:6),function(axes){
    pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=speciesPch[samples$Species],bg=areaCols[samples$area],col=malariaCols[samples$malaria+1],cex=2.2,lwd=2.5)
    points(pos[samples$SIV=='Pos',],col='#FF000099',cex=2.7,lwd=2)
    legend(
      'bottomright',
      c(names(malariaCols),names(areaCols),names(speciesPch),'SIVPos'),
      col=c(malariaCols,rep(malariaCols,c(length(areaCols),length(speciesPch))),'#FF000099'),
      pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch,1),
      pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch)),NA),
      inset=.01,pt.lwd=3,pt.cex=2.5)
    title(main=sprintf('All variables PC %d and %d',axes[1],axes[2]))
    #species
    pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=21,bg=speciesCols[samples$Species],col="#00000077",cex=1.8,lwd=2.5)
    legend('bottomright',names(speciesCols),col='#00000077',pch=21,pt.bg=speciesCols,inset=.01,pt.lwd=3,pt.cex=2)
    title(main=sprintf('Species PC %d and %d',axes[1],axes[2]))
    #malaria
    pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=21,bg=malariaCols2[samples$malaria+1],col="#00000077",cex=1.8,lwd=2.5)
    legend('bottomright',names(malariaCols),col='#00000077',pch=21,pt.bg=malariaCols2,inset=.01,pt.lwd=3,pt.cex=2)
    title(main=sprintf('Malaria PC %d and %d',axes[1],axes[2]))
    pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=21,bg=speciesCols[samples$Species],col=malariaCols[samples$malaria+1],cex=2.25,lwd=4)
    legend('bottomright',as.vector(outer(names(speciesCols),names(malariaCols),paste,sep=' ')),col=as.vector(malariaCols[outer(names(speciesCols),names(malariaCols),function(x,y)y)]),pch=21,pt.bg=as.vector(speciesCols[outer(names(speciesCols),names(malariaCols),function(x,y)x)]),inset=.01,pt.lwd=4,pt.cex=2.5)
    title(main=sprintf('Species/malaria PC %d and %d',axes[1],axes[2]))
    pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=21,bg=splitAreaCols[splitAreas],col=malariaCols[1],cex=2.25,lwd=4)
    legend('bottomright',names(splitAreaCols),col=malariaCols[1],pch=21,pt.bg=splitAreaCols,inset=.01,pt.lwd=4,pt.cex=2.5)
    title(main=sprintf('TL and KR PC %d and %d',axes[1],axes[2]))
    #show sample name text
    #show sample name text
    #bak<-par('cex')
    #par('cex'=.65)
    #biplot.pcoa(uniPca,predictors,plot.axes=axes)
    #title(main='Sample names')
    #par('cex'=bak)
  })
dev.off()
