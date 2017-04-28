if(!exists('swarmData'))source("loadData.R")
library(phyloseq)
library(ape)
library(vipor)
library(Rtsne)
source('16s/myBiplot.R')

primerSet<-sub('^([^0-9]+)[0-9].*','\\1',names(swarmData))

pdf('out/pcoa.pdf',width=10,height=8)
for(ii in unique(primerSet)){
  thisPrimers<-names(swarmData)[primerSet==ii]
  isEnough<-apply(do.call(cbind,lapply(thisPrimers,function(primer)swarmData[[primer]][['isEnough']])),1,all)
  plotProp<-do.call(cbind,lapply(thisPrimers,function(primer)swarmData[[primer]][['props']][isEnough&rownames(swarmData[[primer]][['props']]) %in% rownames(samples),]))
  phyOtuW<-otu_table(plotProp,taxa_are_rows=FALSE)
  qiimeDataW<-phyloseq(otu_table=phyOtuW)
  brayDist<-distance(qiimeDataW,'bray')
  brayPca<-pcoa(brayDist)
  tsneBray<-Rtsne(brayDist,is_distance=TRUE,verbose=TRUE,perplexity=15,max_iter=3000)

  selectSamples<-samples[rownames(plotProp),]
  colorBrew<-c('#e41a1cBB','#377eb8BB','#4daf4aBB','#984ea3BB','#ff7f00BB','#ffff33BB','#a65628BB','#f781bfBB','#999999BB','#88ddffBB')
  nArea<-length(unique(selectSamples$area2))
  if(nArea>length(colorBrew))stop('Need to adjust colors for more areas')
  areaCols<-colorBrew[1:nArea]
  names(areaCols)<-unique(selectSamples$area2)
  areaPch<-sapply(names(areaCols),function(x)mostAbundant(selectSamples$chimpBonobo[selectSamples$area2==x]))
  malariaCols<-c('#00000022','#000000CC')
  malariaCols2<-rainbow.lab(2,alpha=.9,lightMultiple=.7)
  names(malariaCols2)<-names(malariaCols)<-c('Laverania negative','Laverania positive')
  speciesPch<-20+1:length(unique(selectSamples$Species))
  speciesCols<-rainbow.lab(length(unique(selectSamples$Species)),start=-2,end=2,alpha=.9,lightMultiple=.8)
  names(speciesCols)<-names(speciesPch)<-unique(selectSamples$chimpBonobo)

  predictors<-model.matrix(~0+Species+malaria+SIV+area,selectSamples)
  pos<-my.biplot.pcoa(brayPca,predictors,plot.axes=1:2,pch=speciesPch[selectSamples$chimpBonobo],bg=areaCols[selectSamples$area2],col=malariaCols[selectSamples$malaria+1],cex=2.2,lwd=2.5,mar=c(4,4,1.5,10),arrowsFilter=Inf)
  points(pos[selectSamples$SIV=='Pos',],col='#FF000099',cex=2.7,lwd=2)
  legend(
    par('usr')[2]+.01*diff(par('usr')[1:2]), 
    mean(par('usr')[3:4]),
    c(names(malariaCols),names(areaCols),names(speciesPch),'SIVPos'),
    col=c(malariaCols,rep(malariaCols,c(length(areaCols),length(speciesPch))),'#FF000099'),
    pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch,1),
    pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch)),NA),
    inset=.01,pt.lwd=3,pt.cex=2.5,
    xjust=0,xpd=NA
  )
  title(main=sprintf('%s PCoA PC %d and %d',ii,1,2))

  par(mar=c(4,4,1.5,10))
  plot(tsneBray$Y,pch=speciesPch[selectSamples$chimpBonobo],bg=areaCols[selectSamples$area2],col=malariaCols[selectSamples$malaria+1],cex=2.2,lwd=2.5,ylab='t-SNE 2',xlab='t-SNE 1',main=sprintf('%s t-SNE',ii))
  legend(
    par('usr')[2]+.01*diff(par('usr')[1:2]), 
    mean(par('usr')[3:4]),
    c(names(malariaCols),names(areaCols),names(speciesPch)),
    col=c(malariaCols,rep(malariaCols,c(length(areaCols),length(speciesPch)))),
    pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch),
    pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch))),
    inset=.01,pt.lwd=3,pt.cex=2.5,
    xjust=0,xpd=NA
  )
}
dev.off()
