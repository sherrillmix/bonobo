if(!exists('swarmData'))source("loadData.R")
library(phyloseq)
library(ape)
library(vipor)
library(Rtsne)
source('16s/myBiplot.R')

for(ii in names(swarmData)){
  plotProp<-swarmData[[ii]][['props']][swarmData[[ii]][['isEnough']]&rownames(swarmData[[ii]][['props']]) %in% rownames(samples),]
  phyOtuW<-otu_table(plotProp,taxa_are_rows=FALSE)
  qiimeDataW<-phyloseq(otu_table=phyOtuW)
  brayDist<-distance(qiimeDataW,'bray')
  brayPca<-pcoa(brayDist)
  tsneBray<-Rtsne(brayDist,is_distance=TRUE,verbose=TRUE,perplexity=15,max_iter=3000)

  importance<-brayPca$values$Relative_eig
  colnames(brayPca$vectors)<-sprintf('Principal coordinate %d (%d%% of variance)',1:length(importance),round(importance*100))[1:ncol(brayPca$vectors)]
  selectSamples<-samples[rownames(plotProp),]
  selectSamples<-selectSamples[order(selectSamples$bonobo,selectSamples$area2,selectSamples$malaria),]
  colorBrew<-c('#e41a1cBB','#377eb8BB','#4daf4aBB','#984ea3BB','#ff7f00BB','#ffff33BB','#a65628BB','#f781bfBB','#999999BB','#88ddffBB')
  nArea<-length(unique(selectSamples$area2))
  if(nArea>length(colorBrew))stop('Need to adjust colors for more areas')
  areaCols<-colorBrew[1:nArea]
  names(areaCols)<-unique(selectSamples$area2)
  areaPch<-sapply(names(areaCols),function(x)mostAbundant(selectSamples$chimpBonobo[selectSamples$area2==x]))
  #malariaCols3<-c(NA,'#000000CC')
  malariaCols3<-c('#00000022','#000000CC')
  malariaCols<-c('#00000022','#000000CC')
  malariaCols2<-rainbow.lab(2,alpha=.9,lightMultiple=.7)
  names(malariaCols3)<-names(malariaCols2)<-names(malariaCols)<-c('Laverania negative','Laverania positive')
  speciesPch<-20+1:length(unique(selectSamples$Species))
  #speciesCols<-c('#FF0000CC','#0000FFCC')
  speciesCols<-rainbow.lab(length(unique(selectSamples$Species)),start=-2,end=2,alpha=.9,lightMultiple=.8)
  names(speciesCols)<-names(speciesPch)<-sort(unique(selectSamples$chimpBonobo))

  predictors<-model.matrix(~0+Species+malaria+SIV+area,selectSamples)
  #pos<-my.biplot.pcoa(brayPca,predictors,plot.axes=1:2,pch=speciesPch[selectSamples$chimpBonobo],bg=areaCols[selectSamples$area2],col=malariaCols[selectSamples$malaria+1],cex=2.2,lwd=2.5,mar=c(4,4,1.5,10),arrowsFilter=Inf)
  #legend(
    #par('usr')[2]+.01*diff(par('usr')[1:2]), 
    #mean(par('usr')[3:4]),
    #c(names(malariaCols),names(areaCols),names(speciesPch)),
    #col=c(malariaCols,rep(malariaCols,c(length(areaCols),length(speciesPch)))),
    #pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch),
    #pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch))),
    #inset=.01,pt.lwd=3,pt.cex=2.5,
    #xjust=0,xpd=NA
  #)

  pdf(sprintf('out/pcoa_%s.pdf',ii),width=7,height=7)
    pos<-my.biplot.pcoa(brayPca,predictors,plot.axes=1:2,pch=21,bg=speciesCols[selectSamples$chimpBonobo],col=malariaCols3[selectSamples$malaria+1],cex=2.25,lwd=4,arrowsFilter=Inf,las=1,mgp=c(2.75,.75,0))
    legend('bottomright',as.vector(outer(names(speciesCols),names(malariaCols3),paste,sep=' ')),col=as.vector(malariaCols3[outer(names(speciesCols),names(malariaCols3),function(x,y)y)]),pch=21,pt.bg=as.vector(speciesCols[outer(names(speciesCols),names(malariaCols3),function(x,y)x)]),inset=.01,pt.lwd=4,pt.cex=2.5,bty='n')
    title(main=sprintf('%s',ii,1,2))
  dev.off()

  pdf(sprintf('out/tsne_%s.pdf',ii),width=9,height=7)
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
  dev.off()

}
