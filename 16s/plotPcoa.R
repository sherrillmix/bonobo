if(!exists('otuTab'))source('runQiime.R')

library(phyloseq)
library(ape)
selectSamples<-samples[samples$isEnough,]
tree<-multi2di(read_tree('work/qiime/rep_set.tre'))
phyOtu<-otu_table(otuTab[rownames(otuTab) %in% tree$tip.label,selectSamples$name],taxa_are_rows=TRUE)
qiimeData<-phyloseq(otu_table=phyOtu,phy_tree=tree)
#make sure tree is bifurcating or breaks UniFrac without error
uniDist<-UniFrac(qiimeData,weighted=FALSE)
uniDistW<-UniFrac(qiimeData,weighted=TRUE)
uniPca<-pcoa(uniDist)
predictors<-model.matrix(~0+Species+malaria+SIV+area,selectSamples)
colnames(predictors)<-sub('^Species','',colnames(predictors))
colnames(predictors)[colnames(predictors)=='malariaTRUE']<-'malariaPos'


#areaCols<-rainbow.lab(length(unique(samples$area)),alpha=.8)
source('myBiplot.R')
#last color made up
colorBrew<-c('#e41a1cBB','#377eb8BB','#4daf4aBB','#984ea3BB','#ff7f00BB','#ffff33BB','#a65628BB','#f781bfBB','#999999BB','#88ddffBB')
nArea<-length(unique(selectSamples$area2))
if(nArea>length(colorBrew))stop('Need to adjust colors for more areas')
areaCols<-colorBrew[1:nArea]
names(areaCols)<-unique(selectSamples$area2)
areaPch<-sapply(names(areaCols),function(x)mostAbundant(selectSamples$Species[selectSamples$area2==x]))
malariaCols<-c('#00000022','#000000CC')
malariaCols2<-rainbow.lab(2,alpha=.9,lightMultiple=.7)
names(malariaCols2)<-names(malariaCols)<-c('PlasmoNeg','PlasmoPos')
speciesPch<-20+1:length(unique(selectSamples$Species))
speciesCols<-rainbow.lab(length(unique(selectSamples$Species)),start=-2,end=2,alpha=.9,lightMultiple=.8)
names(speciesCols)<-names(speciesPch)<-unique(selectSamples$Species)
#split out TL2-E and -W, bonobos and KR and chimps
splitAreas<-ifelse(selectSamples$area2 %in% c('TL-E','TL-W','KR','TL-NE'),selectSamples$area2,selectSamples$Species)
splitAreaCols<-rainbow.lab(length(unique(splitAreas)))
names(splitAreaCols)<-unique(splitAreas)[order(unique(splitAreas) %in% c('P.t. schweinfurthii'))]
pdf('out/pcoa.pdf',height=9,width=9)
  sapply(list(1:2,3:4,5:6),function(axes){
    pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=speciesPch[selectSamples$Species],bg=areaCols[selectSamples$area2],col=malariaCols[selectSamples$malaria+1],cex=2.2,lwd=2.5)
    points(pos[selectSamples$SIV=='Pos',],col='#FF000099',cex=2.7,lwd=2)
    legend(
      'bottomright',
      c(names(malariaCols),names(areaCols),names(speciesPch),'SIVPos'),
      col=c(malariaCols,rep(malariaCols,c(length(areaCols),length(speciesPch))),'#FF000099'),
      pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch,1),
      pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch)),NA),
      inset=.01,pt.lwd=3,pt.cex=2.5)
    title(main=sprintf('All variables PC %d and %d',axes[1],axes[2]))
    text(pos,selectSamples$Code,cex=.25)
    #species
    pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=21,bg=speciesCols[selectSamples$Species],col="#00000077",cex=1.8,lwd=2.5)
    legend('bottomright',names(speciesCols),col='#00000077',pch=21,pt.bg=speciesCols,inset=.01,pt.lwd=3,pt.cex=2)
    title(main=sprintf('Species PC %d and %d',axes[1],axes[2]))
    #malaria
    pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=21,bg=malariaCols2[selectSamples$malaria+1],col="#00000077",cex=1.8,lwd=2.5)
    legend('bottomright',names(malariaCols),col='#00000077',pch=21,pt.bg=malariaCols2,inset=.01,pt.lwd=3,pt.cex=2)
    title(main=sprintf('Malaria PC %d and %d',axes[1],axes[2]))
    pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=21,bg=speciesCols[selectSamples$Species],col=malariaCols[selectSamples$malaria+1],cex=2.25,lwd=4)
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


#just bonobos
selectSamples2<-samples[samples$bonobo&samples$isEnough,]
phyOtu2<-otu_table(otuTab[rownames(otuTab) %in% tree$tip.label,selectSamples2$name],taxa_are_rows=TRUE)
qiimeData2<-phyloseq(otu_table=phyOtu2,phy_tree=tree)
#make sure tree is bifurcating or breaks UniFrac without error
uniDist2<-UniFrac(qiimeData2,weighted=FALSE)
uniDistW2<-UniFrac(qiimeData2,weighted=TRUE)
selectSamples3<-samples[!samples$bonobo&samples$isEnough,]
phyOtu3<-otu_table(otuTab[rownames(otuTab) %in% tree$tip.label,selectSamples3$name],taxa_are_rows=TRUE)
qiimeData3<-phyloseq(otu_table=phyOtu3,phy_tree=tree)
uniDist3<-UniFrac(qiimeData3,weighted=FALSE)
uniDistW3<-UniFrac(qiimeData3,weighted=TRUE)
library(ape)
uniPca2<-pcoa(uniDist2)
predictors2<-model.matrix(~malaria+area+season,selectSamples2[selectSamples2$bonobo,])
colnames(predictors2)[colnames(predictors2)=='malariaTRUE']<-'malariaPos'

splitAreas<-withAs(selectSamples2=selectSamples2[selectSamples2$bonobo,],ifelse(selectSamples2$area %in% c('TL-E','TL-W'),selectSamples2$area,selectSamples2$Species))
splitAreaCols<-rainbow.lab(length(unique(splitAreas)))
names(splitAreaCols)<-unique(splitAreas)
pdf('out/bonobo_pcoa.pdf',height=9,width=9)
  sapply(list(1:2,3:4,5:6),function(axes){
    pos<-my.biplot.pcoa(uniPca2,predictors2,plot.axes=axes,pch=21,bg=splitAreaCols[splitAreas],col=malariaCols[selectSamples2[selectSamples2$bonobo,'malaria']+1],cex=2.25,lwd=4)
    legend('bottomright',names(splitAreaCols),col=malariaCols[1],pch=21,pt.bg=splitAreaCols,inset=.01,pt.lwd=4,pt.cex=2.5)
    title(main=sprintf('TL and KR PC %d and %d',axes[1],axes[2]))
  })
dev.off()


sampleNs<-apply(otuTab[,samples$name],2,sum)
#apply(otuTab[,samples$name],1,function(x){
  #glm(cbind(x,sampleNs)
#})

splitDist<-function(dists,splits){
  uniqSplits<-unique(splits)
  names(uniqSplits)<-uniqSplits
  splitDists<-lapply(uniqSplits,function(xx){
    lapply(uniqSplits,function(yy){
      out<-dists[splits==xx,splits==yy]
      if(xx==yy)return(out[upper.tri(out)])
      else return(as.vector(out))
    })
  })
  out<-data.frame('x'=NA,'y'=NA,'dist'=NA,stringsAsFactors=FALSE)[0,]
  for(ii in uniqSplits){
    for(jj in uniqSplits){
      out<-rbind(out,data.frame('x'=ii,'y'=jj,'dist'=splitDists[[ii]][[jj]],stringsAsFactors=FALSE))
    }
  }
  return(out)
}

library(vipor)
speciesOtu<-otu_table(otuTab[rownames(otuTab) %in% tree$tip.label,],taxa_are_rows=TRUE)
speciesQiimeData<-phyloseq(otu_table=speciesOtu,phy_tree=tree)
#make sure tree is bifurcating or breaks UniFrac without error
speciesUniDist<-UniFrac(speciesQiimeData,weighted=TRUE)
uniDistMat<-as.matrix(speciesUniDist)
species<-sub('(Toddler)?_.*$','',rownames(uniDistMat))
splits<-splitDist(uniDistMat,species)
meanDists<-withAs(splits=splits[splits$x=="Bonobo"&splits$y!='Primate',],tapply(splits$dist,paste(splits$x,splits$y,sep='\n'),mean))
pdf('out/species.pdf',width=10)
  par(lheight=.7,mar=c(9,4,.1,.1))
  withAs(splits=splits[splits$x=="Bonobo"&splits$y!='Primate',],vpPlot(factor(paste(splits$x,splits$y,sep='\n'),levels=names(sort(meanDists))),splits$dist,las=2,ylab='Pairwise UniFrac distance',las=2))
  heatmap(tapply(splits$dist,list(splits$x,splits$y),mean),scale='none',col=rev(heat.colors(100)))
dev.off()

library(Rtsne)
tsne<-Rtsne(uniDist,is_distance=TRUE,verbose=TRUE,perplexity=10,max_iter=3000)
tsne2<-Rtsne(uniDist2,is_distance=TRUE,verbose=TRUE,perplexity=5,max_iter=3000)
tsne3<-Rtsne(uniDist3,is_distance=TRUE,verbose=TRUE,perplexity=5,max_iter=3000)
#tsne<-Rtsne(t(otuTab[rownames(otuTab) %in% tree$tip.label,selectSamples$name]),verbose=TRUE,perplexity=10)
tsneW<-Rtsne(uniDistW,is_distance=TRUE,verbose=TRUE,perplexity=10,max_iter=3000)
tsneW2<-Rtsne(uniDistW2,is_distance=TRUE,verbose=TRUE,perplexity=5,max_iter=3000)
tsneW3<-Rtsne(uniDistW3,is_distance=TRUE,verbose=TRUE,perplexity=5,max_iter=3000)
pdf('out/tsne.pdf',height=8,width=10)
  inputs<-list('Unweighted'=list(tsne,tsne2,tsne3),'Weighted'=list(tsneW,tsneW2,tsneW3))
  for(tsneType in names(inputs)){
    tsnes<-inputs[[tsneType]]
    par(mar=c(4,4,1.5,9))
    plot(tsnes[[1]]$Y,pch=speciesPch[selectSamples$Species],bg=areaCols[selectSamples$area2],col=malariaCols[selectSamples$malaria+1],cex=2.2,lwd=2.5,ylab='t-SNE 2',xlab='t-SNE 1',main=sprintf('Chimp and bonobo %s',tsneType))
    #points(tsnes[[1]]$Y[selectSamples$SIV=='Pos',],col='#FF000099',cex=2.7,lwd=2)
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
    text(tsnes[[1]]$Y,selectSamples$Code,cex=.25)
    #bonobo only
    plot(tsnes[[2]]$Y,pch=speciesPch[selectSamples2$Species],bg=areaCols[selectSamples2$area2],col=malariaCols[selectSamples2$malaria+1],cex=2.2,lwd=2.5,ylab='t-SNE 2',xlab='t-SNE 1',main='Bonobo')
    points(tsnes[[2]]$Y[selectSamples2$SIV=='Pos',],col='#FF000099',cex=2.7,lwd=2)
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
    text(tsnes[[2]]$Y,selectSamples2$Code,cex=.25)
    #bonobo only
    plot(tsnes[[3]]$Y,pch=speciesPch[selectSamples3$Species],bg=areaCols[selectSamples3$area2],col=malariaCols[selectSamples3$malaria+1],cex=2.2,lwd=2.5,ylab='t-SNE 2',xlab='t-SNE 1',main='Chimp')
    points(tsnes[[3]]$Y[selectSamples3$SIV=='Pos',],col='#FF000099',cex=2.7,lwd=2)
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
    text(tsnes[[3]]$Y,selectSamples3$Code,cex=.25)
  }
dev.off()


