if(!exists('otuTab'))source('runQiime.R')

library(phyloseq)
library(ape)
tree<-multi2di(read_tree('work/qiime/rep_set.tre'))
phyOtu<-otu_table(otuTab[rownames(otuTab) %in% tree$tip.label,samples$name],taxa_are_rows=TRUE)
qiimeData<-phyloseq(otu_table=phyOtu,phy_tree=tree)
#make sure tree is bifurcating or breaks UniFrac without error
uniDist<-UniFrac(qiimeData,weighted=TRUE)
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


#just bonobos
phyOtu2<-otu_table(otuTab[rownames(otuTab) %in% tree$tip.label,samples$name[samples$bonobo]],taxa_are_rows=TRUE)
qiimeData2<-phyloseq(otu_table=phyOtu2,phy_tree=tree)
#make sure tree is bifurcating or breaks UniFrac without error
uniDist2<-UniFrac(qiimeData2,weighted=TRUE)
library(ape)
uniPca2<-pcoa(uniDist2)
predictors2<-model.matrix(~malaria+area+season,samples[samples$bonobo,])
colnames(predictors2)[colnames(predictors2)=='malariaTRUE']<-'malariaPos'

splitAreas<-withAs(samples=samples[samples$bonobo,],ifelse(samples$area %in% c('TL-E','TL-W'),samples$area,samples$Species))
splitAreaCols<-rainbow.lab(length(unique(splitAreas)))
names(splitAreaCols)<-unique(splitAreas)
pdf('out/bonobo_pcoa.pdf',height=9,width=9)
  sapply(list(1:2,3:4,5:6),function(axes){
    pos<-my.biplot.pcoa(uniPca2,predictors2,plot.axes=axes,pch=21,bg=splitAreaCols[splitAreas],col=malariaCols[samples[samples$bonobo,'malaria']+1],cex=2.25,lwd=4)
    legend('bottomright',names(splitAreaCols),col=malariaCols[1],pch=21,pt.bg=splitAreaCols,inset=.01,pt.lwd=4,pt.cex=2.5)
    title(main=sprintf('TL and KR PC %d and %d',axes[1],axes[2]))
  })
dev.off()


sampleNs<-apply(otuProp[,samples$name],2,sum)
#apply(otuTab[,samples$name],1,function(x){
  #glm(cbind(x,sampleNs)
#})
