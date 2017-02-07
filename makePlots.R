
if(!exists('chimpFilter')){
  source('parseBlast.R')
  chimpFilter<-lapply(taxas,function(x)x[x[,'class']!='Mammalia'|is.na(x[,'class']),])
  nChimpReads<-sapply(taxas,nrow)-sapply(chimpFilter,nrow)
  #dont need currently
  rm(taxas)
  rm(taxonomy)
}


#shortNames<-sub('_S[0-9]+_L[0-9]+_R([0-9]).*$','_R\\1',sub('blast_trim_','',names(chimpFilter)))
shortNames<-sub('_S[0-9]+_L[0-9]+_R([0-9]).*$','',sub('blast_trim_','',names(chimpFilter)))
taxaTab<-table(unlist(lapply(chimpFilter,function(x)x$best)),rep(shortNames,sapply(chimpFilter,nrow)))
propTab<-apply(taxaTab,2,function(x)x/sum(x))
cutTab<-propTab[apply(propTab,1,max)>.01,]
print(apply(taxaTab,2,sum))

breaks<-c(0,seq(min(cutTab)*.999,max(cutTab)*1.001,length.out=501))
cols<-c('white',tail(rev(heat.colors(520)),500))

insetLegend<-function(breaks,cols,insetPos=c(.015,.025,.25,.04)){
  insetPos<-c(grconvertX(insetPos[1],'nfc','user'),grconvertY(insetPos[2],'nfc','user'),grconvertX(insetPos[3],'nfc','user'),grconvertY(insetPos[4],'nfc','user'))
  breakPos<-((breaks[-1])-(min(breaks[-1])))/max((breaks[-1])-(min(breaks[-1])))*(insetPos[3]-insetPos[1])+insetPos[1]
  rect(breakPos[-1]+1e-3,insetPos[2],breakPos[-length(breakPos)],insetPos[4],col=cols[-1],xpd=NA,border=NA)
  rect(insetPos[1],insetPos[2],insetPos[3],insetPos[4],xpd=NA)
  prettyLabs<-pretty(breaks)
  prettyLabs<-prettyLabs[prettyLabs<max(breaks)]
  prettyPos<-prettyLabs
  prettyPos<-(prettyLabs-(min(breaks[-1])))/((max(breaks[-1]))-(min(breaks[-1])))*(insetPos[3]-insetPos[1])+insetPos[1]
  segments(prettyPos,insetPos[2],prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.1,xpd=NA)
  text(prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.175,prettyLabs,xpd=NA,adj=c(.5,1),cex=.85)
  text(mean(insetPos[c(1,3)]),insetPos[4]+diff(insetPos[c(2,4)])*.45,"Proportion of sample",xpd=NA,adj=c(.5,0))
}

pdf('out/heat.pdf',height=30,width=30)
  heatmap(cutTab,col=cols,breaks=breaks,margins=c(5,9),scale='none')
  insetLegend(breaks,cols,insetPos=c(.01,.93,.19,.945))
dev.off()

types<-sub('^([^0-9]*[0-9]+).*$','\\1',colnames(cutTab))
types2<-sub("_R[0-9]",'',sub('-bd','',sub('^[^0-9]*[0-9]+','',colnames(cutTab))))
newOrder<-order(types,types2)
cutTab<-cutTab[order(apply(cutTab,1,sum)),newOrder]
types<-types[newOrder]
types2<-types2[newOrder]
#cutTab<-cutTab[order(apply(cutTab,1,sum)),]
pdf('out/heatSort.pdf',width=40,height=30)
  par(mar=c(.1,20,7,.1))
  image(1:ncol(cutTab),1:nrow(cutTab),t(cutTab),col=cols,breaks=breaks,xlab='',ylab='',xaxt='n',yaxt='n')
  axis(3,1:ncol(cutTab),colnames(cutTab),las=2)
  axis(2,1:nrow(cutTab),rownames(cutTab),las=1)
  abline(v=which(types[-1]!=types[-length(types)])+.5)
  abline(h=2:nrow(cutTab)-.5,col='#00000033')
  abline(h=2:nrow(cutTab)-.5,col='#00000033')
  abline(v=2:ncol(cutTab)-.5,col='#00000011')
  box()
  insetLegend(breaks,cols,insetPos=c(.01,.975,.075,.99))
  for(ii in unique(types2)){
    thisCut<-cutTab[,types2==ii]
    image(1:ncol(thisCut),1:nrow(thisCut),t(thisCut),col=cols,breaks=breaks,xlab='',ylab='',xaxt='n',yaxt='n')
    axis(3,1:ncol(thisCut),colnames(thisCut),las=2)
    axis(2,1:nrow(thisCut),rownames(thisCut),las=1)
    abline(v=which(types[types2==ii][-1]!=types[types2==ii][-length(types[types2==ii])])+.5)
    abline(h=2:nrow(thisCut)-.5,col='#00000033')
    abline(h=2:nrow(thisCut)-.5,col='#00000033')
    abline(v=2:ncol(thisCut)-.5,col='#00000011')
  }
dev.off()


table(sub('^(.*[0-9]+).*','\\1',shortNames),sub('^.*[0-9]+','',shortNames))
length(chimpFilter)
taxaColumns<-c("superkingdom","phylum","class","order","family","genus","species")
condenseChimpFilter<-cacheOperation('work/condenseChimp.Rdat',cleanMclapply,chimpFilter,function(xx,taxaColumns)t(apply(xx[,taxaColumns],1,cumpaste,sep='_|_')),taxaColumns,mc.cores=50,extraCode='library(dnar)',EXCLUDE='chimpFilter')
gc()

print(unifracMatrix(condenseChimpFilter[[4]],condenseChimpFilter[[5]],checkUpstream=FALSE))
print(unifracMatrix(condenseChimpFilter[[5]],condenseChimpFilter[[4]],checkUpstream=FALSE))

#abusing cleanMclapply to get smaller enviro
unifracDists<-cleanMclapply(1,function(xx,condenseChimpFilter)cacheOperation('work/unifrac.Rdat',unifracMatrix,condenseChimpFilter,vocal=TRUE,weighted=TRUE,checkUpstream=FALSE,mc.cores=1,EXCLUDE='condenseChimpFilter'),condenseChimpFilter,extraCode='library(dnar)')[[1]]

source('readSamples.R')
sampleNames<-sub('(matK|rbcL).*','',names(condenseChimpFilter))
primers<-sub('.*(matK|rbcL).*_R([0-9]+)_.*','\\1\\2',names(condenseChimpFilter))
if(any(!sampleNames %in% samples$Code & !grepl('monkey[0-9]+',sampleNames)))stop('Problem assigning samples')
rownames(samples)<-samples$Code
metaData<-samples[sampleNames,]
rownames(metaData)<-names(condenseChimpFilter)
metaData$primer<-primers
metaData$primerBase<-sub('[0-9]+$','',metaData$primer)
metaData<-metaData[!grepl('monkey[0-9]+',sampleNames),]


library(ape)
source('16s/myBiplot.R')

colorBrew<-c('#e41a1cBB','#377eb8BB','#4daf4aBB','#984ea3BB','#ff7f00BB','#ffff33BB','#a65628BB','#f781bfBB','#999999BB')
nArea<-length(unique(metaData$area))
if(nArea>length(colorBrew))stop('Need to adjust colors for more areas')
areaCols<-colorBrew[1:nArea]
names(areaCols)<-unique(metaData$area)
areaPch<-sapply(names(areaCols),function(x)mostAbundant(metaData$Species[metaData$area==x]))
malariaCols<-c('#00000022','#000000CC')
malariaCols2<-rainbow.lab(2,alpha=.9,lightMultiple=.7)
primerCols<-rainbow.lab(length(unique(metaData$primer)))
names(primerCols)<-unique(metaData$primer)
names(malariaCols2)<-names(malariaCols)<-c('PlasmoNeg','PlasmoPos')
speciesPch<-20+1:length(unique(metaData$Species))
speciesCols<-rainbow.lab(length(unique(metaData$Species)),start=-2,end=2,alpha=.9,lightMultiple=.8)
names(speciesCols)<-names(speciesPch)<-unique(metaData$Species)

for(ii in unique(metaData$primerBase)){
  thisMetaData<-metaData[metaData$primerBase==ii,]
  uniPca<-pcoa(unifracDists[names(condenseChimpFilter) %in% rownames(thisMetaData),names(condenseChimpFilter) %in% rownames(thisMetaData)])
  predictors<-model.matrix(~0+Species+malaria+SIV+area+primer,thisMetaData)
  colnames(predictors)<-sub('^Species','',colnames(predictors))
  colnames(predictors)[colnames(predictors)=='malariaTRUE']<-'malariaPos'
  pdf(sprintf('out/%s_Pcoa.pdf',ii),height=9,width=11)
    mar<-c(4, 4, 1,9.5)
    sapply(list(1:2,3:4,5:6),function(axes){
      pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=speciesPch[thisMetaData$Species],bg=areaCols[thisMetaData$area],col=malariaCols[thisMetaData$malaria+1],cex=2.2,lwd=2.5,mar=mar)
      points(pos[thisMetaData$SIV=='Pos',],col='#FF000099',cex=2.7,lwd=2)
      legend(
        par('usr')[4]+.01*diff(par('usr')[3:4]),
        mean(par('usr')[3:4]),
        c(names(malariaCols),names(areaCols),names(speciesPch),'SIVPos'),
        col=c(malariaCols,rep(malariaCols,c(length(areaCols),length(speciesPch))),'#FF000099'),
        pch=c(rep(21,length(malariaCols)),speciesPch[areaPch],speciesPch,1),
        pt.bg=c(rep(NA,length(malariaCols)),areaCols,rep(NA,length(speciesPch)),NA),
        inset=.01,pt.lwd=3,pt.cex=2.5,xjust=0,xpd=NA)
      title(main=sprintf('All variables PC %d and %d',axes[1],axes[2]))
      #species
      pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=21,bg=speciesCols[thisMetaData$Species],col="#00000077",cex=1.8,lwd=2.5,mar=mar)
      legend(par('usr')[4]+.01*diff(par('usr')[3:4]), mean(par('usr')[3:4]),names(speciesCols),col='#00000077',pch=21,pt.bg=speciesCols,inset=c(-.01,.01),pt.lwd=3,pt.cex=2,xjust=0,xpd=NA)
      title(main=sprintf('Species PC %d and %d',axes[1],axes[2]))
      #malaria
      pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=21,bg=malariaCols2[thisMetaData$malaria+1],col="#00000077",cex=1.8,lwd=2.5,mar=mar)
      legend(par('usr')[4]+.01*diff(par('usr')[3:4]), mean(par('usr')[3:4]),names(malariaCols),col='#00000077',pch=21,pt.bg=malariaCols2,inset=.01,pt.lwd=3,pt.cex=2,xjust=0,xpd=NA)
      title(main=sprintf('Malaria PC %d and %d',axes[1],axes[2]))
      pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=21,bg=speciesCols[thisMetaData$Species],col=malariaCols[thisMetaData$malaria+1],cex=2.25,lwd=4,mar=mar)
      legend(par('usr')[4]+.01*diff(par('usr')[3:4]), mean(par('usr')[3:4]),as.vector(outer(names(speciesCols),names(malariaCols),paste,sep=' ')),col=as.vector(malariaCols[outer(names(speciesCols),names(malariaCols),function(x,y)y)]),pch=21,pt.bg=as.vector(speciesCols[outer(names(speciesCols),names(malariaCols),function(x,y)x)]),inset=.01,pt.lwd=4,pt.cex=2.5,xjust=0,xpd=NA)
      title(main=sprintf('Species/malaria PC %d and %d',axes[1],axes[2]))
      pos<-my.biplot.pcoa(uniPca,predictors,plot.axes=axes,pch=21,bg=primerCols[thisMetaData$primer],col=primerCols[thisMetaData$primer],cex=2.25,lwd=4,mar=mar)
      legend(par('usr')[4]+.01*diff(par('usr')[3:4]), mean(par('usr')[3:4]),names(primerCols),col=primerCols,pch=21,pt.bg=primerCols,inset=c(.01,-.01),pt.lwd=4,pt.cex=2.5,xjust=0,xpd=NA)
      title(main=sprintf('Primers PC %d and %d',axes[1],axes[2]))
      #show sample name text
      #bak<-par('cex')
      #par('cex'=.65)
      #biplot.pcoa(uniPca,predictors,plot.axes=axes)
      #title(main='Sample names')
      #par('cex'=bak)
    })
  dev.off()
}


