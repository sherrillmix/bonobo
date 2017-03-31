library(dnar)
if(!exists('otuTab'))source('runQiime.R')


medianCI<-function(xx,na.rm=TRUE){
  if(na.rm)xx<-xx[!is.na(xx)]
  return(unname(sort(xx)[qbinom(c(.025,.975), length(xx), 0.5)]))
}

rareSteps<-seq(100,15000,100)
rareCurves<-apply(otuTab[,samples$name[samples$isEnough]],2,rareEquation,rareSteps)

groupings<-paste(samples$Species[samples$isEnough],ifelse(samples$malaria[samples$isEnough],'Plasmodium positive','Plasmodium negative'))
groupings2<-paste(ifelse(grepl('^TL',samples$area[samples$isEnough]),'TL2 ',ifelse(samples$bonobo[samples$isEnough],'Non-endemic ','')),samples$Species[samples$isEnough],sep='')
groupCurves<-lapply(unique(groupings),function(xx)t(apply(rareCurves[,groupings==xx],1,quantile,c(.975,.5,.025),na.rm=TRUE)))
group2Curves<-lapply(unique(groupings2),function(xx)t(apply(rareCurves[,groupings2==xx],1,quantile,c(.975,.5,.025),na.rm=TRUE)))
groupCI<-lapply(unique(groupings),function(xx)t(apply(rareCurves[,groupings==xx],1,medianCI)))
group2CI<-lapply(unique(groupings2),function(xx)t(apply(rareCurves[,groupings2==xx],1,medianCI)))
names(groupCurves)<-names(groupCI)<-unique(groupings)
names(group2Curves)<-names(group2CI)<-unique(groupings2)

#mildly magic number throwing out 3rd color to leave gap in coloring between two groups
groupCols<-rainbow.lab(length(groupCurves)+length(group2Curves)+1,alpha=.9)[1:(length(groupCurves)+1)][-3]
groupCols2<-rainbow.lab(length(groupCurves)+length(group2Curves)+1,alpha=.25)[1:(length(groupCurves)+1)][-3]
names(groupCols)<-names(groupCols2)<-unique(groupings)
group2Cols<-rainbow.lab(length(groupCurves)+length(group2Curves)+1,alpha=.9)[1+length(groupCurves)+1:length(group2Curves)]
group2Cols2<-rainbow.lab(length(groupCurves)+length(group2Curves)+1,alpha=.25)[1+length(groupCurves)+1:length(group2Curves)]
names(group2Cols)<-names(group2Cols2)<-unique(groupings2)


pdf('out/rare.pdf',width=6,height=6)
  par(mar=c(4,4,.1,.2))
  #species and malaria
  plot(1,1,type='n',xlim=range(rareSteps),ylim=range(do.call(rbind,groupCurves)),xlab='Number of reads',ylab='Number of OTUs',las=1)
  for(ii in names(groupCurves)){
    polygon(c(rareSteps,rev(rareSteps)),c(groupCI[[ii]][,1],rev(groupCI[[ii]][,2])),col=groupCols2[ii],border=NA)
  }
  for(ii in names(groupCurves)){
    lines(rareSteps,groupCurves[[ii]][,'50%'],col=groupCols[ii],lwd=3)
  }
  legend('topleft',names(groupCols),lwd=2,col=groupCols,inset=.01)
  #TL vs rest
  plot(1,1,type='n',xlim=range(rareSteps),ylim=range(do.call(rbind,group2Curves)),xlab='Number of reads',ylab='Number of OTUs',las=1)
  for(ii in names(group2Curves)){
    polygon(c(rareSteps,rev(rareSteps)),c(group2CI[[ii]][,1],rev(group2CI[[ii]][,2])),col=group2Cols2[ii],border=NA)
  }
  for(ii in names(group2Curves)){
    lines(rareSteps,group2Curves[[ii]][,'50%'],col=group2Cols[ii],lwd=3)
  }
  legend('topleft',names(group2Cols),lwd=2,col=group2Cols,inset=.01)
dev.off()

shannons<-apply(otuTab[,samples$name[samples$isEnough]],2,shannon)
groupShannon<-sapply(unique(groupings),function(xx)median(shannons[groupings==xx]))
shannonCI<-sapply(unique(groupings),function(xx)medianCI(shannons[groupings==xx]))
group2Shannon<-sapply(unique(groupings2),function(xx)median(shannons[groupings2==xx]))
shannon2CI<-sapply(unique(groupings2),function(xx)medianCI(shannons[groupings2==xx]))

outer(split(shannons,groupings),split(shannons,groupings),function(xx,yy)mapply(function(xxx,yyy)wilcox.test(xxx,yyy)$p.value,xx,yy))
outer(split(shannons,groupings2),split(shannons,groupings2),function(xx,yy)mapply(function(xxx,yyy)wilcox.test(xxx,yyy)$p.value,xx,yy))
outer(split(shannons,sub('Plasmodi.*','',groupings)),split(shannons,sub('Plasmodi.*','',groupings)),function(xx,yy)mapply(function(xxx,yyy)wilcox.test(xxx,yyy)$p.value,xx,yy))

#library(beeswarm)
pdf('out/shannon.pdf',width=8,height=6)
  spacer<-.5
  par(mar=c(6.6,4,.1,3),lheight=.85)
  plot(1,1,type='n',ylab='Shannon diversity',las=2,xlim=c(.5,length(unique(c(groupings,groupings2)))+.5+spacer),ylim=range(shannons),xaxt='n',xlab='',bty='l')
  groupFac<-factor(sub(' Plas','\nPlas',groupings))
  #xPos<-as.numeric(groupFac)+ave(shannons,groupFac,FUN=function(xx)swarmx(0,xx,cex=1.9)$x)
  xPos<-as.numeric(groupFac)+offsetX(shannons,groupFac,width=.3)
  width<-.45
  segments(1:length(levels(groupFac))-width,groupShannon[sub('\n',' ',levels(groupFac))],1:length(levels(groupFac))+width,groupShannon[sub('\n',' ',levels(groupFac))],lwd=3,col=groupCols[sub('\n',' ',levels(groupFac))])
  #rect(1:length(levels(groupFac))-width,shannonCI[1,sub('\n',' ',levels(groupFac))],1:length(levels(groupFac))+width,shannonCI[2,sub('\n',' ',levels(groupFac))],lwd=2,col=NA,border=groupCols[sub('\n',' ',levels(groupFac))])
  rect(1:length(levels(groupFac))-width,shannonCI[1,sub('\n',' ',levels(groupFac))],1:length(levels(groupFac))+width,shannonCI[2,sub('\n',' ',levels(groupFac))],lwd=2,border=NA,col=groupCols2[sub('\n',' ',levels(groupFac))])
  points(xPos,shannons,pch=21,bg=groupCols[groupings],cex=1.7)
    #,shannons,,pch=21,bg=groupCols[groupings],cex=1.5
  #vpPlot(factor(groupings2,levels=unique(groupings2)),shannons,ylab='Shannon diversity',las=2,pch=21,bg=group2Cols[groupings2],cex=1.5)
  slantAxis(1,1:length(levels(groupFac)),levels(groupFac),srt=-30)
  groupFac2<-factor(groupings2,levels=unique(groupings2))
  #xPos<-as.numeric(groupFac2)+ave(shannons,groupFac2,FUN=function(xx)swarmx(0,xx,cex=1.9)$x)
  offset<-max(as.numeric(groupFac))+spacer
  abline(v=offset+.5-spacer/2,lty=2)
  xPos<-offset+as.numeric(groupFac2)+offsetX(shannons,groupFac2,width=.3)
  width<-.45
  segments(offset+1:length(levels(groupFac2))-width,group2Shannon[sub('\n',' ',levels(groupFac2))],offset+1:length(levels(groupFac2))+width,group2Shannon[sub('\n',' ',levels(groupFac2))],lwd=3,col=group2Cols[sub('\n',' ',levels(groupFac2))])
  rect(offset+1:length(levels(groupFac2))-width,shannon2CI[1,sub('\n',' ',levels(groupFac2))],offset+1:length(levels(groupFac2))+width,shannon2CI[2,sub('\n',' ',levels(groupFac2))],lwd=2,border=NA,col=group2Cols2[sub('\n',' ',levels(groupFac2))])
  points(xPos,shannons,pch=21,bg=group2Cols[groupings2],cex=1.7)
  slantAxis(1,offset+1:length(levels(groupFac2)),levels(groupFac2),srt=-30)
dev.off()

naReplace<-function(x,replace){x[is.na(x)]<-replace;return(x)}

samples<-samples[order(samples$Species,samples$area2,samples$malaria,samples$Code),]
moreThanProp<-apply(otuProp[,samples$name[samples$isEnough]],1,max,na.rm=TRUE)>.02
cols<-c('white',tail(rev(heat.colors(110)),99)) 
plotProp<-t(otuProp[moreThanProp,samples$name[samples$isEnough]])
maxProp<-apply(otuProp[moreThanProp,samples$name[samples$isEnough]],1,function(x)x/max(x))
rownames(plotProp)<-rownames(maxProp)<-sprintf('%s%s',ifelse(samples[rownames(plotProp),'malaria'],'+','-'),sub("EasternChimpanzee","Chimp",rownames(plotProp)))
colnames(plotProp)<-colnames(maxProp)<-taxa[colnames(plotProp),'bestId']
#plotTree<-hclust(dist(t(plotProp[,])))
#plotProp<-plotProp[,plotTree$labels[plotTree$order]]
plotProp<-plotProp[,order(apply(plotProp,2,mean),decreasing=TRUE)]
maxTree<-hclust(dist(t(maxProp[,])))
maxProp<-maxProp[,rev(maxTree$labels[maxTree$order])]
breaks<-c(-1e-6,seq(min(plotProp[plotProp>0])-1e-10,max(plotProp)+1e-10,length.out=100))
breaks2<-c(-1e-6,seq(min(maxProp[maxProp>0])-1e-10,max(maxProp)+1e-10,length.out=100))
#plotProp<-plotProp[,order(apply(plotProp,2,mean))]
pdf('out/heat.pdf',height=30,width=20)
  par(mfrow=c(2,1),mar=c(12,5,3.5,.1))
  image(1:ncol(plotProp),1:nrow(plotProp),t(plotProp),col=cols,breaks=breaks,xlab='',ylab='',xaxt='n',yaxt='n')
  text(grconvertX(.005, "nfc", "user"),grconvertY(.995, "nfc", "user"),'A)',xpd=NA,cex=3,adj=0:1)
  box()
  insetScale(round(breaks2,6),cols,c(.97,.75,.98,.99),label='Proportion of sample')
  axis(1,1:ncol(plotProp),colnames(plotProp),cex.axis=.7,las=2,tcl=-.1,mgp=c(0,.3,0))
  #slantAxis(1,1:ncol(plotProp),colnames(plotProp),cex=.7)
  axis(2,1:nrow(plotProp),rownames(plotProp),las=1,tcl=0,mgp=c(0,.2,0),cex.axis=.65)
  abline(h=1:nrow(plotProp)-.5,v=1:ncol(plotProp)+.5,col='#00000011')
  image(1:ncol(maxProp),1:nrow(maxProp),t(maxProp),col=cols,breaks=breaks2,xlab='',ylab='',xaxt='n',yaxt='n')
  text(grconvertX(.005, "nfc", "user"),grconvertY(.995, "nfc", "user"),'B)',xpd=NA,cex=3,adj=0:1)
  box()
  insetScale(round(breaks2,6),cols,c(.97,.75,.98,.99),label='Proportion of OTU maximum')
  axis(1,1:ncol(maxProp),colnames(maxProp),cex.axis=.7,las=2,tcl=-.1,mgp=c(0,.3,0))
  axis(2,1:nrow(maxProp),rownames(maxProp),las=1,tcl=0,mgp=c(0,.2,0),cex.axis=.65)
  abline(h=1:nrow(maxProp)-.5,v=1:ncol(maxProp)+.5,col='#00000011')
  #heatmap(maxProp,col=cols,breaks=breaks2,scale='none',margin=c(8,5),Rowv=NA)
dev.off()


tlPos<-samples$name[samples$isEnough&samples$isTL&samples$malaria]
tlNeg<-samples$name[samples$isEnough&samples$isTL&!samples$malaria]
inTl<-otuProp[apply(otuProp[,c(tlPos,tlNeg)],1,max)>.001,]
pVals<-p.adjust(apply(inTl,1,function(xx)wilcox.test(xx[tlPos],xx[tlNeg])$p.value),'fdr')
effect<-p.adjust(apply(inTl,1,function(xx)wilcox.test(xx[tlPos],xx[tlNeg])$p.value),'fdr')
pVals[is.na(pVals)]<-1
pCut<-.1
selectProp<-apply(inTl[pVals<pCut,samples$name[samples$isEnough]],1,function(x)x/max(x))
hTree<-hclust(dist(t(selectProp[c(tlPos,tlNeg),])))
selectProp<-selectProp[,hTree$labels[hTree$order]]
colnames(selectProp)<-sprintf('%s (q=%0.3f)',taxa[colnames(selectProp),'bestId'],pVals[pVals<pCut])
rownames(selectProp)<-sprintf('%s%s',ifelse(samples[rownames(selectProp),'malaria'],'+','-'),sub("EasternChimpanzee","Chimp",rownames(selectProp)))
#colnames(selectProp)<-sprintf('%s (%s)\nq=%0.3f',colnames(selectProp),naReplace(taxa[colnames(selectProp),'best'],'Unknown'),pVals[pVals<pCut])

allPos<-samples$name[samples$isEnough&samples$malaria]
allNeg<-samples$name[samples$isEnough&!samples$malaria]
inAll<-otuProp[apply(otuProp[,c(allPos,allPos)],1,max)>.001,]
pValsAll<-p.adjust(apply(inAll,1,function(xx)wilcox.test(xx[allPos],xx[allNeg])$p.value),'fdr')
pValsAll[is.na(pValsAll)]<-1
selectPropAll<-apply(inAll[pValsAll<pCut,samples$name[samples$isEnough]],1,function(x)x/max(x))
rownames(selectPropAll)<-sprintf('%s%s',ifelse(samples[rownames(selectPropAll),'malaria'],'+','-'),sub("EasternChimpanzee","Chimp",rownames(selectPropAll)))
colnames(selectPropAll)<-sprintf('%s (%s)\nq=%0.3f',colnames(selectPropAll),taxa[colnames(selectPropAll),'bestId'],pValsAll[pValsAll<pCut])


pdf('out/splitOtus.pdf',height=13,width=12)
  par(mar=c(11.5,.1,3,8.5),lheight=.7)
  image(1:ncol(selectProp),1:nrow(selectProp),t(selectProp),col=cols,breaks=breaks2,xaxt='n',yaxt='n',xlab='',ylab='')
  box()
  insetScale(round(breaks2,6),cols,c(.97,.75,.98,.99),label='Proportion of OTU maximum')
  slantAxis(1,1:ncol(selectProp),colnames(selectProp))
  axis(4,1:nrow(selectProp),rownames(selectProp),las=1,tcl=0,mgp=c(0,.2,0),cex.axis=.7)
  abline(h=1:nrow(selectProp)-.5,v=1:ncol(selectProp)+.5,col='#00000011')
  #heatmap(selectPropAll,col=cols,breaks=breaks2,scale='none',margin=c(16,6),Rowv=NA)
  #insetScale(round(breaks2,6),cols,c(.955,.75,.97,.99),label='Proportion of OTU maximum')
dev.off()


