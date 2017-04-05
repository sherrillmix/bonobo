library(dnar)
if(!exists('otuTab'))source('runQiime.R')

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
  insetScale(round(breaks2,6),cols,c(.97,.75,.98,.99),label='Proportion of OTU within each sample')
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
colnames(selectProp)<-sprintf('%s (q=%0.3f)',taxa[colnames(selectProp),'bestId'],pVals[colnames(selectProp)])
rownames(selectProp)<-sprintf('%s%s',ifelse(samples[rownames(selectProp),'malaria'],'+','-'),sub("EasternChimpanzee","Chimp",rownames(selectProp)))
#colnames(selectProp)<-sprintf('%s (%s)\nq=%0.3f',colnames(selectProp),naReplace(taxa[colnames(selectProp),'best'],'Unknown'),pVals[pVals<pCut])

allPos<-samples$name[samples$isEnough&samples$malaria]
allNeg<-samples$name[samples$isEnough&!samples$malaria]
inAll<-otuProp[apply(otuProp[,c(allPos,allPos)],1,max)>.001,]
pValsAll<-p.adjust(apply(inAll,1,function(xx)wilcox.test(xx[allPos],xx[allNeg])$p.value),'fdr')
pValsAll[is.na(pValsAll)]<-1
selectPropAll<-apply(inAll[pValsAll<pCut,samples$name[samples$isEnough]],1,function(x)x/max(x))
rownames(selectPropAll)<-sprintf('%s%s',ifelse(samples[rownames(selectPropAll),'malaria'],'+','-'),sub("EasternChimpanzee","Chimp",rownames(selectPropAll)))
colnames(selectPropAll)<-sprintf('%s (%s)\nq=%0.3f',colnames(selectPropAll),taxa[colnames(selectPropAll),'bestId'],pValsAll[colnames(selectPropAll)])

pCutTl<-.01
tls<-samples$name[samples$isEnough&samples$isTL&samples$bonobo]
nonTls<-samples$name[samples$isEnough&!samples$isTL&samples$bonobo]
pValsTl<-p.adjust(apply(inAll,1,function(xx)wilcox.test(xx[tls],xx[nonTls])$p.value),'fdr')
pValsTl[is.na(pValsTl)]<-1
selectPropTl<-apply(inAll[pValsTl<pCutTl,samples$name[samples$isEnough]],1,function(x)x/max(x))
hTree<-hclust(dist(t(selectPropTl[c(tls,nonTls),]^.25)))
selectPropTl<-selectPropTl[,hTree$labels[hTree$order]]
rownames(selectPropTl)<-sprintf('%s%s',ifelse(samples[rownames(selectPropTl),'malaria'],'+','-'),sub("EasternChimpanzee","Chimp",rownames(selectPropTl)))
colnames(selectPropTl)<-sprintf('%s (q=%0.3f)',taxa[colnames(selectPropTl),'bestId'],pValsTl[colnames(selectPropTl)])

iks<-samples$name[samples$isEnough&samples$area=='IK'&samples$bonobo]
nonIks<-samples$name[samples$isEnough&samples$area!='IK'&samples$bonobo]
pValsIk<-p.adjust(apply(inAll,1,function(xx)wilcox.test(xx[iks],xx[nonIks])$p.value),'fdr')
pValsIk[is.na(pValsIk)]<-1
selectPropIk<-apply(inAll[pValsIk<pCut,samples$name[samples$isEnough]],1,function(x)x/max(x))
hTree<-hclust(dist(t(selectPropIk[c(iks,nonIks),]^.25)))
selectPropIk<-selectPropIk[,hTree$labels[hTree$order]]
rownames(selectPropIk)<-sprintf('%s%s',ifelse(samples[rownames(selectPropIk),'malaria'],'+','-'),sub("EasternChimpanzee","Chimp",rownames(selectPropIk)))
colnames(selectPropIk)<-sprintf('%s (q=%0.3f)',taxa[colnames(selectPropIk),'bestId'],pValsIk[colnames(selectPropIk)])


plotHeat<-function(selectProp,breaks,cols,xaxt=''){
  image(1:ncol(selectProp),1:nrow(selectProp),t(selectProp),col=cols,breaks=breaks,xaxt='n',yaxt='n',xlab='',ylab='')
  box()
  insetScale(round(breaks,6),cols,c(.97,.75,.98,.99),label='Proportion of OTU maximum')
  if(xaxt!='n')slantAxis(1,1:ncol(selectProp),colnames(selectProp))
  axis(4,1:nrow(selectProp),rownames(selectProp),las=1,tcl=0,mgp=c(0,.2,0),cex.axis=.7)
  abline(h=1:nrow(selectProp)-.5,v=1:ncol(selectProp)+.5,col='#00000011')
}


pdf('out/splitOtus.pdf',height=13,width=12)
  par(mar=c(11.5,.1,3,8.5),lheight=.7)
  plotHeat(selectProp,breaks2,cols)
  title(main='TL + vs TL - (q<.1)')
  plotHeat(selectPropTl,breaks2,cols,xaxt='n')
  axis(1,1:ncol(selectPropTl),colnames(selectPropTl),cex.axis=.7,las=2)
  title(main='TL vs nonTL (q<.01)')
  plotHeat(selectPropIk,breaks2,cols,xaxt='n')
  axis(1,1:ncol(selectPropIk),colnames(selectPropIk),cex.axis=.7,las=2)
  title(main='IK vs nonIK (q<.1)')
  #heatmap(selectPropTl,col=cols,breaks=breaks2,scale='none',margin=c(16,6),Rowv=NA)
  #insetScale(round(breaks2,6),cols,c(.955,.75,.97,.99),label='Proportion of OTU maximum')
dev.off()


