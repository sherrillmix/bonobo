library(dnar)
if(!exists('otuTab'))source('runQiime.R')


samples<-samples[order(samples$Species,samples$area2,samples$malaria,samples$Code),]
moreThanProp<-apply(otuProp[,samples$name[samples$isEnough]],1,max,na.rm=TRUE)>.02
cols<-c('white',tail(rev(heat.colors(110)),99)) 
plotProp<-t(otuProp[moreThanProp,samples$name[samples$isEnough]])
maxProp<-apply(otuProp[moreThanProp,samples$name[samples$isEnough]],1,function(x)x/max(x))
#rownames(plotProp)<-rownames(maxProp)<-sprintf('%s%s',ifelse(samples[rownames(plotProp),'malaria'],'+','-'),sub("EasternChimpanzee","Chimp",rownames(plotProp)))
colnames(plotProp)<-colnames(maxProp)<-taxa[colnames(plotProp),'bestId']
#plotTree<-hclust(dist(t(plotProp[,])))
#plotProp<-plotProp[,plotTree$labels[plotTree$order]]
plotProp<-plotProp[,order(apply(plotProp,2,mean),decreasing=TRUE)]
maxTree<-hclust(dist(t(maxProp[,])))
maxProp<-maxProp[,rev(maxTree$labels[maxTree$order])]
plotProp<-plotProp[,rev(maxTree$labels[maxTree$order])]
breaks<-c(-1e-6,seq(min(plotProp[plotProp>0])-1e-10,max(plotProp)+1e-10,length.out=100))
breaks2<-c(-1e-6,seq(min(maxProp[maxProp>0])-1e-10,max(maxProp)+1e-10,length.out=100))
#plotProp<-plotProp[,order(apply(plotProp,2,mean))]
pdf('out/Fig.S8.pdf',height=30,width=20)
  par(mfrow=c(2,1),mar=c(12,.1,3.5,14))
  image(1:ncol(plotProp),1:nrow(plotProp),t(plotProp),col=cols,breaks=breaks,xlab='',ylab='',xaxt='n',yaxt='n')
  text(grconvertX(.005, "nfc", "user"),grconvertY(.995, "nfc", "user"),'A)',xpd=NA,cex=3,adj=0:1)
  box()
  insetScale(round(breaks2,6),cols,c(.97,.58,.98,.83),label='Proportion of OTU within each sample')
  axis(1,1:ncol(plotProp),colnames(plotProp),cex.axis=.7,las=2,tcl=-.1,mgp=c(0,.3,0))
  #slantAxis(1,1:ncol(plotProp),colnames(plotProp),cex=.7)
  #axis(2,1:nrow(plotProp),rownames(plotProp),las=1,tcl=0,mgp=c(0,.2,0),cex.axis=.65)
  metadata<-samples[rownames(plotProp),c('chimpBonobo','area2','plasmoPM','Code')]
  colnames(metadata)<-c('Species','Site','Laverania','Sample')
  addMetaData(metadata,cex=.75)
  abline(h=1:nrow(plotProp)-.5,v=1:ncol(plotProp)+.5,col='#00000011')
  image(1:ncol(maxProp),1:nrow(maxProp),t(maxProp),col=cols,breaks=breaks2,xlab='',ylab='',xaxt='n',yaxt='n')
  text(grconvertX(.005, "nfc", "user"),grconvertY(.995, "nfc", "user"),'B)',xpd=NA,cex=3,adj=0:1)
  box()
  insetScale(round(breaks2,6),cols,c(.97,.58,.98,.83),label='Proportion of OTU maximum')
  axis(1,1:ncol(maxProp),colnames(maxProp),cex.axis=.7,las=2,tcl=-.1,mgp=c(0,.3,0))
  #axis(2,1:nrow(maxProp),rownames(maxProp),las=1,tcl=0,mgp=c(0,.2,0),cex.axis=.65)
  abline(h=1:nrow(maxProp)-.5,v=1:ncol(maxProp)+.5,col='#00000011')
  metadata<-samples[rownames(maxProp),c('chimpBonobo','area2','plasmoPM','Code')]
  colnames(metadata)<-c('Species','Site','Laverania','Sample')
  addMetaData(metadata,cex=.75)
  #heatmap(maxProp,col=cols,breaks=breaks2,scale='none',margin=c(8,5),Rowv=NA)
dev.off()


sampleOrder<-withAs(s=samples[samples$isEnough,],s$name[order(!s$bonobo,ifelse(s$area=='KR',FALSE,s$malaria),s$area2,s$malaria)])
sampleOrder2<-withAs(s=samples[samples$isEnough,],s$name[order(!s$bonobo,s$area2,s$malaria)])

tlPos<-samples$name[samples$isEnough&samples$isTL&samples$malaria]
tlNeg<-samples$name[samples$isEnough&samples$isTL&!samples$malaria]
tmp<-setupHeat(tlPos,tlNeg,otuProp[,sampleOrder],taxa,pCut=.05)
selectProp<-tmp[[1]]
effectSplit<-tmp[[2]]
#
allPos<-samples$name[samples$isEnough&samples$malaria]
allNeg<-samples$name[samples$isEnough&!samples$malaria]
tmp<-setupHeat(allPos,allNeg,otuProp[,sampleOrder],taxa)
selectPropAll<-tmp[[1]]
effectSplitAll<-tmp[[2]]
#
pCutTl<-.01
tls<-samples$name[samples$isEnough&samples$isTL&samples$bonobo]
nonTls<-samples$name[samples$isEnough&!samples$isTL&samples$bonobo]
tmp<-setupHeat(tls,nonTls,otuProp[,sampleOrder2],taxa,pCut=pCutTl)
selectPropTl<-tmp[[1]]
effectSplitTl<-tmp[[2]]
#
iks<-samples$name[samples$isEnough&samples$area=='IK'&samples$bonobo]
nonIks<-samples$name[samples$isEnough&samples$area!='IK'&samples$bonobo]
tmp<-setupHeat(iks,nonIks,otuProp[,sampleOrder2],taxa)
selectPropIk<-tmp[[1]]
effectSplitIk<-tmp[[2]]
#
pdf('out/splitOtus.pdf',height=13,width=15)
  par(mar=c(11.5,.1,3,13.5),lheight=.7)
  plotHeat(selectProp,breaks2,cols,yaxt='n')
  #title(main='TL + vs TL - (q<.05)')
  abline(v=effectSplit-.5)
  metadata<-samples[rownames(selectProp),c('chimpBonobo','area2','plasmoPM','Code')]
  colnames(metadata)<-c('Species','Site','Laverania','Sample')
  addMetaData(metadata,cex=.75)
  plotHeat(selectPropTl,breaks2,cols,xaxt='n',yaxt='n')
  axis(1,1:ncol(selectPropTl),colnames(selectPropTl),cex.axis=.7,las=2)
  title(main='TL vs nonTL (q<.01)')
  abline(v=effectSplitTl-.5)
  metadata<-samples[rownames(selectPropTl),c('chimpBonobo','area2','plasmoPM','Code')]
  colnames(metadata)<-c('Species','Site','Laverania','Sample')
  addMetaData(metadata,cex=.75)
  plotHeat(selectPropIk,breaks2,cols,xaxt='n',yaxt='n')
  axis(1,1:ncol(selectPropIk),colnames(selectPropIk),cex.axis=.7,las=2)
  title(main='IK vs nonIK (q<.1)')
  addMetaData(metadata,cex=.75)
  abline(v=effectSplitIk-.5)
  #heatmap(selectPropTl,col=cols,breaks=breaks2,scale='none',margin=c(16,6),Rowv=NA)
  #insetScale(round(breaks2,6),cols,c(.955,.75,.97,.99),label='Proportion of OTU maximum')
dev.off()
system('pdftk out/splitOtus.pdf cat 1 output out/Fig.S10.pdf')


