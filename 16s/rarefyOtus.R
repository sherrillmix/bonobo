if(!exists('otuTab'))source('runQiime.R')


medianCI<-function(xx,na.rm=TRUE){
  if(na.rm)xx<-xx[!is.na(xx)]
  return(sort(xx)[qbinom(c(.025,.975), length(xx), 0.5)])
}

rareSteps<-seq(100,15000,100)
rareCurves<-apply(otuTab[,samples$name[samples$isEnough]],2,rareEquation,rareSteps)
groupings<-paste(samples$Species[samples$isEnough],ifelse(samples$malaria[samples$isEnough],'Plasmodium positive','Plasmodium negative'))
groupCurves<-lapply(unique(groupings),function(xx)t(apply(rareCurves[,groupings==xx],1,quantile,c(.975,.5,.025),na.rm=TRUE)))
groupCI<-lapply(unique(groupings),function(xx)t(apply(rareCurves[,groupings==xx],1,medianCI)))
names(groupCurves)<-names(groupCI)<-unique(groupings)
groupCols<-rainbow.lab(length(groupCurves),alpha=.9)
groupCols2<-rainbow.lab(length(groupCurves),alpha=.2)
names(groupCols)<-names(groupCols2)<-unique(groupings)
pdf('out/rare.pdf',width=6,height=6)
  par(mar=c(4,4,.1,.2))
  plot(1,1,type='n',xlim=range(rareSteps),ylim=range(do.call(rbind,groupCurves)),xlab='Number of reads',ylab='Number of OTUs',las=1)
  for(ii in names(groupCurves)){
    polygon(c(rareSteps,rev(rareSteps)),c(groupCI[[ii]][,1],rev(groupCI[[ii]][,2])),col=groupCols2[ii],border=NA)
  }
  for(ii in names(groupCurves)){
    #polygon(c(rareSteps,rev(rareSteps)),c(groupCurves[[ii]][,'97.5%'],rev(groupCurves[[ii]][,'2.5%'])),col=groupCols2[ii],border=NA)
    lines(rareSteps,groupCurves[[ii]][,'50%'],col=groupCols[ii],lwd=3)
  }
  legend('topleft',names(groupCols),lwd=2,col=groupCols,inset=.01)
dev.off()



moreThanProp<-apply(otuProp[,samples$name[samples$isEnough]],1,max,na.rm=TRUE)>.02
cols<-c('white',tail(rev(heat.colors(110)),99)) 
plotProp<-t(otuProp[moreThanProp,samples$name[samples$isEnough]])
breaks<-c(-1e-6,seq(min(plotProp[plotProp>0])-1e-10,max(plotProp)+1e-10,length.out=100))
rownames(plotProp)<-sub("EasternChimpanzee","Chimp",rownames(plotProp))
colnames(plotProp)<-taxa[colnames(plotProp),'best']
pdf('out/heat.pdf',height=20,width=20)
  heatmap(plotProp,col=cols,breaks=breaks,scale='none',margin=c(8,5),Rowv=NA)
  insetScale(round(breaks,6),cols,c(.955,.75,.97,.99),label='Proportion')
dev.off()
