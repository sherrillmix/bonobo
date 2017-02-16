if(!exists('otuTab'))source('runQiime.R')

samples$isEnough<-apply(otuTab[,samples$name],2,sum)>20000
tail(sort(apply(otuTab[,samples$name],2,sum)))

medianCI<-function(xx,na.rm=TRUE){
  if(na.rm)xx<-xx[!is.na(xx)]
  return(sort(xx)[qbinom(c(.025,.975), length(xx), 0.5)])
}

rareSteps<-seq(100,20000,100)
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
