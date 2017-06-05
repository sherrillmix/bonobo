library(dnar)
if(!exists('otuTab'))source('runQiime.R')


rareSteps<-seq(100,15000,100)
rareCurves<-apply(otuTab[,samples$name[samples$isEnough]],2,rareEquation,rareSteps)

groupings<-paste(ifelse(samples$bonobo[samples$isEnough],'Bonobo','Chimpanzee'),ifelse(samples$malaria[samples$isEnough],'Laverania positive','Laverania negative'))
groupings2<-paste(ifelse(samples$bonobo[samples$isEnough],'Bonobo ','Chimpanzee '),ifelse(grepl('^TL',samples$area[samples$isEnough]),'TL2 site',ifelse(samples$bonobo[samples$isEnough],'non-endemic sites','all sites')),sep='')
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
pdf('out/Fig.S9A.pdf',width=8,height=6)
  spacer<-.5
  par(mar=c(6.6,4,.3,5),lheight=.85)
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
  slantAxis(1,1:length(levels(groupFac)),sub('(Chimpanzee|Bonobo)','\\1 samples',sub(' Laverania','\nLaverania',levels(groupFac))),srt=-30)
  groupFac2<-factor(groupings2,levels=unique(groupings2))
  #xPos<-as.numeric(groupFac2)+ave(shannons,groupFac2,FUN=function(xx)swarmx(0,xx,cex=1.9)$x)
  offset<-max(as.numeric(groupFac))+spacer
  abline(v=offset+.5-spacer/2,lty=2)
  xPos<-offset+as.numeric(groupFac2)+offsetX(shannons,groupFac2,width=.3)
  width<-.45
  segments(offset+1:length(levels(groupFac2))-width,group2Shannon[sub('\n',' ',levels(groupFac2))],offset+1:length(levels(groupFac2))+width,group2Shannon[sub('\n',' ',levels(groupFac2))],lwd=3,col=group2Cols[sub('\n',' ',levels(groupFac2))])
  rect(offset+1:length(levels(groupFac2))-width,shannon2CI[1,sub('\n',' ',levels(groupFac2))],offset+1:length(levels(groupFac2))+width,shannon2CI[2,sub('\n',' ',levels(groupFac2))],lwd=2,border=NA,col=group2Cols2[sub('\n',' ',levels(groupFac2))])
  points(xPos,shannons,pch=21,bg=group2Cols[groupings2],cex=1.7)
  #replacing first space with \n
  slantAxis(1,offset+1:length(levels(groupFac2)),sub('(Chimpanzee|Bonobo)','\\1 samples',sub(' ','\n',levels(groupFac2))),srt=-30)
dev.off()

