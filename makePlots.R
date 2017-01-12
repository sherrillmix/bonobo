
if(!exists('taxas'))source('parseBlast.R')

chimpFilter<-lapply(taxas,function(x)x[x[,'class']!='Mammalia'|is.na(x[,'class']),])
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
  for(ii in types2){
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
