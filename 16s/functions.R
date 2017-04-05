naReplace<-function(x,replace){x[is.na(x)]<-replace;return(x)}

plotHeat<-function(selectProp,breaks,cols,xaxt=''){
  image(1:ncol(selectProp),1:nrow(selectProp),t(selectProp),col=cols,breaks=breaks,xaxt='n',yaxt='n',xlab='',ylab='')
  box()
  insetScale(round(breaks,6),cols,c(.97,.75,.98,.99),label='Proportion of OTU maximum')
  if(xaxt!='n')slantAxis(1,1:ncol(selectProp),colnames(selectProp))
  axis(4,1:nrow(selectProp),rownames(selectProp),las=1,tcl=0,mgp=c(0,.2,0),cex.axis=.7)
  abline(h=1:nrow(selectProp)-.5,v=1:ncol(selectProp)+.5,col='#00000011')
}


fishers<-function(ps,correct=1){
  chi2<- -2*sum(log(ps))
  out<-1-pchisq(chi2/correct,2*length(ps)/correct)
  return(out)
}
