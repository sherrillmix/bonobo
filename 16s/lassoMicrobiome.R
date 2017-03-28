library(rpart)
library(randomForest)
library(glmnet)
library(glmnetPlotR)
if(!exists('otuProp'))source('runQiime.R')
tl<-samples[samples$isEnough&samples$isTL,]
tlOtu<-t(otuProp[,rownames(tl)])

fit<-cv.glmnet(tlOtu,tl$malaria,nfolds=nrow(tl),grouped=FALSE,type.measure='class',family='binomial')
coefs<-coef(fit)
coefs<-coefs[abs(coefs[,1])>0,]
coefs<-coefs[names(coefs)!="(Intercept)"]
taxa[names(coefs),]
pdf('out/lasso.pdf')
plotGlmnet(fit,markBest1SE=TRUE)
#plotBetas(fit$glmnet.fit,fit$lambda.1se)
for(ii in names(coefs)){
  #vpPlot(ifelse(tl$malaria,'Positive','Negative'),tlOtu[,ii],ylab=sprintf('Proportion OTU %s (%s)',ii,taxa[ii,'best']),bg=rainbow.lab(2)[tl$malaria+1],pch=21,cex=2,offsetXArgs=list(width=.25),xlim=c(.5,2.5))
  labels<-withAs(s=samples[samples$isEnough,],sprintf('%s %s%s',ifelse(s$bonobo,'Bonobo','Chimp'),ifelse(s$isTL,'TL ',ifelse(s$bonobo,'non-endemic ','')),ifelse(s$malaria,'Positive','Negative')))
  par(mar=c(13,4.5,.1,.1))
  vpPlot(labels,otuProp[ii,rownames(samples)[samples$isEnough]],ylab=sprintf('Proportion OTU %s (%s)',ii,taxa[ii,'best']),bg=rainbow.lab(2)[samples[samples$isEnough,'malaria']+1],pch=21,cex=2,offsetXArgs=list(width=.4),las=2,mgp=c(3.3,1,0))
}
dev.off()

forest<-rpart(tl$malaria~.,data=as.data.frame(tlOtu[,apply(tlOtu,2,max)>.01]))
forest2<-randomForest(as.factor(tl$malaria)~.,data=as.data.frame(tlOtu[,apply(tlOtu,2,max)>.0001]),importance=TRUE)
pdf('test.pdf');varImpPlot(forest2);dev.off()
