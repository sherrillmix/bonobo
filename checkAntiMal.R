library(taxonomizr)
source('readSamples.R')
library(vipor)

antiMal<-read.csv('AntiMalPlant_List_021817.csv',stringsAsFactors=FALSE)
getNamesAndNodes('~/db/taxo/')
taxaNodes<-read.nodes('~/db/taxo/nodes.dmp')
taxaNamesAmbig<-read.names('~/db/taxo/names.dmp',FALSE)
taxaNames<-read.names('~/db/taxo/names.dmp')
antiMal$clean<-sub(' sp$','',sub('([A-Za-z]+ [a-z]+).*','\\1',antiMal$Plant))
antiMal$clean[antiMal$clean=="AAloe parvibracteata"]<-"Aloe parvibracteata"
antiMal$clean[antiMal$clean=="Curcuma aromatic"]<-"Curcuma aromatica"
antiMal$clean[antiMal$clean=="Uvariopsis congensis"]<-"Uvariopsis"
antiMal$isGenus<-!grepl(' ',antiMal$clean)
antiMal$id<-getId(antiMal$clean,taxaNamesAmbig)
antiMal$species<-getTaxonomy(antiMal$id,taxaNodes,taxaNames,'species')
antiMal$genus<-getTaxonomy(antiMal$id,taxaNodes,taxaNames,'genus')
write.csv(antiMal,'work/antiMal.csv',row.names=FALSE)

if(!exists('chimpFilter')){
  source('parseBlast.R')
  chimpFilter<-lapply(taxas,function(x)x[x[,'class']!='Mammalia'|is.na(x[,'class']),])
  codes<-sub('(rbcL|matK)(-bd)?_S[0-9]+_L[0-9]+_R([0-9]).*$','',sub('blast_trim_','',names(chimpFilter)))
  primers<-sub('.*(rbcL|matK)(-bd)?_S[0-9]+_L[0-9]+_R([0-9]).*$','\\1\\3',sub('blast_trim_','',names(chimpFilter)))
  bds<-sub('.*(rbcL|matK)(-bd)?_S[0-9]+_L[0-9]+_R([0-9]).*$','\\2',sub('blast_trim_','',names(chimpFilter)))=='-bd'
  #dont need currently
  rm(taxas)
}

propAntiMal<-sapply(chimpFilter,function(xx)mean(xx$species %in% antiMal$species[!is.na(antiMal$species)] | xx$genus %in% antiMal$genes[antiMal$isGenus]))

sampleInfo<-samples[codes,]
sampleInfo$primer<-primers
sampleInfo$bd<-bds
#assuming chow
sampleInfo[is.na(sampleInfo$Code),'Code']<-'Chow'
sampleInfo$speciesMal<-ifelse(sampleInfo$Code=='Chow','Chow',paste(sampleInfo$Species,ifelse(sampleInfo$malaria,'Pos','Neg')))
sampleInfo$file<-names(chimpFilter)

cols<-rainbow.lab(length(unique(sampleInfo$speciesMal)),alpha=.7)
names(cols)<-unique(sampleInfo$speciesMal)

tapply(propAntiMal,list(sampleInfo$Species,sampleInfo$malaria,sampleInfo$primer),median)
tapply(propAntiMal,list(sampleInfo$area,sampleInfo$primer),mean)
pdf('out/allAntiMal.pdf',width=10)
for(ii in unique(sampleInfo$primer)){
  pTL<-wilcox.test(
    propAntiMal[grepl('TL',sampleInfo$area)&!is.na(sampleInfo$area)&sampleInfo$bonobo&sampleInfo$primer==ii&sampleInfo$bd],
    propAntiMal[!grepl('TL',sampleInfo$area)&!is.na(sampleInfo$area)&sampleInfo$bonobo&sampleInfo$primer==ii&sampleInfo$bd],
    alternative='less'
  )$p.value
  withAs(si=sampleInfo[sampleInfo$primer==ii&(sampleInfo$bd|sampleInfo$speciesMal=='Chow'),],vpPlot(factor(ifelse(is.na(si$area),'Chow',si$area),levels=c('Chow',sort(unique(sampleInfo$area)))),propAntiMal[sampleInfo$primer==ii&(sampleInfo$bd|sampleInfo$speciesMal=='Chow')],ylab='Proportion antimalarial',las=1,main=sprintf('%s p(TL>=other bonobos)=%0.3f',ii,pTL),bg=cols[si$speciesMal],pch=21,cex=1.5))
  legend('topright',names(cols),pch=21,pt.bg=cols,pt.cex=1.5)
}
dev.off()


antiMalPlants<-sapply(chimpFilter,function(xx){
  ifelse(xx$species %in% antiMal$species[!is.na(antiMal$species)],xx$species,
    ifelse(xx$genus %in% antiMal$genes[antiMal$isGenus],xx$genus, 'Other')
  )
})

allFound<-unique(unlist(lapply(antiMalPlants,unique)))
antiMalProps<-apply(do.call(cbind,lapply(antiMalPlants,function(xx)table(c(xx,allFound))-1)),2,function(xx)xx/sum(xx))
pdf('out/antiMal.pdf',width=11)
  for(ii in rownames(antiMalProps)[rownames(antiMalProps)!='Other']){
    par(mfrow=c(4,1))
    for(jj in unique(sampleInfo$primer)){
      thisDat<-antiMalProps[ii,sampleInfo$primer==jj]
      names(thisDat)<-sampleInfo$Code[sampleInfo$primer==jj]
      pMal<-withAs(si=sampleInfo[sampleInfo$primer==jj,],wilcox.test(thisDat[si$malaria&!is.na(si$malaria)],thisDat[si$malaria&!is.na(si$malaria)]))$p.value
      if(is.na(pMal)|is.nan(pMal))pMal<-1
      if(pMal<.05)message('Sig: ',ii,' ',jj)
      pTL<-withAs(si=sampleInfo[sampleInfo$primer==jj,],wilcox.test(thisDat[grepl('TL',si$area)&si$bonobo],thisDat[!grepl('TL',si$area)&si$bonobo],alternative='less'))$p.value
      if(is.na(pTL)|is.nan(pTL))pTL<-1
      if(pTL<.05/nrow(antiMalProps))message('Sig: ',ii,' ',jj)
      barInfo<-barplot(thisDat,las=2,main=sprintf('%s %s (malaria p=%0.3f, TL p=%0.3f)',ii,jj,pMal,pTL),col=cols[sampleInfo[sampleInfo$primer==jj,'speciesMal']])[,1]
      axis(1,barInfo,names(thisDat),las=2)
    }
  }
dev.off()



ssReads<-lapply(chimpFilter,function(x)unique(rownames(x)[!is.na(x$species)&x$species=='Strychnos spinosa']))
names(ssReads)<-names(taxonomy)

pairFiles<-withAs(si=sampleInfo[grepl('matK',sampleInfo$primer)&!is.na(sampleInfo$Species)&sampleInfo$bd,],tapply(si$file,si$Code,c))
lapply(pairFiles['PA1044'],function(files){
  taxa1<-taxonomy[[files[1]]]
  chimp1<-chimpFilter[[files[1]]]
  chimp2<-chimpFilter[[files[2]]]
  taxa2<-taxonomy[[files[2]]]
  #ssReads<-unique(taxa2[taxa2$species=='Strychnos spinosa','qName'])
  ssReads2<-ssReads[[files[2]]]
  taxa1[taxa1$qName %in% ssReads2[1],]
  chimp2[rownames(chimp2) %in% ssReads2[10],]
  browser()
})
