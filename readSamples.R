library(lubridate)
samples<-read.csv('SampleList_microbiome.csv',stringsAsFactors=FALSE)
samples$malaria<-!grepl('[Nn]eg',samples$Plasmodium)
samples$bonobo<-samples$Species=='Pan paniscus'
samples$area<-sub('[0-9]+$','',samples$Code)
samples$area[samples$Site=='TL']<-'TL-E'
samples$area[samples$Site=='TL (West)']<-'TL-W'
samples$isTL<-samples$area %in% c('TL-E','TL-W')
samples$rDate<-mdy(samples$Date)
samples$month<-month(samples$rDate)
samples$lat<-NA
samples$lat[samples$isTL]<-as.numeric(gsub('[^0-9.]','',samples$Lat[samples$isTL]))
samples$tlNorth<-samples$isTL&samples$lat<2.55
samples$area2<-samples$area
samples$area2[samples$tlNorth]<-'TL-NE'
samples$area2<-sub('TL-','TL2-',samples$area2)
samples$plasmoPM<-ifelse(samples$malaria,'+','-')
samples$chimpBonobo<-ifelse(samples$bonobo,'Bonobo','Chimp')
#arbitrary seasons
samples$season<-as.character(ceiling(samples$month/3))
rownames(samples)<-samples$Code
#only infected with p vivax so not clearly negative or positive
samples<-samples[samples$Code!='UB2041',]

#clean up lat lon
samples$cleanLat<-trimws(sub('^0','',sub('([NS]) *([0-9.]+).*','\\2 \\1',samples$Lat)))
#all other IK are S so assuming missing S
samples[samples$area=='IK'&!grepl('S',samples$cleanLat),'cleanLat']<-sprintf('%sS',samples[samples$area=='IK'&!grepl('S',samples$cleanLat),'cleanLat'])
probs<-samples$cleanLat[!grepl('[0-9.]+ [NS]',samples$cleanLat)&samples$cleanLat!='']
probHemisphere<-substring(probs,nchar(probs))
convertDeg<-function(xx){
  if(length(xx)==1)return(xx)
  if(length(xx)==2)return(as.numeric(xx[1])+as.numeric(xx[2])/60)
  if(length(xx)==3)return(as.numeric(xx[1])+as.numeric(xx[2])/60+as.numeric(xx[2])/60/60)
}
samples[!grepl('[0-9.]+ [NS]',samples$cleanLat)&samples$cleanLat!='','cleanLat']<-sprintf('%0.5f %s',sapply(strsplit(probs,'[^0-9.]+'),convertDeg),probHemisphere)
samples$cleanLon<-sub(';','.',sub('^0','',sub('([EW]) *([0-9.]+).*','\\2 \\1',trimws(samples$Lon)))) 
#all sites are eastern hemisphere
samples[!grepl('E',samples$cleanLon)&samples$cleanLon!='','cleanLon']<-sprintf('%sE',samples[!grepl('E',samples$cleanLon)&samples$cleanLon!='','cleanLon'])
probs<-samples$cleanLon[!grepl('[0-9.]+ [EW]',samples$cleanLon)&samples$cleanLon!='']
probHemisphere<-substring(probs,nchar(probs))
samples[!grepl('[0-9.]+ [EW]',samples$cleanLon)&samples$cleanLon!='','cleanLon']<-sprintf('%0.5f %s',sapply(strsplit(probs,'[^0-9.]+'),convertDeg),probHemisphere)


primerSeqs<-list('rbcl'=c('rbcL634F'='ATGCGTTGGAGAGACCGTTTC','rbcLbR'='TCGGTYAGAGCRGGCATRTGCCA'),'matk'=c('matK472F'='CCCRTYCATCTGGAAATCTTGGTTC','matK1248R'='GCTRTRATAATGAGAAAGATTTCTGC'))
