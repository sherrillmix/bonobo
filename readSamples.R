samples<-read.csv('SampleList_microbiome.csv',stringsAsFactors=FALSE)
samples$malaria<-!grepl('[Nn]eg',samples$Plasmodium)
samples$bonobo<-samples$Species=='Pan paniscus'
samples$area<-sub('[0-9]+$','',samples$Code)
samples$area[samples$Site=='TL']<-'TL-E'
samples$area[samples$Site=='TL (West)']<-'TL-W'
samples$isTL<-samples$area %in% c('TL-E','TL-W')
samples$lat<-NA
samples$lat[samples$isTL]<-as.numeric(gsub('[^0-9.]','',samples$Lat[samples$isTL]))
samples$tlNorth<-samples$isTL&samples$lat<2.55
samples$area2<-samples$area
samples$area2[samples$tlNorth]<-'TL-NE'
samples$area2<-sub('TL-','TL2-',samples$area2)
samples$plasmoPM<-ifelse(samples$malaria,'+','-')
samples$chimpBonobo<-ifelse(samples$bonobo,'Bonobo','Chimp')
rownames(samples)<-samples$Code

primerSeqs<-list('rbcl'=c('rbcL634F'='ATGCGTTGGAGAGACCGTTTC','rbcLbR'='TCGGTYAGAGCRGGCATRTGCCA'),'matk'=c('matK472F'='CCCRTYCATCTGGAAATCTTGGTTC','matK1248R'='GCTRTRATAATGAGAAAGATTTCTGC'))
