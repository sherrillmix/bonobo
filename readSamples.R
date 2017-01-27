samples<-read.csv('SampleList_microbiome.csv',stringsAsFactors=FALSE)
samples$malaria<-!grepl('[Nn]eg',samples$Plasmodium)
samples$bonobo<-samples$Species=='Pan paniscus'
samples$area<-sub('[0-9]+$','',samples$Code)
samples$area[samples$Site=='TL']<-'TL-E'
samples$area[samples$Site=='TL (West)']<-'TL-W'
