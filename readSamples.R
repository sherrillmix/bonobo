samples<-read.csv('SampleList_microbiome.csv')
samples$malaria<-!grepl('[Nn]eg',samples$Plasmodium)
