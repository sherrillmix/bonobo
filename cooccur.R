tl<-read.table('tlSamples.csv',sep='\t',header=TRUE)
tl$b1<-grepl('B1',tl$plas)
tl$c2<-grepl('C2',tl$plas)
tl$animal<-sub('[^(]+\\(([^)]+)\\)','\\1',tl$id)
b1c2<-data.frame('B1'=tapply(tl$b1,tl$animal,any),'C2'=tapply(tl$c2,tl$animal,any))
cooccur<-table(ifelse(b1c2$B1,'B1','noB1'),ifelse(b1c2$C2,'C2','noC2'))
nBonobo<-63
cooccur['noB1','noC2']<-nBonobo-sum(cooccur)
fisher.test(cooccur)

