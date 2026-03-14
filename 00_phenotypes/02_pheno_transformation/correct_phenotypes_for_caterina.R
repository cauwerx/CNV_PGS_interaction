


source('src/genABEL.R')

FILE='data/pheno_continuous_train_raw.v2.tsv.gz'
d<-as.data.frame(data.table::fread(FILE, hea=T, sep='\t'))
#n<-unlist(lapply(colnames(d), FUN=function(x){unlist(strsplit(x,'-'))[1]}))
#colnames(d)<-n

dd<-data.frame(FID=d$FID, IID=d$IID)




for (i in 3:length(colnames(d))){
	PHN=colnames(d)[i]
	sub<-d[colnames(d)==PHN]
#	p<-apply(sub, 1, FUN=function(x){mean(x, na.rm=T)})
	p<-rntransform(sub[,1])
	dd$pheno<-rntransform(p)
	colnames(dd)[i]<-PHN

}

write.table(dd, 'data/pheno_continuous_train_raw.v2.irnt.tsv', quote=F, col.names=T, row.names=F, sep='\t')
