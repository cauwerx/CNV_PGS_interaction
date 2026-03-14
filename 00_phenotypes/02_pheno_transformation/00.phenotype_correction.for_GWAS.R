source('../software/genABEL.R')

FILE='./project/data/phenotypes/pheno_continuous_train_raw.tsv.gz'

d<-as.data.frame(data.table::fread(FILE, hea=T, sep='\t'))

dd<-data.frame(FID=d$FID, IID=d$IID)

for (i in 3:length(colnames(d))){
        PHN=colnames(d)[i]
        sub<-d[colnames(d)==PHN]
        dd$pheno<-rntransform(sub[,1])
        colnames(dd)[i]<-PHN
}

write.table(dd, './project/data/phenotypes/pheno_continuous_train_raw.irnt.tsv', quote=F, col.names=T, row.names=F, sep='\t')

