source('../prepare.data.R')
if(! 'noodles.M.methylation' %in% ls()) stop('data is not prepared!')

threshold<-1; #inclusive
test.typenames<-typenames;
test.typenames[test.typenames=='T8' | test.typenames=='T9']='T89'

Scol=c('red','blue','green','magenta','orange')
names(Scol)<-unique(test.typenames)

source('../calc.significance.R')
source('../PCA.R')
source('../boruta.R')
