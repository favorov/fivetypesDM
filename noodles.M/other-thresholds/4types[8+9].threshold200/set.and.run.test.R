source('../prepare.data.R')
if(! 'noodles.M.methylation' %in% ls()) stop('data is not prepared!')

threshold<-200; #inclusive
test.typenames<-typenames;
test.typenames[test.typenames=='T8' | test.typenames=='T9']='T89'

source('../calc.significance.R')
source('../PCA.R')
source('../boruta.R')
