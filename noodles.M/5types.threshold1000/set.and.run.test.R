source('../prepare.data.R')
if(! 'noodles.M.methylation' %in% ls()) stop('data is not prepared!')

threshold<-1000;
test.typenames<-typenames;

source('../calc.significance.R')
source('../PCA.R ')
source('../boruta.R ')
