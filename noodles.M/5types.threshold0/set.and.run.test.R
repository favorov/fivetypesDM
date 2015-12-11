source('../prepare.data.R')
if(! 'noodles.M.methylation' %in% ls()) stop('data is not prepared!')

threshold<-0;
test.typenames<-typenames;

source('../calc.significance.R ')
