if(! 'noodles.M.methylation' %in% ls()) load('5types.meth.data.Rda')
if(! 'noodles.M.p.values'%in% ls()) load('5types.p.values.Rda')
significant.bonf.corr<-which(noodles.M.p.values*length(noodles.M.p.values)<0.05)
noodles.M.significant.bonf.corr<-noodles.M[significant.bonf.corr]
noodles.M.methylation.significant.bonf.corr<-noodles.M.methylation[significant.bonf.corr]
save(file='noodles.M.significant.bonf.corr.Rda',list=c('noodles.M.significant.bonf.corr','noodles.M.methylation.significant.bonf.corr'))
