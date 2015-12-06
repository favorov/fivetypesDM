source('calc.fisher.R')
#it will load or calc everything
significant.bonf.corr<-which(noodles.M.p.values*length(noodles.M.p.values)<0.05)
noodles.M.significant.bonf.corr<-noodles.M[significant.bonf.corr]
noodles.M.methylation.significant.bonf.corr<-noodles.M.methylation[significant.bonf.corr,]
noodles.M.p.values.significant.bonf.corr<-noodles.M.p.values[significant.bonf.corr]
save(file='noodles.M.significant.bonf.corr.Rda',list=c('noodles.M.significant.bonf.corr','noodles.M.methylation.significant.bonf.corr','noodles.M.p.values.significant.bonf.corr'))
