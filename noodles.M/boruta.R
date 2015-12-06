source('significant.probes.list.R')
library(Boruta)

Scol=c('red','blue','green','magenta','orange')
names(Scol)<-unique(typenames)

boruta.result<-Boruta(t(as.matrix(noodles.M.methylation.significant.bonf.corr)),as.factor(typenames))

noodles.M.methylation.significant.bonf.corr.ind<-0+(as.matrix(noodles.M.methylation.significant.bonf.corr)>0)

boruta.result.ind<-Boruta(t(noodles.M.methylation.significant.bonf.corr.ind),as.factor(typenames))

features<-which(boruta.result$finalDecision=='Confirmed')
features.ind<-which(boruta.result.ind$finalDecision=='Confirmed')

noodles.M.boruta<-noodles.M.significant.bonf.corr[features]
noodles.M.boruta.ind<-noodles.M.significant.bonf.corr[features.ind]

noodles.M.methylation.boruta<-
	as.matrix(noodles.M.methylation.significant.bonf.corr[features,])
noodles.M.methylation.boruta.ind<-
	noodles.M.methylation.significant.bonf.corr.ind[features.ind,]

pdf('heatmap.boruta.bin.filtered.pdf')
heatmap(noodles.m.methylation.boruta,ColSideColors = Scol[typenames])
dev.off()

pdf('heatmap.boruta.filtered.pdf')
heatmap(noodles.M.methylation.boruta.ind,ColSideColors = Scol[typenames])
dev.off()

save(file='5types.boruta.Rda',list=c('noodles.M.boruta','noodles.M.boruta.ind','noodles.M.methylation.boruta','noodles.M.methylation.boruta.ind','Scol','typenames'))

