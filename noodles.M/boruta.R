library(differential.coverage)
load("significant.DM.Rda")

boruta.selected.loaded<-FALSE
if(file.exists('boruta.selected.Rda'))
{
	loaded<-load('boruta.selected.Rda')
	if ('boruta.selected.methylation' %in% loaded) 
		if (class(boruta.selected.methylation)=='matrix')
			boruta.selected.loaded<-TRUE
}

if(!boruta.selected.loaded)
{
	if (!suppressWarnings(require('Boruta')))
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("Boruta")
		library("Boruta")
	}

	#boruta.result<-Boruta(t(significant.DM.methylation),as.factor(test.typenames))

	boruta.result.bin<-Boruta(t(significant.DM.methylation.binarised),as.factor(test.typenames))

	#features<-which(boruta.result$finalDecision=='Confirmed')
	features.bin<-which(boruta.result.bin$finalDecision=='Confirmed')

	#boruta.selected.probes<-significant.DM.probes[features]
	boruta.bin.selected.probes<-significant.DM.probes[features.bin]

	#boruta.selected.methylation<-significant.DM.methylation[features,]
	#rownames(boruta.selected.methylation)<-sub('-(.*)','',as.character(boruta.selected.probes))
	
	boruta.bin.selected.methylation<-significant.DM.methylation.binarised[features.bin,]
	rownames(boruta.bin.selected.methylation)<-sub('-(.*)','',as.character(boruta.bin.selected.probes))

	save(file='boruta.selected.Rda',list=c(
		#'features',
		features.bin,
		#'boruta.selected.probes',
		'boruta.bin.selected.probes',
		#'boruta.selected.methylation',
		'boruta.bin.selected.methylation',
		'Scol',
		'test.typenames'))
}

pdf('heatmap.boruta.filtered.pdf')
heatmap(boruta.selected.methylation,ColSideColors = Scol[test.typenames])
dev.off()

pdf('heatmap.boruta.bin.filtered.pdf')
heatmap(boruta.bin.selected.methylation,ColSideColors = Scol[test.typenames])
dev.off()

#sink('boruta.selected.probes.txt')
#print(as.data.frame(closest.gene.start.by.interval(noodles = boruta.selected.probes)))
#sink()
boruta.bin.selected.probes$'p-value'<-significant.p.values[features.bin]
sink('boruta.bin.selected.probes.txt')
print(as.data.frame(closest.gene.start.by.interval(noodles = boruta.bin.selected.probes)))
sink()

