library(differential.coverage)
load("significant.DM.Rda")

boruta.selected.loaded<-FALSE
if(file.exists('boruta.selected.Rda'))
{
	loaded<-load('boruta.selected.Rda')
	if ('boruta.bin.selected.methylation' %in% loaded) 
		if (class(boruta.bin.selected.methylation)=='matrix')
			boruta.selected.loaded<-TRUE
}

if(!boruta.selected.loaded)
{
	if (dim(significant.DM.methylation.binarised)[2] < 2)
	{
		features.bin<-c(1)
	} else {
		if (!suppressWarnings(require('Boruta')))
		{
			source("http://bioconductor.org/biocLite.R")
			biocLite("Boruta")
			library("Boruta")
		}

		#boruta.result<-Boruta(t(significant.DM.methylation),as.factor(test.typenames))
		set.seed("1248312")
		
		boruta.result.bin<-Boruta(t(significant.DM.methylation.binarised),as.factor(test.typenames),maxRuns = 1000)

		#features<-which(boruta.result$finalDecision=='Confirmed')
		features.bin<-which(boruta.result.bin$finalDecision=='Confirmed')
	}
	#boruta.selected.probes<-significant.DM.probes[features]
	boruta.bin.selected.probes<-significant.DM.probes[features.bin]

	#boruta.selected.methylation<-significant.DM.methylation[features,]
	#rownames(boruta.selected.methylation)<-sub('-(.*)','',as.character(boruta.selected.probes))
	
	boruta.bin.selected.methylation<-significant.DM.methylation.binarised[features.bin,]
	if(length(features.bin)>1) rownames(boruta.bin.selected.methylation)<-sub('-(.*)','',as.character(boruta.bin.selected.probes))

	save(file='boruta.selected.Rda',list=c(
		#'features',
		'features.bin',
		#'boruta.selected.probes',
		'boruta.bin.selected.probes',
		#'boruta.selected.methylation',
		'boruta.bin.selected.methylation',
		'Scol',
		'test.typenames'))
}

#pdf('heatmap.boruta.filtered.pdf')
#heatmap(boruta.selected.methylation,ColSideColors = Scol[test.typenames])
#dev.off()

if(min(boruta.bin.selected.methylation)>1){
	pdf('heatmap.boruta.bin.confirmed.pdf')
	heatmap(boruta.bin.selected.methylation,ColSideColors = Scol[test.typenames])
	dev.off()
}

boruta.annotated.loaded<-FALSE
if(file.exists('boruta.annotated.Rda'))
{
	loaded<-load('boruta.annotated.Rda')
	if ('boruta.bin.selected.annotated.probes' %in% loaded) 
		if (class(boruta.bin.selected.annotated.probes)=='data.frame')
			boruta.annotated.loaded<-TRUE
}

if(!boruta.annotated.loaded)
{
	#boruta.selected.probes is a GRanges; we start to annotate it
	#first we add Fisher test p-values
	boruta.bin.selected.probes$'p.value'<-significant.p.values[features.bin]
	#now, we add closest genes (do not forget, hg 18!)
	boruta.bin.selected.probes<-
		closest.gene.start.by.interval(noodles = boruta.bin.selected.probes,genome.id='hg18')
	#now, we need a data frame
	boruta.bin.selected.annotated.probes<-as.data.frame(boruta.bin.selected.probes)
	#let's create a noddles x samples table with the highes MACS peak scores
	# for each noodle and sample
	peaks<-max.peak.score.for.each.noodle(noodles = boruta.bin.selected.probes,paste0(beddir,bedfiles),bed.ids)
	# and now, we aggregate min/med/max for all the samples of the a type for each noodle
	#types are original 5 types
	score.annotation<-t(apply(peaks,1,function(a) unlist(tapply(a,typenames,function(a) c(min(a),median(a),max(a))))))
	#and add pretty colnames , min.M, med.T8, etc
	colnames(score.annotation)<-do.call('paste',c(expand.grid(c('min','med','max'),unique(typenames)),sep='.'))
	#and we bind it to the 	boruta.bin.selected.annotated.probes
	boruta.bin.selected.annotated.probes<-cbind(boruta.bin.selected.annotated.probes,score.annotation)
	#but, we can have some other that typenames types (like M+T6, etc)
	#they are coded in test.typenames
	if(!min(test.typenames==typenames))
	{
		score.annotation<-data.frame(t(apply(peaks,1,function(a) unlist(tapply(a,test.typenames,function(a) c(min(a),median(a),max(a)))))))
		#and add pretty colnames , min.M, med.T8, etc
		colnames(score.annotation)<-do.call('paste',c(expand.grid(c('min','med','max'),unique(test.typenames)),sep='.'))
		#we need only those columns that are not already in boruta.bin.selected.annotated.probes  
		score.annotation<-score.annotation[,!(colnames(score.annotation) %in% colnames(boruta.bin.selected.annotated.probes))]
		#and we bind it to the 	boruta.bin.selected.annotated.probes
		boruta.bin.selected.annotated.probes<-cbind(boruta.bin.selected.annotated.probes,score.annotation)
	}

	save(file='boruta.annotated.Rda',list=c('boruta.bin.selected.annotated.probes'))
}

#sink('boruta.selected.probes.txt')
#print(as.data.frame(closest.gene.start.by.interval(noodles = boruta.selected.probes)))
#sink()

options(width = 500) 
sink('boruta.bin.selected.probes.txt')
write.table(boruta.bin.selected.annotated.probes,quote = FALSE,row.names = FALSE,sep='\t')
sink()

