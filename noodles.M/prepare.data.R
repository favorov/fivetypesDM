if (!suppressWarnings(require('differential.coverage')))
{
	if (!suppressWarnings(require('devtools')))
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("devtools")
		library("devtools")
	}
	install_github('favorov/Differential.Coverage')
	#load_all('../../../../differential.coverage/')
	library('differential.coverage')
}

noodles.M.loaded<-FALSE
# we prepare: noodles and all the bedfile stuff and methylation - to 5types.meth.data.Rda

#they are already here... ok
if ('noodles.M.methylation' %in% ls()) 
	if ('dgCMatrix' %in% class(noodles.M.methylation))
		if ('noodles.M' %in% ls())
			if(class(noodles.M)=='GRanges')
				noodles.M.loaded<-TRUE

#we run it from one of the location that are prticular cases of the pioeline,
#but this part is common
if(!noodles.M.loaded && file.exists('../5types.meth.data.Rda'))
{
	loaded<-load('../5types.meth.data.Rda')
	if ('noodles.M.methylation' %in% loaded) 
		if ('dgCMatrix' %in% class(noodles.M.methylation))
			if ('noodles.M' %in% loaded)
				if(class(noodles.M)=='GRanges')
			noodles.M.loaded<-TRUE
}

if(!noodles.M.loaded)
{
	beddir<-'../../cancer_bed_files/'
	noodle.length<-1000
	chrs<-nucl.chromosomes.hg19()
	if(!'noodles.M' %in% ls()) 
		noodles.M<-prepare.covering.noodles(chrs,noodle.length)
	bedfiles<-dir(beddir)
	typenames<-substr(bedfiles,0,2)
	typenames<-replace(typenames,typenames=='M_','M')
	typenames<-replace(typenames,typenames=='N_','N')
	bed.ids<-sapply(strsplit(split = '_',bedfiles),function(s){paste0(s[1],'_',s[2])})
	if(length(unique(bed.ids))!=length(bed.ids))
		stop('Non-unique bed.ids. So what?')
	noodles.M.methylation<-count.coverage.of.noodles(noodles.M,paste0(beddir,bedfiles),bed.ids)
	save(file='../5types.meth.data.Rda',list=c('noodles.M','noodles.M.methylation','typenames','bed.ids','noodle.length'))
}
