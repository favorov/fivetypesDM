if (!'threshold' %in% ls()) stop('I am calc.singicance. I need threshold to run')

if (!'test.typenames' %in% ls()) stop('I am calc.singicance. I need test.typenames to run')

significant.DM.loaded<-FALSE

if(file.exists('significant.DM.Rda'))
{
	loaded<-load('significant.DM.Rda')
	if ('significant.p.values' %in% loaded) 
		if (class('significant.p.values')=='numeric')
			if (length(noodles.M)==length(significant.p.values))
				significant.DM.loaded<-TRUE
}


if(!significant.DM.loaded)
{
	noodles.M.p.values.loaded<-FALSE

	if(file.exists('p.values.Rda'))
	{
		loaded<-load('p.values.Rda')
		if ('noodles.M.p.values' %in% loaded) 
			if (class(noodles.M.p.values)=='numeric')
				if (length(noodles.M)==length(noodles.M.p.values))
					noodles.M.p.values.loaded<-TRUE
	}

	if(!noodles.M.p.values.loaded)
	{
		noodles.M.p.values<-apply(
			noodles.M.methylation,1,
			function(met) 
			{
				t<-table(test.typenames,met>=threshold)
				if(dim(t)[2]==2)
					fisher.test(t)$p.value
				else
					1.
			}
		)
		save(file='p.values.Rda',list=c('noodles.M.p.values','test.typenames'))
	}
	which.probes.are.significant<-which(noodles.M.p.values*length(noodles.M.p.values)<0.05) #bonferroni correction
	significant.DM.probes<-noodles.M[which.probes.are.significant]
	significant.DM.methylation<-as.matrix(noodles.M.methylation[which.probes.are.significant,])
	significant.DM.methylation.binarised<-0+(as.matrix(significant.DM.methylation)>=threshold) # 0/1
	significant.p.values<-noodles.M.p.values[which.probes.are.significant]
	save(file='significant.DM.Rda',list=c('significant.DM.probes','significant.DM.methylation','significant.DM.methylation.binarised','significant.p.values','test.typenames'))
}
