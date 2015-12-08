source('../prepare.data.R')
if(! 'noodles.M.methylation' %in% ls()) stop('data is not prepared!')

noodles.M.p.values.loaded<-FALSE
#they are already here... ok
if ('noodles.M.p.values' %in% ls()) 
	if (class(noodles.M.p.values)=='numeric')
		if (length(noodles.M)==length(noodles.M.p.values))
			noodles.M.p.values.loaded<-TRUE

if(!noodles.M.p.values.loaded && file.exists('5types.p.values.Rda'))
{
	loaded<-load('5types.p.values.Rda')
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
			t<-table(typenames,as.logical(met))
			if(dim(t)[2]==2)
				fisher.test(t)$p.value
			else
				1.
		}
	)
	save(file='p.values.Rda',list=c('noodles.M.p.values','typenames'))
}
