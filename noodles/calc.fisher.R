if(! 'noodles.M.methylation' %in% ls()) load('5types.meth.data.Rda')
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
save(file='5types.p.values.Rda',list=c('noodles.M.p.values'))

