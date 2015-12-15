if (!suppressWarnings(require('scatterplot3d')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("scatterplot3d")
	library("scatterplot3d")
}

PCA3Dplot<-function(pca,col='blue',pdfname=paste0('PCA_3D ',gsub(':',' ',date()),'.pdf')){
	pdf(pdfname, width = 14, height = 14)
        layout(matrix(c(1,2,3,4), 2, 2))
        par(oma = c(.1,1,.1,.1))
				varsum<-sum(pca$sdev)
        pr1 = round(pca$sdev[1]/varsum*100,0)
        pr2 = round(pca$sdev[2]/varsum*100,0)
        pr3 = round(pca$sdev[3]/varsum*100,0)
        plot(pca$load[,1],pca$load[,3], col = col, xlab = paste0('PC1 (', pr1, '% explained variance)'),
                ylab = paste0('PC3 (', pr3, '% explained variance)'), pch = 16, cex = 2, cex.lab = 1.5, cex.axis = 1.5)
        plot(pca$load[,1],pca$load[,2], col = col, xlab = paste0('PC1 (', pr1, '% explained variance)'),
                ylab = paste0('PC2 (', pr2, '% explained variance)'), pch = 16, cex = 2, cex.lab = 1.5, cex.axis = 1.5)
        plot(pca$load[,2],pca$load[,3], col = col, xlab = paste0('PC2 (', pr2, '% explained variance)'),
                ylab = paste0('PC3 (', pr3, '% explained variance)'), pch = 16, cex = 2, cex.lab = 1.5, cex.axis = 1.5)

        scatterplot3d(pca$load[,1],pca$load[,2],pca$load[,3],color = col, pch=16,
                xlab = 'PC1',ylab = 'PC2',zlab = 'PC3', zlim = c(-0.4,0.4), angle = 30, cex.symbols = 2, cex.lab = 1.2, cex.axis = 1.2)
dev.off()
}

PCA2Dplot<-function(pca,col='blue',n=3,pdfname=paste0('PCA_2D ',gsub(':',' ',date()),'.pdf')){
	pdf(pdfname)
	varsum<-sum(pca$sdev)
	expvar = round(pca$sdev/varsum*100,0)
	for(i in 1:n)
		for(j in i:n)
		{
			if (i==j) next;
			plot(pca$load[,i],pca$load[,j], pch = 16,col=col,
				xlab = paste0('PC',i,'(', expvar[i], '% explained variance)'),
				ylab = paste0('PC',j,'(', expvar[j], '% explained variance)')
				)
		}
	dev.off()
}

load("significant.DM.Rda")
Scol=c('red','blue','green','magenta','orange')
names(Scol)<-unique(test.typenames)
if (dim(significant.DM.methylation)[0]<dim(significant.DM.methylation)[1])
	message('We omit PCA cause number of features is not enough') 
else
{
	pca<-princomp(significant.DM.methylation,center=FALSE)
	PCA2Dplot(pca,col=Scol[test.typenames],pdfname='PCA_62_samples_2D.pdf')
	PCA3Dplot(pca,col=Scol[test.typenames],pdfname='PCA_62_samples_3D.pdf')

	pcaind<-princomp(significant.DM.methylation.binarised,center=FALSE)
	PCA2Dplot(pca,col=Scol[test.typenames],pdfname='PCA_62_samples_bin_2D.pdf')
	PCA3Dplot(pca,col=Scol[test.typenames],pdfname='PCA_62_samples_bin_3D.pdf')
}
