'plotGAPdensity' <- function(NET, contourPlot=FALSE, Centr_CL=NULL)
{
	if (!contourPlot) {
		NR = dim(NET)[1]
		NF = dim(NET)[2]
		SR = as.numeric(rownames(NET))
		SF = as.numeric(colnames(NET))
		plot(c(SF[1],SF[NF]), c(max(-1,SR[1]),SR[NR]), type="n", main="Genome Alteration Print\n", xlab="Major B allele Frequency", ylab="Log2 Ratio", xlim=c(0.5,1), cex.main=1.3)
		title(main="\nDensity Estimate", cex.main=1.1)
		NETM = as.numeric(NET)
		NETM = cbind(NETM,NETM,NETM)
		NETM[,1] = rep(seq(1,length(SR)),length(SF))
		NETM[,2] = rep(seq(1,length(SF)),each=length(SR))
		NETM = NETM[-which(NETM[,3]<10),]
		NETM[,3] = log2(NETM[,3])
		NETM = cbind(NETM,NETM[,3])
		NETM = NETM[order(NETM[,3],decreasing=FALSE),]
		hh = hist(NETM[,3], breaks=64, plot=FALSE)$breaks
		p = 1
		k = 1
		while (p<=length(hh)) {
			while ((k<dim(NETM)[[1]])&&(NETM[k,3]<hh[p])) {
				NETM[k,4] = p
				k = k+1
			}
			p = p+1
		}
		tt = which(NETM[1,]==NR)
		if (length(tt)>0) {
			NETM[1,tt] = NR-1
		}
		tt = which(NETM[2,]==NF)
		if (length(tt)>0) {
			NETM[2,tt] = NF-1
		}
		RB = rainbow(length(hh))
		for (p in 1:dim(NETM)[1]) {
			rect(SF[NETM[p,2]], SR[NETM[p,1]], SF[NETM[p,2]+1], SR[NETM[p,1]+1], col=RB[NETM[p,4]], border=RB[NETM[p,4]])
		}
	} else {
		NR = dim(NET)[1]
		NF = dim(NET)[2]
		SR = as.numeric(rownames(NET))
		SF = as.numeric(colnames(NET))
		plot(c(SF[1],SF[NF]), c(max(-1,SR[1]),SR[NR]), type="n", main="Genome Alteration Print\n", xlab="Major B allele Frequency", ylab="Log2 Ratio", xlim=c(0.5,1), cex.main=1.3)
		title(main="\nDensity Estimate", cex.main=1.1)
		contour(x=as.numeric(colnames(NET)), y=as.numeric(rownames(NET)), z=t(NET), nlevels=100, col=rainbow(100), add=TRUE, drawlabels=FALSE)
		for (j in 1:dim(Centr_CL)[2]) {
			text(Centr_CL[4,j],Centr_CL[5,j],paste(as.character(Centr_CL[1,j]),"/",as.character(Centr_CL[2,j]),sep=""),cex=.9, col="grey40")
		}
	}
}
