'GAPplot_preliminary' <- function(tmp, Centr_CL, p_BAF, q_LRR, NETM, sampleName="")
{
	plot(tmp[,GAP::.GapEnv$c_baf],tmp[,GAP::.GapEnv$c_lrr], main="Genome Alteration Print\n", xlim=c(0,1), ylim=c(-1,1), type="n", xlab="B Allele Frequency", ylab="Log2 Ratio", cex.main=1.3)
	title(main=paste("\n", sampleName, sep=""), cex.main=1.1)
	CEX = tmp[,GAP::.GapEnv$c_len]/quantile(tmp[,GAP::.GapEnv$c_len],prob=0.7)
	tt = which(CEX>4)
	if (length(tt)>0) {
		CEX[tt] = 4
	}
	tt = which(CEX<0.03)
	if (length(tt)>0) {
		CEX[tt] = 0
	}
	points(1-tmp[,GAP::.GapEnv$c_baf],tmp[,GAP::.GapEnv$c_lrr],cex=CEX,col="black")
	points(tmp[,GAP::.GapEnv$c_baf],tmp[,GAP::.GapEnv$c_lrr],cex=CEX,col="grey")
	abline(h=0)
	abline(v=0.5)
	RB = rainbow(40)
 	for (p in 1:25) {
		tt1 = which(NETM[,3]==p)
		if (length(tt1)>0) {
			rect(min(NETM[tt1,2])-.002,min(NETM[tt1,1]),max(NETM[tt1,2])+.002,max(NETM[tt1,1]), col=RB[p], border=RB[p])
		}
	}
	for (p in 31:38) {
		tt1 = which(NETM[,3]==p)
		if (length(tt1)>0) {
			rect(min(NETM[tt1,2]),min(NETM[tt1,1]),max(NETM[tt1,2]),max(NETM[tt1,1]), col=RB[p],border=RB[p])
		}
	}
	for (j in 1:dim(Centr_CL)[2]) {
		if (Centr_CL[1,j]<8) {
			if (Centr_CL[1,j]==Centr_CL[2,j]*2) {
				text(0.54,Centr_CL[5,j],paste(as.character(Centr_CL[1,j]),"/",as.character(Centr_CL[2,j]),sep=""), cex=.8)
			} else {
				text(Centr_CL[4,j],Centr_CL[5,j],paste(as.character(Centr_CL[1,j]),"/",as.character(Centr_CL[2,j]),sep=""), cex=.8)
			}
		}
	}
	text(0,-0.8,paste("p_BAF = ",as.character(round(p_BAF,2)),sep=""),cex=0.9, pos=4)
	text(0,-0.9,paste("q_LRR = ",as.character(round(q_LRR,2)),sep=""),cex=0.9, pos=4)
	text(0,-1,paste("2 Copy LRR = ",as.character(Centr_CL[5,4])),cex=0.9, pos=4)
}
