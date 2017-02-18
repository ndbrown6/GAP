'GAPplot' <- function(tmp, sampleName, SNP_annot, Centr_CL, NETM, p_BAF, q_LRR)
{
	tt = sort(union(which(tmp[,GAP::.GapEnv$c_len]>15),which(tmp[,GAP::.GapEnv$c_chr]==21)))
	tt = sort(union(tt,which(trunc(tmp[,GAP::.GapEnv$c_chr])==24)))
	if (length(tt)>0) {
		tmp1 = tmp[tt,]
	}
	Chr_counts = 0
	Centr_counts = 0
	centr = tmp1[1,]
	for (k in 1:(dim(tmp1)[1]-1)) {
		if ((tmp1[k+1,GAP::.GapEnv$c_chr]>tmp1[k,GAP::.GapEnv$c_chr])&&(round((tmp1[k+1,GAP::.GapEnv$c_chr]-0.2),0)==tmp1[k,GAP::.GapEnv$c_chr])) {
			centr = rbind(centr,tmp1[k:(k+1),])
			Chr_counts = Chr_counts+(tmp1[k,GAP::.GapEnv$c_CN]+tmp1[k+1,GAP::.GapEnv$c_CN])/2
			Centr_counts = Centr_counts+tmp1[k,GAP::.GapEnv$c_CN]+tmp1[k+1,GAP::.GapEnv$c_CN]
		} else {
			tt = round((tmp1[k+1,GAP::.GapEnv$c_chr]-0.2),0)
			if (tmp1[k+1,GAP::.GapEnv$c_chr]>tmp1[k,GAP::.GapEnv$c_chr]) {
				if ((tt==13)||(tt==14)||(tt==15)||(tt==21)||(tt==22)) {
					centr = rbind(centr,tmp1[k+1,])
					Chr_counts = Chr_counts+tmp1[k+1,GAP::.GapEnv$c_CN]
					Centr_counts = Centr_counts+tmp1[k+1,GAP::.GapEnv$c_CN]
				}
			}
		}
	}
	centr = centr[-1,]
	centr = cbind(centr[,1:7],centr[,16:17])
	C_Ind = 1
	C_Chr = 2
	C_Pos = 4
	DNAi = 0
	tt = 0
	k = 1
	for (k in 1:dim(tmp)[1]) {
		DNAi = DNAi+(SNP_annot[tmp[k,GAP::.GapEnv$c_if],C_Pos]-SNP_annot[tmp[k,GAP::.GapEnv$c_is],C_Pos])/1000*tmp[k,GAP::.GapEnv$c_CN]
		tt = tt+(SNP_annot[tmp[k,GAP::.GapEnv$c_if],C_Pos]-SNP_annot[tmp[k,GAP::.GapEnv$c_is],C_Pos])/1000
	}
	DNAi = round(DNAi/tt/2,2)
	plot(tmp[,GAP::.GapEnv$c_baf],tmp[,GAP::.GapEnv$c_lrr],main="Genome Alteration Print\n",xlim=c(0,1),ylim=c(-1,1),type="n",xlab="B Allele Frequency",ylab="Log2 Ratio", cex.main=1.3)
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
			rect(min(NETM[tt1,2])-0.002,min(NETM[tt1,1]),max(NETM[tt1,2])+0.002,max(NETM[tt1,1]), col=RB[p],border=RB[p])
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
				text(0.54,Centr_CL[5,j],paste(as.character(Centr_CL[1,j]),"/",as.character(Centr_CL[2,j]),sep=""))
			} else {
				text(Centr_CL[4,j],Centr_CL[5,j],paste(as.character(Centr_CL[1,j]),"/",as.character(Centr_CL[2,j]),sep=""))
			}
		}
	}
	tt = which(tmp[,GAP::.GapEnv$c_len]<30)
	if (length(tt)>0) {
		tmp1 = tmp[-tt,]
	}
	tmp1 = shrinkReprTMP(tmp1)
	text(0,-0.8,paste("p_BAF = ",as.character(round(p_BAF,2)),sep=""),cex=0.9, pos=4)
	text(0,-0.9,paste("q_LRR = ",as.character(round(q_LRR,2)),sep=""),cex=0.9, pos=4)
	text(0,-1,paste("DNAind = ",as.character(DNAi)),cex=0.9, pos=4)
	text(0.24,-0.8,paste("Breaks = ",as.character(dim(tmp1)[[1]]-41)),cex=0.9, pos=4)
	text(0.24,-0.9,paste("Chrs = ",as.character(Chr_counts)),cex=0.9, pos=4)
	text(0.24,-1,paste("Centers = ",as.character(Centr_counts)),cex=0.9, pos=4)
}
