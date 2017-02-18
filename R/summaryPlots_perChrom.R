'summaryPlots_perChrom' <- function(tmp, sampleName="", SNP_annot, datA, datL, Centr_CL)
{
	C_Ind = 1
	C_Chr = 2
	C_Pos = 4
	for (k in 1:22) {
		pdf(file=paste(sampleName,"_Chr_",as.character(k),".pdf",sep=""))
		KN = which(SNP_annot[,C_Chr]==k)
		KK = KN[length(KN)]
		KN = KN[1]
		plot(c(1,SNP_annot[KK,C_Pos]/1000),c(0,0), ylim=c(-15,15), main="", type="n", xlab="Position (Mb)", ylab="", xaxt="n", yaxt="n")
		title(main=paste(sampleName, "\n", sep=""), cex.main=1.1)
		title(main=paste("\nChromosome ", as.character(k), sep=""), cex.main=0.9)
		labs = pretty(seq(from=1,to=SNP_annot[KK,C_Pos]/1000, length=5))
		axis(side=1, at=labs, labels=ceiling(labs/1000))
		points(SNP_annot[KN:KK,C_Pos]/1000,datL[KN:KK]*4-9, col="lightgrey", pch=".")
		KN = which(SNP_annot[,C_Chr]==k)
		KK = KN[length(KN)]
		KN = KN[1]
		points(SNP_annot[KN:KK,C_Pos]/1000, 6*(0.5+(datA[KN:KK]-0.5)*sample(c(-1,1),(KK-KN+1),replace=TRUE))+7, col="lightgrey", pch=".")
		abline(h=(-2),col="grey20", lty=3)
		abline(h=0,col="grey20", lty=3)
		axis(side=2, at=c(0,-2), labels=c("4","2"), tick=FALSE, las=1)
		KN = which(trunc(tmp[,GAP::.GapEnv$c_chr])==k)
		KK = KN[length(KN)]
		KN = KN[1]
		p = KN
		for (p in KN:KK) {
			if (tmp[p,GAP::.GapEnv$c_len]>20) {
				lines(c(SNP_annot[tmp[p,GAP::.GapEnv$c_is],C_Pos]/1000,SNP_annot[tmp[p,GAP::.GapEnv$c_if],C_Pos]/1000),c(tmp[p,GAP::.GapEnv$c_lrr]*4-9,tmp[p,GAP::.GapEnv$c_lrr]*4-9),col="blue",lwd=3)
			}
		}
		KN = which(trunc(tmp[,GAP::.GapEnv$c_chr])==k)
		KK = KN[length(KN)]
		KN = KN[1]
		p = KN
		for (p in KN:KK) {
			lines(c(SNP_annot[tmp[p,GAP::.GapEnv$c_is],C_Pos]/1000,SNP_annot[tmp[p,GAP::.GapEnv$c_if],C_Pos]/1000),c(tmp[p,GAP::.GapEnv$c_lrr]*4-9,tmp[p,GAP::.GapEnv$c_lrr]*4-9),col="red",type="l",lwd=3)
			lines(c(SNP_annot[tmp[p,GAP::.GapEnv$c_is],C_Pos]/1000,SNP_annot[tmp[p,GAP::.GapEnv$c_if],C_Pos]/1000),c(tmp[p,GAP::.GapEnv$c_baf],tmp[p,GAP::.GapEnv$c_baf])*6+7,col="red",lwd=3)
			lines(c(SNP_annot[tmp[p,GAP::.GapEnv$c_is],C_Pos]/1000,SNP_annot[tmp[p,GAP::.GapEnv$c_if],C_Pos]/1000),c(1-tmp[p,GAP::.GapEnv$c_baf],1-tmp[p,GAP::.GapEnv$c_baf])*6+7,col="red",lwd=3)
			lines(c(SNP_annot[tmp[p,GAP::.GapEnv$c_is],C_Pos]/1000,SNP_annot[tmp[p,GAP::.GapEnv$c_if],C_Pos]/1000),c(tmp[p,GAP::.GapEnv$c_CN],tmp[p,GAP::.GapEnv$c_CN])-4,col="blue",lwd=3)
		}
		if (KK>(KN+1)) {
			for (p in (KN+1):KK) {
				lines(c(SNP_annot[tmp[p-1,GAP::.GapEnv$c_if],C_Pos]/1000,SNP_annot[tmp[p,GAP::.GapEnv$c_is],C_Pos]/1000),c(tmp[p-1,GAP::.GapEnv$c_CN],tmp[p,GAP::.GapEnv$c_CN])-4,col="blue",lwd=1)
			}
		}
		for (p in KN:KK) {
			if (tmp[p,GAP::.GapEnv$c_CN]==tmp[p,GAP::.GapEnv$c_BA]) {
				lines(c(SNP_annot[tmp[p,GAP::.GapEnv$c_is],C_Pos]/1000,SNP_annot[tmp[p,GAP::.GapEnv$c_if],C_Pos]/1000),c(tmp[p,GAP::.GapEnv$c_CN],tmp[p,GAP::.GapEnv$c_CN])-4,col="red",lwd=3)
			}
		}
		for (p in KN:KK) {
			if (tmp[p,GAP::.GapEnv$c_BA]==0) {
				lines(c(SNP_annot[tmp[p,GAP::.GapEnv$c_is],C_Pos]/1000,SNP_annot[tmp[p,GAP::.GapEnv$c_if],C_Pos]/1000),c(tmp[p,GAP::.GapEnv$c_CN],tmp[p,GAP::.GapEnv$c_CN])-4,col="green",lwd=3)
			}
		}
		for (p in KN:KK) {
			bb = ""
			if ((tmp[p,GAP::.GapEnv$c_BA]==0)||(tmp[p,GAP::.GapEnv$c_BA]==tmp[p,GAP::.GapEnv$c_CN])) {
				if (tmp[p,GAP::.GapEnv$c_CN]>1) {
					for (s in 1:tmp[p,GAP::.GapEnv$c_CN]) {
						bb = paste(bb,"A",sep="")
					}
				} else {
					if (tmp[p,GAP::.GapEnv$c_CN]==1) {
						bb = "A"
					} else {
						bb = "0"
					}
				}
			} else {
				if (tmp[p,GAP::.GapEnv$c_BA]>1) {
					for (s in 1:tmp[p,GAP::.GapEnv$c_BA]) {
						bb = paste(bb,"A",sep="")
					}
				} else {
					if (tmp[p,GAP::.GapEnv$c_BA]==1) {
						bb = "A"
					}
				}
				if (tmp[p,GAP::.GapEnv$c_CN]>tmp[p,GAP::.GapEnv$c_BA]+1) {
					for (s in (tmp[p,GAP::.GapEnv$c_BA]+1):tmp[p,GAP::.GapEnv$c_CN]) {
						bb = paste(bb,"B",sep="")
					}
				} else {
					if (tmp[p,GAP::.GapEnv$c_CN]==tmp[p,GAP::.GapEnv$c_BA]+1) {
						bb = paste(bb,"B",sep="")
					}
				}
			}
			if ((tmp[p,GAP::.GapEnv$c_if]-tmp[p,GAP::.GapEnv$c_is])/2>200) {
				if (tmp[p,GAP::.GapEnv$c_is]==1) {
					text((SNP_annot[tmp[p,GAP::.GapEnv$c_is],C_Pos]+SNP_annot[tmp[p,GAP::.GapEnv$c_if],C_Pos])/2000,tmp[p,GAP::.GapEnv$c_CN]-4.5,bb,cex=0.3, col="grey30")
				} else {
					text((SNP_annot[tmp[p,GAP::.GapEnv$c_is],C_Pos]+SNP_annot[tmp[p,GAP::.GapEnv$c_if],C_Pos])/2000,tmp[p,GAP::.GapEnv$c_CN]-3.3,bb,cex=0.3, col="grey30")
				}
			}
		}
		dev.off()
	}
}
