'summaryPlots' <- function(tmp, sampleName, SNP_annot, datA, datL, Centr_CL, NETM, p_BAF, q_LRR)
{
	C_Ind = 1
	C_Chr = 2
	C_Pos = 4
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
				if ((tt==13)||(tt==14)||(tt==15)||(tt==22)) {
					centr = rbind(centr,tmp1[k+1,])
					Chr_counts = Chr_counts+tmp1[k+1,GAP::.GapEnv$c_CN]
					Centr_counts = Centr_counts+tmp1[k+1,GAP::.GapEnv$c_CN]
				}
			}
		}
	}
	centr = centr[-1,]
	tt = table(SNP_annot[,C_Chr],SNP_annot[,C_Chr+1])
	CHR = c(0)
	tCHR = c(0)
	for (k in 1:dim(tt)[1]) {
		CHR = c(CHR,CHR[2*k-1]+tt[k,1])
		CHR = c(CHR,CHR[2*k]+tt[k,2])
		tCHR = c(tCHR,as.numeric(rownames(tt))[k])
		tCHR = c(tCHR,as.numeric(rownames(tt))[k]+0.5)
	}
	names(CHR) = tCHR
	rm(tCHR)
	k = 2
	while (k<=(length(CHR)-1)) {
		if ((CHR[k]-CHR[k-1]<3)||(CHR[k+1]-CHR[k]<3)) {
			CHR[k] = CHR[k+1]
			names(CHR)[k+1] = 0
		}
		k = k+2
	}
	if (length(which(names(CHR)=="0"))) {
		CHR = CHR[-which(names(CHR)=="0")]
	}
	CHR = c(0,CHR)
	names(CHR)[1] = "0"
	split.screen(c(2,1)) 
	screen(1)
	plot(c(-10000, max(tmp[,GAP::.GapEnv$c_if])+2000), c(-3,3), type="n", main=sampleName, xlab="", ylab="" , lab=c(7,1,7), frame.plot=FALSE, yaxt="n", xaxt="n")
	for (k in 2:(length(CHR)-2)) {
		if (round((as.numeric(names(CHR))[k]),0)==as.numeric(names(CHR))[k]) {
			text((CHR[k-1]+CHR[k])/2,3.1,names(CHR)[k],col="black",cex=.8)
			text((CHR[k-1]+CHR[k])/2,-2.6,names(CHR)[k],col="black",cex=.8)
			
		}
	}
	RB = rainbow(40)
	for (p in 1:dim(tmp)[1]) {
		cl = which(abs(Centr_CL[1,]-tmp[p,GAP::.GapEnv$c_CN])==0)
		ba = which(abs(Centr_CL[2,cl]-tmp[p,GAP::.GapEnv$c_BA])==0)
		if (length(ba)>0) {
			tmp[p,GAP::.GapEnv$c_sb] = Centr_CL[3,cl[ba[1]]]
		} else {
			tmp[p,GAP::.GapEnv$c_sb] = 30+tmp[p,GAP::.GapEnv$c_CN]
		}
	}
	tt = which(SNP_annot[,C_Ind]%%10==0)
	Cols = rep(rainbow_hcl(23, start=30, end=300)[c(2,17)], 12)
	points(SNP_annot[tt,C_Ind],datL[tt]/2-1.5,col=Cols[SNP_annot[tt,2]],pch=3, cex=.1)
	tt = which(datA<0.97)
	tt = tt[which(tt%%2==0)]
	points(SNP_annot[tt,C_Ind], (0.5+(datA[tt]-0.5)*sample(c(-1,1),(length(datA[tt])),replace=TRUE))+1.5, col=Cols[SNP_annot[tt,2]], pch=3, cex=.1)
	lines(c(tmp[1,GAP::.GapEnv$c_is],tmp[dim(tmp)[[1]],GAP::.GapEnv$c_if]),c(1.5,1.5), col="black",type="l",lwd=2)
	lines(c(tmp[1,GAP::.GapEnv$c_is],tmp[dim(tmp)[[1]],GAP::.GapEnv$c_if]),c(2.5,2.5), col="black",type="l",lwd=2)
	for (p in 1:dim(tmp)[1]) {
		rect(tmp[p,GAP::.GapEnv$c_is],-2.8,tmp[p,GAP::.GapEnv$c_if],-3, col=RB[tmp[p,GAP::.GapEnv$c_sb]],border=RB[tmp[p,GAP::.GapEnv$c_sb]])
		if (tmp[p,GAP::.GapEnv$c_len]>50) {
			lines(c(tmp[p,GAP::.GapEnv$c_is],tmp[p,GAP::.GapEnv$c_if]),c(tmp[p,GAP::.GapEnv$c_lrr],tmp[p,GAP::.GapEnv$c_lrr])/2-1.5,col="red",type="l",lwd=2)
			lines(c(tmp[p,GAP::.GapEnv$c_is],tmp[p,GAP::.GapEnv$c_if]),c(tmp[p,GAP::.GapEnv$c_baf],tmp[p,GAP::.GapEnv$c_baf])+1.5,col="red",lwd=2)
			lines(c(tmp[p,GAP::.GapEnv$c_is],tmp[p,GAP::.GapEnv$c_if]),1-c(tmp[p,GAP::.GapEnv$c_baf],tmp[p,GAP::.GapEnv$c_baf])+1.5,col="red",lwd=2)
		}
	}
	lines(c(tmp[1,GAP::.GapEnv$c_is],tmp[dim(tmp)[[1]],GAP::.GapEnv$c_if]), c(-1,-1), col="grey30", type="l", lty=2)
	lines(c(tmp[1,GAP::.GapEnv$c_is],tmp[dim(tmp)[[1]],GAP::.GapEnv$c_if]), c(-1.5,-1.5), col="grey50", type="l", lty=2)
	lines(c(tmp[1,GAP::.GapEnv$c_is],tmp[dim(tmp)[[1]],GAP::.GapEnv$c_if]), c(-2,-2), col="grey30", type="l", lty=2)
	Line_H = c(1, 0.8, 0.6, 0.4, 0.2, 0, -0.2, -0.4)
	for (i in 1:length(Line_H)) {
		lines(c(tmp[1,GAP::.GapEnv$c_is],tmp[dim(tmp)[[1]],GAP::.GapEnv$c_if]), rep(Line_H[i],2),col=terrain_hcl(8,c=c(65,0),l=c(45, 90),power=c(1/2,1.5))[i],type="l", lty=1)
	}
	text(-9000,2.0, "B Allele\nFrequency", col="black", cex=.75, srt=90)
	text(-9000,0.35, "DNA Copy\nNumber", col="black", cex=.75, srt=90)
	text(-9000,-1.5, "Log 2 Ratio\n",col="black", cex=.95, srt=90)
	Line_A = c("0", ".2", ".4", ".6", ".8", "1")
	Line_H = c(1.5, 1.7, 1.9, 2.1, 2.3, 2.5)	
	for (i in 1:length(Line_H)) {
		text(-3500,Line_H[i],as.character(Line_A[i]), col="grey30", cex=.5, pos=4)
	}
	Line_A = 0:7
	Line_H = seq(from=-.4, to=1.0, by=.2)
	for (i in 1:length(Line_H)) {
		text(-3500,Line_H[i],as.character(Line_A[i]), col="grey30", cex=.5, pos=4)
	}
	text(-3500,-1.5, "0",col="grey30", cex=.5, pos=4)
	text(-3500,-1, "1",col="grey30", cex=.5, pos=4)
	text(-3500,-2, "-1",col="grey30", cex=.5, pos=4)
	tt = which(tmp[,GAP::.GapEnv$c_len]<100);if(length(tt)>0){tmp1 = tmp[-tt,]}
	for (p in 1:dim(tmp1)[1]) {
		lines(c(tmp1[p,GAP::.GapEnv$c_is],tmp1[p,GAP::.GapEnv$c_if]),c(tmp1[p,GAP::.GapEnv$c_CN],tmp1[p,GAP::.GapEnv$c_CN])/5-0.4,col="navy",type="l",lwd=4)
	}
	for (p in 2:dim(tmp1)[1]) {
		lines(c(tmp1[p-1,GAP::.GapEnv$c_if],tmp1[p,GAP::.GapEnv$c_is]),c(tmp1[p-1,GAP::.GapEnv$c_CN],tmp1[p,GAP::.GapEnv$c_CN])/5-0.4,col="navy",lwd=1)
	}
	for (p in 2:dim(tmp1)[1]) {
		if ((tmp1[p,GAP::.GapEnv$c_CN]==tmp1[p,GAP::.GapEnv$c_BA])||(tmp1[p,GAP::.GapEnv$c_BA]==0)) {
			lines(c(tmp1[p,GAP::.GapEnv$c_is],tmp1[p,GAP::.GapEnv$c_if]),c(tmp1[p,GAP::.GapEnv$c_CN],tmp1[p,GAP::.GapEnv$c_CN])/5-0.4,col="red",lwd=4)
		}
	}
	screen(2)
	split.screen(c(1,2))
	screen(3)
	plot(tmp[,GAP::.GapEnv$c_baf],tmp[,GAP::.GapEnv$c_lrr],main="Genome Alteration Print",xlim=c(0,1),ylim=c(-1,1),type="n",xlab="B Allele Frequency",ylab="Log 2 Ratio")
	CEX = tmp[,GAP::.GapEnv$c_len]/quantile(tmp[,GAP::.GapEnv$c_len],prob=0.7)
	tt = which(CEX>8)
	if (length(tt)>0) {
		CEX[tt] = 8
	}
	tt = which(CEX<0.04)
	if (length(tt)>0) {
		CEX[tt] = 0
	}
	points(1-tmp[,GAP::.GapEnv$c_baf],tmp[,GAP::.GapEnv$c_lrr],cex=CEX,col="black")
	points(tmp[,GAP::.GapEnv$c_baf],tmp[,GAP::.GapEnv$c_lrr],cex=CEX,col="grey")
	abline(h=0)
	abline(v=0.5)
	tt = which(tmp[,GAP::.GapEnv$c_len]<50)
	if (length(tt)>0) {
		tmp1 = tmp[-tt,]
	}
	RB = rainbow(40)
	for (p in 1:25) {
		tt1 = which(NETM[,3]==p)
		if (length(tt1)>0) {
			rect(min(NETM[tt1,2]),min(NETM[tt1,1]),max(NETM[tt1,2]),max(NETM[tt1,1]), col=RB[p],border=RB[p],density=80)
		}
	}
	for (p in 31:38) {
		tt1 = which(NETM[,3]==p)
		if (length(tt1)>0) {
			rect(min(NETM[tt1,2]),min(NETM[tt1,1]),max(NETM[tt1,2]),max(NETM[tt1,1]), col=RB[p],border=RB[p],density=80)
		}
	}
	for (j in 1:dim(Centr_CL)[2]) {
		if (Centr_CL[1,j]<=7) {
			if (Centr_CL[1,j]==Centr_CL[2,j]*2) {
				text(0.525,Centr_CL[5,j],paste(as.character(Centr_CL[1,j]),"/",as.character(Centr_CL[2,j]),sep=""),font=1)
			} else {
				text(Centr_CL[4,j],Centr_CL[5,j],paste(as.character(Centr_CL[1,j]),"/",as.character(Centr_CL[2,j]),sep=""),font=1)
			}
		}
	}
	text(0,-0.8,paste("p_BAF = ",as.character(round(p_BAF,2)),sep=""),cex=0.8, pos=4)
	text(0,-0.9,paste("q_LRR = ",as.character(round(q_LRR,2)),sep=""),cex=0.8, pos=4)
	text(0,-1,paste("2 Copy LRR = ",as.character(round(Centr_CL[5,4],3)),sep=""),cex=0.8, pos=4)
	DNAi = 0
	tt = 0
	k = 1
	for (k in 1:dim(tmp)[1]) {
		DNAi = DNAi+(SNP_annot[tmp[k,GAP::.GapEnv$c_if],C_Pos]-SNP_annot[tmp[k,GAP::.GapEnv$c_is],C_Pos])/1000*tmp[k,GAP::.GapEnv$c_CN]
		tt = tt+(SNP_annot[tmp[k,GAP::.GapEnv$c_if],C_Pos]-SNP_annot[tmp[k,GAP::.GapEnv$c_is],C_Pos])/1000
	}
	DNAi = round(DNAi/tt/2,2)
	HomoP = 0
	for (k in 1:dim(tmp)[1]) {
	  if ((tmp[k,GAP::.GapEnv$c_CN]==tmp[k,GAP::.GapEnv$c_BA])||(tmp[k,GAP::.GapEnv$c_BA]==0)) {
		HomoP = HomoP+(SNP_annot[tmp[k,GAP::.GapEnv$c_if],C_Pos]-SNP_annot[tmp[k,GAP::.GapEnv$c_is],C_Pos])/1000
	  }
	}
	HomoP = round(HomoP/tt*100,0)
	BaseP = 0
	BaseLossP = 0
	if (Chr_counts>60) {
		Ploidy = 2
	} else {
		Ploidy = 1
	}
	for (k in 1:dim(tmp)[1]) {
		if ((tmp[k,GAP::.GapEnv$c_CN]==2*Ploidy)&&(tmp[k,GAP::.GapEnv$c_BA]==1*Ploidy)) {
			BaseP = BaseP+(SNP_annot[tmp[k,GAP::.GapEnv$c_if],C_Pos]-SNP_annot[tmp[k,GAP::.GapEnv$c_is],C_Pos])/1000
	  	}
	  	if ((tmp[k,GAP::.GapEnv$c_CN]==1*Ploidy)&&(tmp[k,GAP::.GapEnv$c_BA]==1*Ploidy)) {
			BaseLossP = BaseLossP+(SNP_annot[tmp[k,GAP::.GapEnv$c_if],C_Pos]-SNP_annot[tmp[k,GAP::.GapEnv$c_is],C_Pos])/1000
	  	}
	}
	BaseP = round(BaseP/tt*100,0)
	BaseLossP = round(BaseLossP/tt*100,0)
	tt = which(tmp[,GAP::.GapEnv$c_len]<20)
	if (length(tt)>0) {
		tmp1 = tmp[-tt,]
	}
	tmp1 = shrinkReprTMP(tmp1)
	screen(4)
	plot(c(1,56),c(1,65),type="n",main="Summary and Centromeres",cex.main=1.3,xlab="",ylab="",col.axis="white",xaxt="n",yaxt="n")
	text(1,65,paste("Segment counts (>50)  = ",as.character(dim(tmp1)[[1]]-41)),cex=0.8,pos=4)
	text(1,62,paste("Chromosome counts = ",as.character(Chr_counts)),cex=0.8,pos=4)
	text(1,59,paste("Centromere counts = ",as.character(Centr_counts)),cex=0.8,pos=4)
	text(1,56,paste("DNA Index = ",as.character(DNAi)),cex=0.8,pos=4)
	text(30,65,paste("Ploidy detected (n) = ",as.character(Ploidy*2)),cex=0.8,pos=4)
	text(30,62,paste("Base state (%) = ",as.character(BaseP)),cex=0.8,pos=4)
	text(30,59,paste("'1' Copy loss (%) = ",as.character(BaseLossP)),cex=0.8,pos=4)
	text(30,56,paste("Acquired homozygous (%) = ",as.character(HomoP)),cex=0.8,pos=4)
	rect(0,46,57,50, col="grey80", border=NA)
	rect(0,0,57,50,col=NA,border="black")
	
	for (p in 1:2) {
		lines(c(7+(p-1)*30,7+(p-1)*30),c(0,50),col="black")
		lines(c(14+(p-1)*30,14+(p-1)*30),c(0,50),col="black")
		lines(c(21+(p-1)*30,21+(p-1)*30),c(0,50),col="black")
		lines(c(28+(1-1)*30,28+(1-1)*30),c(0,50),col="black")
		text(3+(p-1)*30,48,"Chr",cex=0.8)
		text(10+(p-1)*30,48,"Copy",cex=0.8)
		text(17+(p-1)*30,48,"  B Allele",cex=0.8)
		text(24+(p-1)*30,48,"SNPs",cex=0.8)
	}
	Centro_breaks = 0
	for (k in 1:dim(centr)[1]) {
		if (k<23) {
			p = 1
		} else {
			p = 2
		}
		if (trunc(centr[k,GAP::.GapEnv$c_chr])==centr[k,GAP::.GapEnv$c_chr]) {
			text(3+(p-1)*30,46*p-2*k-(p-1)*2,paste(as.character(trunc(centr[k,GAP::.GapEnv$c_chr])),"p"),cex=0.7)
			if ((k<dim(centr)[[1]])&&(trunc(centr[k+1,GAP::.GapEnv$c_chr])==centr[k,GAP::.GapEnv$c_chr])&&(centr[k,GAP::.GapEnv$c_CN]==centr[k+1,GAP::.GapEnv$c_CN])&&((centr[k,GAP::.GapEnv$c_BA]==centr[k+1,GAP::.GapEnv$c_BA])||(centr[k,GAP::.GapEnv$c_BA]*centr[k+1,GAP::.GapEnv$c_BA]==0))) {
				Centro_breaks = Centro_breaks+1
			}
		} else {
			text(3+(p-1)*30,46*p-2*k-(p-1)*2,paste(as.character(trunc(centr[k,GAP::.GapEnv$c_chr])),"q"),cex=0.7)
		}
		text(10+(p-1)*30,46*p-2*k-(p-1)*2,as.character(centr[k,GAP::.GapEnv$c_CN]),cex=0.7)
		text(17+(p-1)*30,46*p-2*k-(p-1)*2,as.character(centr[k,GAP::.GapEnv$c_BA]),cex=0.7)
		text(24+(p-1)*30,46*p-2*k-(p-1)*2,as.character(centr[k,GAP::.GapEnv$c_len]),cex=0.7)
	}
	close.screen(all = TRUE)
	summary_info = c(p_BAF,
					 q_LRR,
					 Centr_CL[5,4],
					 (dim(tmp1)[[1]]-41),
					 Chr_counts,
					 Centr_counts,
					 DNAi,
					 Ploidy*2)
	return(invisible(summary_info))
}
