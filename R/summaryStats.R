'summaryStats' <- function(tmp, SNP_annot)
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
	for (k in 1:(dim(tmp1)[1]-1)) {
		if ((tmp1[k+1,GAP::.GapEnv$c_chr]>tmp1[k,GAP::.GapEnv$c_chr])&&(round((tmp1[k+1,GAP::.GapEnv$c_chr]-0.2),0)==tmp1[k,GAP::.GapEnv$c_chr])) {
			Chr_counts = Chr_counts+(tmp1[k,GAP::.GapEnv$c_CN]+tmp1[k+1,GAP::.GapEnv$c_CN])/2
			Centr_counts = Centr_counts+tmp1[k,GAP::.GapEnv$c_CN]+tmp1[k+1,GAP::.GapEnv$c_CN]
		} else {
			tt = round((tmp1[k+1,GAP::.GapEnv$c_chr]-0.2),0)
			if (tmp1[k+1,GAP::.GapEnv$c_chr]>tmp1[k,GAP::.GapEnv$c_chr]) {
				if ((tt==13)||(tt==14)||(tt==15)||(tt==22)) {
					Chr_counts = Chr_counts+tmp1[k+1,GAP::.GapEnv$c_CN]
					Centr_counts = Centr_counts+tmp1[k+1,GAP::.GapEnv$c_CN]
				}
			}
		}
	}
	DNAi = 0
	tt = 0
	k = 1
	for (k in 1:dim(tmp)[1]) {
		DNAi = DNAi+(SNP_annot[tmp[k,GAP::.GapEnv$c_if],C_Pos]-SNP_annot[tmp[k,GAP::.GapEnv$c_is],C_Pos])/1000*tmp[k,GAP::.GapEnv$c_CN]
		tt = tt+(SNP_annot[tmp[k,GAP::.GapEnv$c_if],C_Pos]-SNP_annot[tmp[k,GAP::.GapEnv$c_is],C_Pos])/1000
	}
	DNAi = round(DNAi/tt/2,2)
	if (Chr_counts>60) {
		Ploidy = 2
	} else {
		Ploidy = 1
	}
	summary_stats = c(Chr_counts, Centr_counts, DNAi, Ploidy*2)
	return(invisible(summary_stats))
}
