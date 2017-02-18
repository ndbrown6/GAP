'fitSmallGermHomozyg' <- function(tmp, Length)
{
	FL = TRUE
	while (FL) {
		FL = FALSE
		tmp = UniteSEG_TMP(tmp)
		for (k in 2:(dim(tmp)[1]-1)) {
			if ((tmp[k-1,GAP::.GapEnv$c2_chr]==tmp[k+1,GAP::.GapEnv$c2_chr])&&(tmp[k,GAP::.GapEnv$c2_BA]==0)&&(tmp[k-1,GAP::.GapEnv$c2_CN]==tmp[k+1,GAP::.GapEnv$c2_CN])&&(tmp[k-1,GAP::.GapEnv$c2_BA]==tmp[k+1,GAP::.GapEnv$c2_BA])&&(tmp[k+1,GAP::.GapEnv$c2_BA]!=0)&&(tmp[k,GAP::.GapEnv$c2_CN]==tmp[k+1,GAP::.GapEnv$c2_CN])&&(tmp[k,GAP::.GapEnv$c2_len]<=Length)) {
				tmp[k,GAP::.GapEnv$c2_BA] = tmp[k+1,GAP::.GapEnv$c2_BA]
				FL = TRUE
			}
		}
		for (k in 2:(dim(tmp)[1]-1)) {
			if (tmp[k,GAP::.GapEnv$c2_BA]==0) {
				if ((tmp[k,GAP::.GapEnv$c2_chr]==tmp[k+1,GAP::.GapEnv$c2_chr])&&(tmp[k,GAP::.GapEnv$c2_CN]==tmp[k+1,GAP::.GapEnv$c2_CN])&&(tmp[k+1,GAP::.GapEnv$c2_BA]!=0)&&(tmp[k,GAP::.GapEnv$c2_len]<=Length)&&((tmp[k-1,GAP::.GapEnv$c2_CN]!=tmp[k,GAP::.GapEnv$c2_CN])||(tmp[k-1,GAP::.GapEnv$c2_chr]!=tmp[k,GAP::.GapEnv$c2_chr]))) {
					tmp[k,GAP::.GapEnv$c2_BA] = tmp[k+1,GAP::.GapEnv$c2_BA]
					FL = TRUE
				}
				if ((tmp[k,GAP::.GapEnv$c2_chr]==tmp[k-1,GAP::.GapEnv$c2_chr])&&(tmp[k,GAP::.GapEnv$c2_CN]==tmp[k-1,GAP::.GapEnv$c2_CN])&&(tmp[k-1,GAP::.GapEnv$c2_BA]!=0)&&(tmp[k,GAP::.GapEnv$c2_len]<=Length)&&((tmp[k+1,GAP::.GapEnv$c2_CN]!=tmp[k,GAP::.GapEnv$c2_CN])||(tmp[k+1,GAP::.GapEnv$c2_chr]!=tmp[k,GAP::.GapEnv$c2_chr]))) {
					tmp[k,GAP::.GapEnv$c2_BA] = tmp[k-1,GAP::.GapEnv$c2_BA]
					FL = TRUE
				}
			}
		}
	}
	return(invisible(tmp))
}
