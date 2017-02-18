'checkSmallDiscrepancies' <- function(tmp, Length)
{
	for (k in 2:(dim(tmp)[1]-1)) {
		if ((tmp[k,GAP::.GapEnv$c2_len]<Length)&&(tmp[k,GAP::.GapEnv$c2_CN]>0)&&(max(tmp[k+1,GAP::.GapEnv$c2_len],tmp[k-1,GAP::.GapEnv$c2_len])>=2*Length)) {
			if ((tmp[k-1,GAP::.GapEnv$c2_chr]==tmp[k+1,GAP::.GapEnv$c2_chr])&&(tmp[k-1,GAP::.GapEnv$c2_CN]==tmp[k+1,GAP::.GapEnv$c2_CN])&&(tmp[k-1,GAP::.GapEnv$c2_BA]==tmp[k+1,GAP::.GapEnv$c2_BA])) {
				if((abs(tmp[k,GAP::.GapEnv$c2_nC]-tmp[k+1,GAP::.GapEnv$c2_nC])<1)||(abs(tmp[k,GAP::.GapEnv$c2_nC]-tmp[k-1,GAP::.GapEnv$c2_nC])<1)) {
					tmp[k,GAP::.GapEnv$c2_CN] = tmp[k-1,GAP::.GapEnv$c2_CN]
					tmp[k,GAP::.GapEnv$c2_BA] = tmp[k-1,GAP::.GapEnv$c2_BA]
				}
			}
		}
	}
	for (k in 2:(dim(tmp)[1]-1)) {
  		if ((tmp[k,GAP::.GapEnv$c2_len]<Length)&&(tmp[k,GAP::.GapEnv$c2_CN]>0)&&(tmp[k-1,GAP::.GapEnv$c2_chr]==tmp[k+1,GAP::.GapEnv$c2_chr])&&(tmp[k-1,GAP::.GapEnv$c2_CN]!=tmp[k+1,GAP::.GapEnv$c2_CN])) {
			if ((tmp[k+1,GAP::.GapEnv$c2_len]>=2*Length)&&(abs(tmp[k,GAP::.GapEnv$c2_nC]-tmp[k+1,GAP::.GapEnv$c2_nC])<=1)&&(abs(tmp[k,GAP::.GapEnv$c2_nC]-tmp[k-1,GAP::.GapEnv$c2_nC])>1)) {
				tmp[k,GAP::.GapEnv$c2_CN] = tmp[k+1,GAP::.GapEnv$c2_CN]
				tmp[k,GAP::.GapEnv$c2_BA] = tmp[k+1,GAP::.GapEnv$c2_BA]
			}
			if ((tmp[k-1,GAP::.GapEnv$c2_len]>=2*Length)&&(abs(tmp[k,GAP::.GapEnv$c2_nC]-tmp[k+1,GAP::.GapEnv$c2_nC])>1)&&(abs(tmp[k,GAP::.GapEnv$c2_nC]-tmp[k-1,GAP::.GapEnv$c2_nC])<=1)) {
				tmp[k,GAP::.GapEnv$c2_CN] = tmp[k-1,GAP::.GapEnv$c2_CN]
				tmp[k,GAP::.GapEnv$c2_BA] = tmp[k-1,GAP::.GapEnv$c2_BA]
			}
		}
	}
	for (k in 2:(dim(tmp)[1]-1)) {
		if ((tmp[k+1,GAP::.GapEnv$c2_chr]!=tmp[k,GAP::.GapEnv$c2_chr])&&(tmp[k,GAP::.GapEnv$c2_len]<Length)&&(tmp[k,GAP::.GapEnv$c2_CN]>0)) {
			if ((tmp[k-1,GAP::.GapEnv$c2_chr]==tmp[k,GAP::.GapEnv$c2_chr])&&(tmp[k-1,GAP::.GapEnv$c2_len]>=2*Length)&&(abs(tmp[k,GAP::.GapEnv$c2_nC]-tmp[k-1,GAP::.GapEnv$c2_nC])<1)) {
				tmp[k,GAP::.GapEnv$c2_CN] = tmp[k-1,GAP::.GapEnv$c2_CN]
				tmp[k,GAP::.GapEnv$c2_BA] = tmp[k-1,GAP::.GapEnv$c2_BA]
			}
		}
	}
	for (k in 2:(dim(tmp)[1]-1)) {
		if ((tmp[k-1,GAP::.GapEnv$c2_chr]!=tmp[k,GAP::.GapEnv$c2_chr])&&(tmp[k,GAP::.GapEnv$c2_len]<Length)&&(tmp[k,GAP::.GapEnv$c2_CN]>0)) {
			if ((tmp[k+1,GAP::.GapEnv$c2_chr]==tmp[k,GAP::.GapEnv$c2_chr])&&(tmp[k+1,GAP::.GapEnv$c2_len]>=2*Length)&&(abs(tmp[k,GAP::.GapEnv$c2_nC]-tmp[k+1,GAP::.GapEnv$c2_nC])<1)) {
				tmp[k,GAP::.GapEnv$c2_CN] = tmp[k+1,GAP::.GapEnv$c2_CN]
				tmp[k,GAP::.GapEnv$c2_BA] = tmp[k+1,GAP::.GapEnv$c2_BA]
			}
		}
	}
	tmp = UniteSEG_TMP(tmp)
	tmp[,GAP::.GapEnv$c2_conf] = 0
	tmp[1,GAP::.GapEnv$c2_conf] = 1
	for (k in 2:dim(tmp)[1]) {
		if ((tmp[k,GAP::.GapEnv$c2_chr]==tmp[k-1,GAP::.GapEnv$c2_chr])&&(tmp[k,GAP::.GapEnv$c2_CN]==tmp[k-1,GAP::.GapEnv$c2_CN])&&(tmp[k,GAP::.GapEnv$c2_BA]==tmp[k-1,GAP::.GapEnv$c2_BA])) {
			tmp[k,GAP::.GapEnv$c2_conf] = tmp[k-1,GAP::.GapEnv$c2_conf]
		} else {
			tmp[k,GAP::.GapEnv$c2_conf] = tmp[k-1,GAP::.GapEnv$c2_conf]+1
		}
	}
	return(invisible(tmp))
}
