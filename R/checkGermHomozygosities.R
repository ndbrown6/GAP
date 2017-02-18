'checkGermHomozygosities' <- function(tmp, Length)
{
	for (k in 2:(dim(tmp)[1]-1)) {
		if ((tmp[k,GAP::.GapEnv$c2_BA]==0)&&(tmp[k-1,GAP::.GapEnv$c2_chr]==tmp[k+1,GAP::.GapEnv$c2_chr])&&(tmp[k-1,GAP::.GapEnv$c2_CN]!=tmp[k,GAP::.GapEnv$c2_CN])&&(tmp[k+1,GAP::.GapEnv$c2_CN]!=tmp[k,GAP::.GapEnv$c2_CN])) {
			tmp[k,GAP::.GapEnv$c2_BA] = tmp[k,GAP::.GapEnv$c2_CN]
		}
	}
	for (k in 2:(dim(tmp)[1]-1)) {
		if ((tmp[k,GAP::.GapEnv$c2_BA]==0)&&(tmp[k,GAP::.GapEnv$c2_len]>=Length)) {
			tmp[k,GAP::.GapEnv$c2_BA] = tmp[k,GAP::.GapEnv$c2_CN]
		}
	}
	return(invisible(tmp))
}
