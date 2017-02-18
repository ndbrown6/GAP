'UniteSEG_TMP' <- function(tmp)
{
	FL = TRUE
	while (FL) {
		FL = FALSE
		for (k in 1:(dim(tmp)[[1]]-1)) {
			if ((tmp[k,GAP::.GapEnv$c_chr]==tmp[k+1,GAP::.GapEnv$c_chr])&&(trunc(tmp[k,GAP::.GapEnv$c_CN])==tmp[k,GAP::.GapEnv$c_CN])&&(tmp[k,GAP::.GapEnv$c_CN]==tmp[k+1,GAP::.GapEnv$c_CN])&&(tmp[k,GAP::.GapEnv$c_BA]==tmp[k+1,GAP::.GapEnv$c_BA])) {
				tmp[k+1,GAP::.GapEnv$c_is] = tmp[k,GAP::.GapEnv$c_is]
				w = c(tmp[k+1,GAP::.GapEnv$c_len],tmp[k,GAP::.GapEnv$c_len])
				tmp[k+1,c_cn] = weighted.mean(c(tmp[k+1,c_cn],tmp[k,c_cn]),w)
				tmp[k+1,GAP::.GapEnv$c_nC] = weighted.mean(c(tmp[k+1,GAP::.GapEnv$c_nC],tmp[k,GAP::.GapEnv$c_nC]),w)
				tmp[k+1,GAP::.GapEnv$c_CN1:(GAP::.GapEnv$c_CN1+8)] = (-1)
				tmp[k+1,GAP::.GapEnv$c_CN1] = tmp[k+1,GAP::.GapEnv$c_CN]
				tmp[k+1,GAP::.GapEnv$c_BA1] = tmp[k+1,GAP::.GapEnv$c_BA]
				tmp[k+1,GAP::.GapEnv$c_rnk1] = 1
				tmp[k+1,GAP::.GapEnv$c_len] = tmp[k+1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
				tmp[k+1,GAP::.GapEnv$c_conf] = max(tmp[k+1,GAP::.GapEnv$c_conf],tmp[k,GAP::.GapEnv$c_conf])	
				tmp[k,GAP::.GapEnv$c_ind] = 0
			} else {
				if ((tmp[k,GAP::.GapEnv$c_chr]==tmp[k+1,GAP::.GapEnv$c_chr])&&(tmp[k,GAP::.GapEnv$c_CN]==0)&&(tmp[k+1,GAP::.GapEnv$c_CN]==0)) {
					tmp[k+1,GAP::.GapEnv$c_is] = tmp[k,GAP::.GapEnv$c_is]
					w = c(tmp[k+1,GAP::.GapEnv$c_len],tmp[k,GAP::.GapEnv$c_len])
					tmp[k+1,c_cn] = weighted.mean(c(tmp[k+1,c_cn],tmp[k,c_cn]),w)
					tmp[k+1,GAP::.GapEnv$c_nC] = weighted.mean(c(tmp[k+1,GAP::.GapEnv$c_nC],tmp[k,GAP::.GapEnv$c_nC]),w)
					tmp[k+1,GAP::.GapEnv$c_CN1:(GAP::.GapEnv$c_CN1+8)] = (-1)
					tmp[k+1,GAP::.GapEnv$c_CN1] = tmp[k+1,GAP::.GapEnv$c_CN]
					tmp[k+1,GAP::.GapEnv$c_BA1] = tmp[k+1,GAP::.GapEnv$c_BA]
					tmp[k+1,GAP::.GapEnv$c_rnk1] = 1
					tmp[k+1,GAP::.GapEnv$c_len] = tmp[k+1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
					tmp[k+1,GAP::.GapEnv$c_conf] = max(tmp[k+1,GAP::.GapEnv$c_conf],tmp[k,GAP::.GapEnv$c_conf])	
					tmp[k,GAP::.GapEnv$c_ind] = 0
				}
			}
		}
		tt = which(tmp[,GAP::.GapEnv$c_ind]==0)
		if (length(tt)!=0) {
			tmp = tmp[-tt,]
			FL = TRUE
		}
	}
	tmp = round(tmp,4)
	tmp[,GAP::.GapEnv$c_nC] = round(tmp[,GAP::.GapEnv$c_nC],1)
	return(invisible(tmp))
}