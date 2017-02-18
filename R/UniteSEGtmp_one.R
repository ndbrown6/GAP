'UniteSEGtmp_one' <- function(tmp)
{
	FL = TRUE
	while (FL) {
		FL = FALSE
		for (k in 1:(dim(tmp)[1]-1)) {
			if ((tmp[k,GAP::.GapEnv$c_chr]==tmp[k+1,GAP::.GapEnv$c_chr])&&(tmp[k,GAP::.GapEnv$c_CN]==tmp[k+1,GAP::.GapEnv$c_CN])&&(tmp[k,GAP::.GapEnv$c_BA]==tmp[k+1,GAP::.GapEnv$c_BA])) {
				tmp[k+1,GAP::.GapEnv$c_is] = tmp[k,GAP::.GapEnv$c_is]
				tmp[k+1,GAP::.GapEnv$c_lrr] = weighted.mean(c(tmp[k+1,GAP::.GapEnv$c_lrr],tmp[k,GAP::.GapEnv$c_lrr]),c(tmp[k+1,GAP::.GapEnv$c_len],tmp[k,GAP::.GapEnv$c_len]))
				tmp[k+1,GAP::.GapEnv$c_nC] = tmp[k+1,GAP::.GapEnv$c_CN]
				tmp[k+1,GAP::.GapEnv$c_len] = tmp[k+1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
				tmp[k,1] = 0
				FL = TRUE
			} else {
				if ((tmp[k,GAP::.GapEnv$c_chr]==tmp[k+1,GAP::.GapEnv$c_chr])&&(tmp[k,GAP::.GapEnv$c_nC]==0)&&(tmp[k+1,GAP::.GapEnv$c_nC]==0)) {
					tmp[k+1,GAP::.GapEnv$c_is] = tmp[k,GAP::.GapEnv$c_is]
					tmp[k+1,GAP::.GapEnv$c_lrr] = weighted.mean(c(tmp[k+1,GAP::.GapEnv$c_lrr],tmp[k,GAP::.GapEnv$c_lrr]),c(tmp[k+1,GAP::.GapEnv$c_len],tmp[k,GAP::.GapEnv$c_len]))
					tmp[k+1,GAP::.GapEnv$c_nC] = tmp[k+1,GAP::.GapEnv$c_CN]
					tmp[k+1,GAP::.GapEnv$c_len] = tmp[k+1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
					tmp[k,1] = 0
					FL = TRUE
				}
			}
		}
		if (length(which(tmp[,1]==0))!=0) {
			tmp = tmp[-which(tmp[,1]==0),]
		}
	}
	return(invisible(tmp))
}
