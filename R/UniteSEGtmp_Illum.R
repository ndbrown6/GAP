'UniteSEGtmp_Illum' <- function(tmp, p_BAF)
{
	tmp = UniteSEGtmp_one(tmp)
	if (p_BAF>0.05) {
		FL = TRUE
		while (FL) {
			FL = FALSE
			for (k in 2:(dim(tmp)[1]-1)) {
				if ((tmp[k-1,GAP::.GapEnv$c_chr]==tmp[k+1,GAP::.GapEnv$c_chr])&&(tmp[k-1,GAP::.GapEnv$c_CN]==tmp[k+1,GAP::.GapEnv$c_CN])&&(tmp[k-1,GAP::.GapEnv$c_BA]==tmp[k+1,GAP::.GapEnv$c_BA])&&(tmp[k+1,GAP::.GapEnv$c_BA]!=0)&&(tmp[k,GAP::.GapEnv$c_CN]==tmp[k+1,GAP::.GapEnv$c_CN])&&(tmp[k,GAP::.GapEnv$c_BA]==0)&&(tmp[k,GAP::.GapEnv$c_len]<100)) {
					tmp[k,GAP::.GapEnv$c_BA] = tmp[k+1,GAP::.GapEnv$c_BA]
					FL = TRUE
				}
			}
			for (k in 2:(dim(tmp)[1]-1)) {
				if ((tmp[k,GAP::.GapEnv$c_BA]==0)&&(tmp[k,GAP::.GapEnv$c_len]<100)) {
					if ((tmp[k,GAP::.GapEnv$c_chr]==tmp[k+1,GAP::.GapEnv$c_chr])&&(tmp[k,GAP::.GapEnv$c_CN]==tmp[k+1,GAP::.GapEnv$c_CN])&&(tmp[k+1,GAP::.GapEnv$c_BA]!=0)&&((tmp[k-1,GAP::.GapEnv$c_CN]!=tmp[k,GAP::.GapEnv$c_CN])||(tmp[k-1,GAP::.GapEnv$c_chr]!=tmp[k,GAP::.GapEnv$c_chr]))) {
						tmp[k,GAP::.GapEnv$c_BA] = tmp[k+1,GAP::.GapEnv$c_BA]
						FL = TRUE
					}
					if ((tmp[k,GAP::.GapEnv$c_chr]==tmp[k-1,GAP::.GapEnv$c_chr])&&(tmp[k,GAP::.GapEnv$c_CN]==tmp[k-1,GAP::.GapEnv$c_CN])&&(tmp[k-1,GAP::.GapEnv$c_BA]!=0)&&((tmp[k+1,GAP::.GapEnv$c_CN]!=tmp[k,GAP::.GapEnv$c_CN])||(tmp[k+1,GAP::.GapEnv$c_chr]!=tmp[k,GAP::.GapEnv$c_chr]))) {
						tmp[k,GAP::.GapEnv$c_BA] = tmp[k-1,GAP::.GapEnv$c_BA]
						FL = TRUE
					}
				}
			}
			if ((tmp[1,GAP::.GapEnv$c_BA]==0)&&(tmp[k,GAP::.GapEnv$c_len]<100)) {
				if ((tmp[1,GAP::.GapEnv$c_chr]==tmp[2,GAP::.GapEnv$c_chr])&&(tmp[1,GAP::.GapEnv$c_CN]==tmp[2,GAP::.GapEnv$c_CN])&&(tmp[2,GAP::.GapEnv$c_BA]!=0)) {
					tmp[1,GAP::.GapEnv$c_BA] = tmp[2,GAP::.GapEnv$c_BA]
					FL = TRUE
				}
			}
			k = dim(tmp)[1]
			if ((tmp[k,GAP::.GapEnv$c_BA]==0)&&(tmp[k,GAP::.GapEnv$c_len]<100)) {
				if ((tmp[k,GAP::.GapEnv$c_chr]==tmp[k-1,GAP::.GapEnv$c_chr])&&(tmp[k,GAP::.GapEnv$c_CN]==tmp[k-1,GAP::.GapEnv$c_CN])&&(tmp[k-1,GAP::.GapEnv$c_BA]!=0)) {
					tmp[k,GAP::.GapEnv$c_BA] = tmp[k-1,GAP::.GapEnv$c_BA]
					FL = TRUE
				}
			}
			tmp = UniteSEGtmp_one(tmp)
		}
	}
	return(invisible(tmp))
}
