'shrinkReprTMP2' <- function(tmp)
{
	tmp[1,GAP::.GapEnv$c2_cn] = tmp[1,GAP::.GapEnv$c2_len]*tmp[1,GAP::.GapEnv$c2_cn]
	for (k in 1:(dim(tmp)[1]-1)) {
		if (tmp[k,GAP::.GapEnv$c2_conf]==tmp[k+1,GAP::.GapEnv$c2_conf]) {
			tmp[k+1,GAP::.GapEnv$c2_is] = tmp[k,GAP::.GapEnv$c2_is]
			tmp[k+1,GAP::.GapEnv$c2_CN1:(GAP::.GapEnv$c2_CN1+8)] = (-1)
			tmp[k+1,GAP::.GapEnv$c2_CN1] = tmp[k+1,GAP::.GapEnv$c2_CN]
			tmp[k+1,GAP::.GapEnv$c2_BA1] = tmp[k+1,GAP::.GapEnv$c2_BA]
			tmp[k+1,GAP::.GapEnv$c2_rnk1] = 1
			tmp[k,1] = 0
			FL = TRUE
		}
	}
	k = dim(tmp)[1]
	tt = which(tmp[,1]==0)
	if (length(tt)!=0) {
		tmp = tmp[-tt,]
	}
	return(invisible(tmp))
}
