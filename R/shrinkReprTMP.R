'shrinkReprTMP' <- function(tmp)
{
	tmp[1,GAP::.GapEnv$c_lrr] = tmp[1,GAP::.GapEnv$c_len]*tmp[1,GAP::.GapEnv$c_lrr]
	for (k in 1:(dim(tmp)[1]-1)) {
		if (tmp[k,GAP::.GapEnv$c_conf]==tmp[k+1,GAP::.GapEnv$c_conf]) {
			tmp[k+1,GAP::.GapEnv$c_is] = tmp[k,GAP::.GapEnv$c_is]
			tmp[k+1,GAP::.GapEnv$c_lrr] = tmp[k,GAP::.GapEnv$c_lrr]+tmp[k+1,GAP::.GapEnv$c_len]*tmp[k+1,GAP::.GapEnv$c_lrr]
			tmp[k+1,GAP::.GapEnv$c_len] = tmp[k+1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
			tmp[k,1] = 0
			FL = TRUE
		} else {
			tmp[k,GAP::.GapEnv$c_lrr] = tmp[k,GAP::.GapEnv$c_lrr]/tmp[k,GAP::.GapEnv$c_len]
			tmp[k+1,GAP::.GapEnv$c_lrr] = tmp[k+1,GAP::.GapEnv$c_len]*tmp[k+1,GAP::.GapEnv$c_lrr]
		}
	}
	k = dim(tmp)[1]
	tmp[k,GAP::.GapEnv$c_lrr] = tmp[k,GAP::.GapEnv$c_lrr]/tmp[k,GAP::.GapEnv$c_len]
	tt = which(tmp[,1]==0)
	if (length(tt)!=0) {
		tmp = tmp[-tt,]
	}
	return(invisible(tmp))
}
