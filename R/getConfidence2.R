'getConfidence2' <- function(tmp)
{
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
