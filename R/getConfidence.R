'getConfidence' <- function(tmp)
{
	tmp[,GAP::.GapEnv$c_conf] = 0
	tmp[1,GAP::.GapEnv$c_conf] = 1
	for (k in 2:dim(tmp)[1]) {
		if ((tmp[k,GAP::.GapEnv$c_chr]==tmp[k-1,GAP::.GapEnv$c_chr])&&(tmp[k,GAP::.GapEnv$c_CN]==tmp[k-1,GAP::.GapEnv$c_CN])&&(tmp[k,GAP::.GapEnv$c_BA]==tmp[k-1,GAP::.GapEnv$c_BA])) {
			tmp[k,GAP::.GapEnv$c_conf] = tmp[k-1,GAP::.GapEnv$c_conf]
		} else {
			tmp[k,GAP::.GapEnv$c_conf] = tmp[k-1,GAP::.GapEnv$c_conf]+1
		}
	}
	return(invisible(tmp))
}
