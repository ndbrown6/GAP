'iliminateShort_FinalPolish' <- function(tmp)
{
	for (k in 2:(dim(tmp)[1]-1)) {
		if ((tmp[k,GAP::.GapEnv$c_conf]!=0)&&(tmp[k,GAP::.GapEnv$c_len]<10)&&(tmp[k-1,GAP::.GapEnv$c_chr]==tmp[k+1,GAP::.GapEnv$c_chr])) {
			if ((tmp[k,GAP::.GapEnv$c_conf]==tmp[k-1,GAP::.GapEnv$c_conf])&&(tmp[k,GAP::.GapEnv$c_conf]==tmp[k+1,GAP::.GapEnv$c_conf])) {
				if ((is.na(tmp[k,GAP::.GapEnv$c_lrr]))||(abs(tmp[k,GAP::.GapEnv$c_lrr]-tmp[k+1,GAP::.GapEnv$c_lrr])>abs(tmp[k,GAP::.GapEnv$c_lrr]-tmp[k-1,GAP::.GapEnv$c_lrr]))) {
					tmp[k-1,GAP::.GapEnv$c_if] = tmp[k,GAP::.GapEnv$c_if]
					tmp[k-1,GAP::.GapEnv$c_len] = tmp[k-1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
					tmp[k,] = tmp[k-1,]
					tmp[k-1,1] = 0
				} else {
					if (abs(tmp[k,GAP::.GapEnv$c_lrr]-tmp[k+1,GAP::.GapEnv$c_lrr])<=abs(tmp[k,GAP::.GapEnv$c_lrr]-tmp[k-1,GAP::.GapEnv$c_lrr])) {
						tmp[k+1,GAP::.GapEnv$c_is] = tmp[k,GAP::.GapEnv$c_is]
						tmp[k+1,GAP::.GapEnv$c_len] = tmp[k+1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
						tmp[k,1] = 0
					}
				}
			}
		}
	}
	tt = which(tmp[,1]==0)
	if (length(tt)!=0) {
		tmp = tmp[-tt,]
	}
	for (k in 2:dim(tmp)[1]) {
		if ((tmp[k,GAP::.GapEnv$c_conf]!=0)&&(tmp[k,GAP::.GapEnv$c_len]<10)&&(tmp[k,GAP::.GapEnv$c_chr]==tmp[k-1,GAP::.GapEnv$c_chr])&&(tmp[k,GAP::.GapEnv$c_conf]==tmp[k-1,GAP::.GapEnv$c_conf])) {
			tmp[k-1,GAP::.GapEnv$c_if] = tmp[k,GAP::.GapEnv$c_if]
			tmp[k-1,GAP::.GapEnv$c_len] = tmp[k-1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
			tmp[k,] = tmp[k-1,]
			tmp[k-1,1] = 0
		}
	}
	tt = which(tmp[,1]==0)
	if (length(tt)!=0) {
		tmp = tmp[-tt,]
	}
	for (k in 1:(dim(tmp)[1]-1)) {
		if ((tmp[k,GAP::.GapEnv$c_conf]!=0)&&(tmp[k,GAP::.GapEnv$c_len]<10)&&(tmp[k,GAP::.GapEnv$c_chr]==tmp[k+1,GAP::.GapEnv$c_chr])&&(tmp[k,GAP::.GapEnv$c_conf]==tmp[k+1,GAP::.GapEnv$c_conf])) {
			tmp[k+1,GAP::.GapEnv$c_is] = tmp[k,GAP::.GapEnv$c_is]
			tmp[k+1,GAP::.GapEnv$c_len] = tmp[k+1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
			tmp[k,1] = 0
		}
	}
	tt = which(tmp[,1]==0)
	if (length(tt)!=0) {
		tmp = tmp[-tt,]
	}
	k = 412
	for (k in 2:(dim(tmp)[1]-1)) {
		if ((tmp[k,GAP::.GapEnv$c_conf]!=0)&&(tmp[k,GAP::.GapEnv$c_len]<10)&&(tmp[k-1,GAP::.GapEnv$c_chr]==tmp[k+1,GAP::.GapEnv$c_chr])) {
			if (is.na(tmp[k+1,GAP::.GapEnv$c_lrr])) {
				tmp[k,GAP::.GapEnv$c_if] = tmp[k+1,GAP::.GapEnv$c_if]
				tmp[k,GAP::.GapEnv$c_len] = tmp[k+1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
				tmp[k+1,] = tmp[k,]
				tmp[k,] = tmp[k-1,]
				tmp[k-1,1] = 0
			} else {
				if ((is.na(tmp[k,GAP::.GapEnv$c_lrr]))||(abs(tmp[k,GAP::.GapEnv$c_lrr]-tmp[k+1,GAP::.GapEnv$c_lrr])>abs(tmp[k,GAP::.GapEnv$c_lrr]-tmp[k-1,GAP::.GapEnv$c_lrr]))) {
					tmp[k-1,GAP::.GapEnv$c_if] = tmp[k,GAP::.GapEnv$c_if]
					tmp[k-1,GAP::.GapEnv$c_len] = tmp[k-1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
					tmp[k,] = tmp[k-1,]
					tmp[k-1,1] = 0
				} else {
					if (abs(tmp[k,GAP::.GapEnv$c_lrr]-tmp[k+1,GAP::.GapEnv$c_lrr])<=abs(tmp[k,GAP::.GapEnv$c_lrr]-tmp[k-1,GAP::.GapEnv$c_lrr])) {
						tmp[k+1,GAP::.GapEnv$c_is] = tmp[k,GAP::.GapEnv$c_is]
						tmp[k+1,GAP::.GapEnv$c_len] = tmp[k+1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
						tmp[k,1] = 0
					}
				}
			}
		}
	}
	tt = which(tmp[,1]==0)
	if (length(tt)!=0) {
		tmp = tmp[-tt,]
	}
	tmp[,1] = seq(1,dim(tmp)[1])
	tt = which(tmp[,GAP::.GapEnv$c_conf]<=0)
	if (length(tt)>0) {
		tmp[tt,GAP::.GapEnv$c_conf] = -seq(1,length(tt))
	}
	return(invisible(tmp))
}
