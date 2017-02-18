'fitShortSubclones' <- function(tmp, Centr_CL, NETM, p_BAF, Length)
{
	TMP = tmp
	tt = which(TMP[,GAP::.GapEnv$c_len]>=Length)
	tt1 = table(TMP[,GAP::.GapEnv$c_chr])
	if (length(which(tt1==1))>0) {
		tt1 = tt1[which(tt1==1)]
		tt1 = as.numeric(rownames(tt1))
		for (j in 1:length(tt1)) {
			tt = c(tt,which(TMP[,GAP::.GapEnv$c_chr]==tt1[j]))
		}
		tt = as.numeric(rownames(table(tt)))
	}
	if (length(tt)>0) {
		tMP = TMP[-tt,]
		TMP = TMP[tt,]
	}
	TMP = UniteSEGtmp_Illum(TMP, p_BAF)
	TMP = adjustTMPstableSubclone(TMP,Centr_CL,NETM,p_BAF)
	TMP = UniteSEGtmp_Illum(TMP,p_BAF)
	tmp[,GAP::.GapEnv$c_conf] = 0
	for (k in 1:dim(TMP)[1]) {
		tt = which(tmp[,GAP::.GapEnv$c_is]>=TMP[k,GAP::.GapEnv$c_is])
		if (length(tt)>0) {
			kk = which(tmp[tt,GAP::.GapEnv$c_if]<=TMP[k,GAP::.GapEnv$c_if])
			if (length(kk)>0) {
				tmp[tt[kk],GAP::.GapEnv$c_CN] = TMP[k,GAP::.GapEnv$c_CN]
				tmp[tt[kk],GAP::.GapEnv$c_BA] = TMP[k,GAP::.GapEnv$c_BA]
				tmp[tt[kk],GAP::.GapEnv$c_conf] = k
			}
		}
	}
	for (k in 1:dim(tMP)[1]) {
		tt = which(tmp[,GAP::.GapEnv$c_ind]==tMP[k,GAP::.GapEnv$c_ind])
		if (length(tt)>0) {
			tmp[tt,GAP::.GapEnv$c_conf] = -k
		}
	}
	return(invisible(tmp))
}
