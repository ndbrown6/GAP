'adjustTMPstableSubclone' <- function(tmp, Centr_CL, NETM, p_BAF)
{
	tt = which(trunc(tmp[,GAP::.GapEnv$c_nC])-tmp[,GAP::.GapEnv$c_nC]!=0)
	if (length(tt)>0) {
		for (j in 1:length(tt)) {
			if (tmp[tt[j],GAP::.GapEnv$c_len]>500) {
				tmp[tt[j],GAP::.GapEnv$c_nC] = tmp[tt[j],GAP::.GapEnv$c_CN]
			}
			if (tt[j]==dim(tmp)[[1]]) {
				if (tmp[tt[j],GAP::.GapEnv$c_chr]!=tmp[tt[j]-1,GAP::.GapEnv$c_chr]) {
					tmp[tt[j],GAP::.GapEnv$c_nC] = tmp[tt[j],GAP::.GapEnv$c_CN]
				}
			} else {
				if (tt[j]==1) {
					if (tmp[tt[j],GAP::.GapEnv$c_chr]!=tmp[tt[j]+1,GAP::.GapEnv$c_chr]) {
						tmp[tt[j],GAP::.GapEnv$c_nC] = tmp[tt[j],GAP::.GapEnv$c_CN]
					}
				} else {
					if ((tmp[tt[j],GAP::.GapEnv$c_chr]!=tmp[tt[j]+1,GAP::.GapEnv$c_chr])&&(tmp[tt[j],GAP::.GapEnv$c_chr]!=tmp[tt[j]-1,GAP::.GapEnv$c_chr])) {
						tmp[tt[j],GAP::.GapEnv$c_nC] = tmp[tt[j],GAP::.GapEnv$c_CN]
					}
				}
			}
		}
	}
	tt = which(trunc(tmp[,GAP::.GapEnv$c_nC])-tmp[,GAP::.GapEnv$c_nC]!=0)
	if (length(tt)>0) {
		tmpB = tmp[-tt,]
		if (length(tt)==1) {
			tmpS = rbind(tmp[tt,],tmp[tt,])
		} else {
			tmpS = tmp[tt,]
		}
		tmpB = UniteSEGtmp_Illum(tmpB,p_BAF)
		tmpS[,GAP::.GapEnv$c_bafH] = 0
		j = 1
		for (j in 1:dim(tmpS)[1]) {
			tt = which(tmpB[,GAP::.GapEnv$c_is]<tmpS[j,GAP::.GapEnv$c_is])
			if (length(tt)>0) {
				tt = tt[length(tt)]
				if (tmpB[tt,GAP::.GapEnv$c_if]>tmpS[j,GAP::.GapEnv$c_if]) {
					tmpS[j,GAP::.GapEnv$c_bafH] = tt
				} else {
					if (tt==dim(tmpB)[[1]]) {
						tmpS[j,GAP::.GapEnv$c_bafH] = tt
					} else {
						tmpS[j,GAP::.GapEnv$c_bafH] = tt+0.5
					}
				}
			} else {
				tt = 0
				tmpS[j,GAP::.GapEnv$c_bafH] = tt+0.5
			}
		}
		tt = table(tmpS[,GAP::.GapEnv$c_bafH])
		if (length(tt)==1) {
			tt = c(tt,tt)
			tt = rbind(as.numeric(names(tt)),tt)
			ONE = TRUE
		} else {
			tt = rbind(as.numeric(rownames(tt)),tt)
			ONE = FALSE
		}
		colnames(tt) = seq(1,dim(tt)[2])	
		TT = 1
		j = 1
		if (ONE) {
			dim_tt = 1
		} else {
			dim_tt = dim(tt)[2]
		}
		for (j in 1:dim_tt) {
			if (trunc(tt[1,j])==tt[1,j]) {
				if (tt[2,j]==1) {
					if (abs(tmpS[TT,GAP::.GapEnv$c_nC]-tmpB[tt[1,j],GAP::.GapEnv$c_CN])<0.9) {
						tmpS[TT,GAP::.GapEnv$c_nC] = tmpB[tt[1,j],GAP::.GapEnv$c_CN]
					} else {
						tmpS[TT,GAP::.GapEnv$c_nC] = tmpS[TT,GAP::.GapEnv$c_CN]
					}
				} else {
					if (max(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpB[tt[1,j],GAP::.GapEnv$c_CN]))<0.9) {
						tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC] = tmpB[tt[1,j],GAP::.GapEnv$c_CN]
					} else {
						tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC] = tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_CN]
					}
				}
				TT = TT+tt[2,j]
			} else {
				kk = trunc(tt[1,j])
				if (kk==0) {
					kk = 1
					if (tt[2,j]==1) {
						if (abs(tmpS[TT,GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN])<0.9) {
							tmpS[TT,GAP::.GapEnv$c_nC] = tmpB[kk,GAP::.GapEnv$c_CN]
						} else {
							tmpS[TT,GAP::.GapEnv$c_nC] = tmpS[TT,GAP::.GapEnv$c_CN]
						}
					} else {
						if (max(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN]))<0.9) {
							tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC] = tmpB[kk,GAP::.GapEnv$c_CN]
						} else {
							tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC] = tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_CN]
						}
					}
					TT = TT+tt[2,j]
				} else {
					if (trunc(tmpB[kk,GAP::.GapEnv$c_chr])==trunc(tmpB[kk+1,GAP::.GapEnv$c_chr])) {
						if (tt[2,j]==1) {
							kk = which(abs(tmpS[TT,GAP::.GapEnv$c_nC]-tmpB[kk:(kk+1),GAP::.GapEnv$c_CN])<0.9)
							if (length(kk)==0) {
								tmpS[TT,GAP::.GapEnv$c_nC] = tmpS[TT,GAP::.GapEnv$c_CN]
							}
							if (length(kk)==1) {
								tmpS[TT,GAP::.GapEnv$c_nC] = tmpB[trunc(tt[1,j])+kk-1,GAP::.GapEnv$c_CN]
							}
							if (length(kk)==2) {
								tmpS[TT,GAP::.GapEnv$c_nC] = tmpS[TT,GAP::.GapEnv$c_CN]
							}
						} else {
							if (max(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN]))<0.9) {
								if (min(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpB[kk+1,GAP::.GapEnv$c_CN]))>0.9) {
									tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC] = tmpB[kk,GAP::.GapEnv$c_CN]
								} else {
									tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC] = tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_CN]
								}
							} else {
								if (min(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN]))>0.9) {
									if (max(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpB[kk+1,GAP::.GapEnv$c_CN]))<0.9) {
										tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC] = tmpB[kk+1,GAP::.GapEnv$c_CN]
									} else {
										tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC] = tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_CN]
									}
								} else {
									tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC] = tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_CN]
								}
							}
						}
						TT = TT+tt[2,j]
					} else {
						ttchr = numeric(0)
						for (p in 1:tt[2,j]) {
							ttchr = c(ttchr,trunc(tmpS[TT+p-1,GAP::.GapEnv$c_chr]))
						}
						ttchr = table(ttchr)
						p = 1
						for (p in 1:length(ttchr)) {
							kkchr = as.numeric(rownames(ttchr)[p])
							tt2j = ttchr[p]
							kk1 = -1
							if (trunc(tmpB[kk,GAP::.GapEnv$c_chr])==kkchr) {
								kk1 = kk
							}
							if (trunc(tmpB[kk+1,GAP::.GapEnv$c_chr])==kkchr) {
								kk1 = kk+1
							}
							if (kk1==-1) {
								
							} else {
								if (tt2j==1) {
									if (abs(tmpS[TT,GAP::.GapEnv$c_nC]-tmpB[kk1,GAP::.GapEnv$c_CN])<0.9) {
										tmpS[TT,GAP::.GapEnv$c_nC] = tmpB[kk1,GAP::.GapEnv$c_CN]
									} else {
										tmpS[TT,GAP::.GapEnv$c_nC] = tmpS[TT,GAP::.GapEnv$c_CN]
									}
								} else {
									if (max(abs(tmpS[TT:(TT+tt2j-1),GAP::.GapEnv$c_nC]-tmpB[kk1,GAP::.GapEnv$c_CN]))<0.9) {
										tmpS[TT:(TT+tt2j-1),GAP::.GapEnv$c_nC] = tmpB[kk1,GAP::.GapEnv$c_CN]
									} else {
										tmpS[TT:(TT+tt2j-1),GAP::.GapEnv$c_nC] = tmpS[TT:(TT+tt2j-1),GAP::.GapEnv$c_CN]
									}
								}
							}
							TT = TT+tt2j
						}
					}
				}
			}
		}
		tmpS = getGTwithRank(tmpS,Centr_CL,NETM)
		tmpS[,GAP::.GapEnv$c_CN] = tmpS[,GAP::.GapEnv$c_CN1]
		tmpS[,GAP::.GapEnv$c_BA] = tmpS[,GAP::.GapEnv$c_BA1]
		if ((dim(tmpS)[1]==2)&&(tmpS[1,GAP::.GapEnv$c_ind]==tmpS[2,GAP::.GapEnv$c_ind])) {
			dim_tmpS = 1
		} else {
			dim_tmpS = dim(tmpS)[1]
		}
		for (k in 1:dim_tmpS) {
			kk = which(tmp[,GAP::.GapEnv$c_ind]==tmpS[k,GAP::.GapEnv$c_ind])
			tmp[kk,] = tmpS[k,]
		}
	}
	return(invisible(tmp))
}
