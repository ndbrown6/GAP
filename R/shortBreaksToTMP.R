'shortBreaksToTMP' <- function(tmp, CHR, Centr_CL, NETM)
{
	ss = data.matrix(table(tmp[,GAP::.GapEnv$c_chr],sign(tmp[,GAP::.GapEnv$c_conf])))
	tt = as.numeric(which(ss[,2]==0))
	if (length(tt)>0) {
		ss = as.numeric(rownames(ss))[tt]
		k = 1
		for (k in 1:length(ss)) {
			kk = which(tmp[,GAP::.GapEnv$c_chr]==ss[k])
			dd = max(tmp[1:kk[1],GAP::.GapEnv$c_conf])
			for (p in 1:length(kk)){
				tmp[kk[p],GAP::.GapEnv$c_conf] = dd+p
			}
			if (kk[length(kk)]<dim(tmp)[1]){ 
				tmp[(kk[length(kk)]+1):dim(tmp)[1],GAP::.GapEnv$c_conf] = tmp[(kk[length(kk)]+1):dim(tmp)[1],GAP::.GapEnv$c_conf]+length(kk)
			}
		}
	}
	tt = matrix(0,2,2)
	t_pos = 0
	t_chr = 1
	k_neg = 0
	k = 64
	for (k in 1:dim(tmp)[1]) {
		if (tmp[k,GAP::.GapEnv$c_chr]==t_chr) {
			if (tmp[k,GAP::.GapEnv$c_conf]<0) {
				k_neg = k_neg + 1
			} else {
				if (k_neg>0) {
					if (is.na(t_pos)) {
						tt = cbind(tt,c(tmp[k,GAP::.GapEnv$c_conf],k_neg));t_pos = tmp[k,GAP::.GapEnv$c_conf]
					} else {
						if (tmp[k,GAP::.GapEnv$c_conf]==t_pos) {
							tt = cbind(tt,c(t_pos,k_neg))
						} else {
							tt = cbind(tt,c((t_pos+tmp[k,GAP::.GapEnv$c_conf])/2,k_neg))
							t_pos = tmp[k,GAP::.GapEnv$c_conf]
						}
					}
					k_neg = 0
				} else {
					t_pos = tmp[k,GAP::.GapEnv$c_conf]
				}
			}
		} else {
			if (k_neg>0) {
				tt = cbind(tt,c(t_pos,k_neg))
				k_neg = 0
			}
			t_pos = NA
			t_chr = tmp[k,GAP::.GapEnv$c_chr]
			if (tmp[k,GAP::.GapEnv$c_conf]<0) {
				k_neg = k_neg+1
			} else {
				t_pos = tmp[k,GAP::.GapEnv$c_conf]
			}
		}
	}
	Rest0 = 0
	if (dim(tt)[2]>4) {
		tt = tt[,-c(1,2)]
		Rest0 = 2
	} else {
		tt = cbind(tt[,-c(1,2)],tt[,-c(1,2)])
		Rest0 = 1
	}
	if (Rest0>0) {
		tt = data.matrix(tt)
		colnames(tt) = seq(1,dim(tt)[2])
		kk = which(tmp[,GAP::.GapEnv$c_conf]>0)
		if (length(kk)>0) {
			tmpS = tmp[-kk,]
			if ((dim(tmp)[1]-length(kk)==1)||(Rest0==1)) {
				tmpS = rbind(tmpS,tmpS)
			}
			tmpB = tmp[kk,]
		}
		tmpB[1,GAP::.GapEnv$c_is] = 1
		for (k in 2:dim(tmpB)[1]) {
			if (tmpB[k-1,GAP::.GapEnv$c_chr]!=tmpB[k,GAP::.GapEnv$c_chr]) {
				kk = which(as.numeric(names(CHR))==tmpB[k-1,GAP::.GapEnv$c_chr])
				tmpB[k-1,GAP::.GapEnv$c_if] = CHR[kk]
				tmpB[k,GAP::.GapEnv$c_is] = CHR[kk]+1
			}
		}
		tmpB = shrinkReprTMP(tmpB)
		tt1 = which(tt[1,]==0.5)
		if (length(tt1)>0) {
			tt[1,tt1] = 1
		}
		tmpS[,GAP::.GapEnv$c_nS] = 1
		tt_end = length(which(tt[1,]<=dim(tmpB)[1]))
		TT = 1
		for (j in 1:tt_end) {
			if (tt[2,j]>1) {
				tt1 = TT:(TT+tt[2,j]-1)
			} else {
				tt1 = TT
			}
			tmpS[tt1,GAP::.GapEnv$c_bafH] = tt[1,j]
			if (trunc(tt[1,j])==tt[1,j]) {
				tmpS[tt1,GAP::.GapEnv$c_CN] = tmpB[tt[1,j],GAP::.GapEnv$c_CN]
				tmpS[tt1,GAP::.GapEnv$c_BA] = tmpB[tt[1,j],GAP::.GapEnv$c_BA]
			}
			TT = TT+abs(tt[2,j])
		}
		tmpS[,GAP::.GapEnv$c_nS] = 1
		tt_end = length(which(tt[1,]<=dim(tmpB)[1]))
		if ((Rest0==1)&&(tt_end==2)) {
			tt_end = 1
		}
		TT = 1
		j = 26
		for (j in 1:tt_end) {
			if (trunc(tt[1,j])==tt[1,j]) {
				if (tt[2,j]==1) {
					if (tmpS[TT,GAP::.GapEnv$c_len]<10) {
						tmpS[TT,GAP::.GapEnv$c_nS] = 0
						tt[2,j] = (-tt[2,j])
					} else {
						if (abs(tmpS[TT,GAP::.GapEnv$c_nC]-tmpS[TT,GAP::.GapEnv$c_CN])<1) {
							if (tmpS[TT,GAP::.GapEnv$c_len]<20) {
								tmpS[TT,GAP::.GapEnv$c_nS] = 0
								tt[2,j] = -tt[2,j]
							} else {
								if (((tmpS[TT,GAP::.GapEnv$c_CN]==tmpS[TT,GAP::.GapEnv$c_CN1])&&((tmpS[TT,GAP::.GapEnv$c_BA1]==0)||(tmpS[TT,GAP::.GapEnv$c_BA1]==tmpS[TT,GAP::.GapEnv$c_BA])||(tmpS[TT,GAP::.GapEnv$c_BA]==0)))
								   ||((tmpS[TT,GAP::.GapEnv$c_CN]==tmpS[TT,GAP::.GapEnv$c_CN1+3])&&((tmpS[TT,GAP::.GapEnv$c_BA1+3]==0)||(tmpS[TT,GAP::.GapEnv$c_BA1+3]==tmpS[TT,GAP::.GapEnv$c_BA])||(tmpS[TT,GAP::.GapEnv$c_BA1+3]==0)))) {
										tmpS[TT,GAP::.GapEnv$c_nS] = 0
										tt[2,j] = -tt[2,j]
								}
							}
	  					} else {
	  						if (min(tmpS[TT,GAP::.GapEnv$c_nC],tmpS[TT,GAP::.GapEnv$c_CN])>5) {
	  							tmpS[TT,GAP::.GapEnv$c_nS] = 0
	  							tt[2,j] = -tt[2,j]
	  						}
						}
					}
				}
			}
			TT = TT + abs(tt[2,j])
		}
		Rest1 = 0
		kk = which(tt[2,]<0)
		if (length(kk)>0) {
			if (dim(tt)[2]>length(kk)+1) {
				tt = tt[,-kk]
				Rest1 = 2
			} else {
				if (dim(tt)[2]==length(kk)+1) {
					if (Rest0==2) {
						tt = cbind(tt[,-kk],tt[,-kk])
						Rest1 = 1
					}
				}
			}
		}
		kk = which(tmpS[,GAP::.GapEnv$c_nS]==0)
		if (length(kk)>0) {
			if (dim(tmpS)[1]>length(kk)+1) {
				tmpS = tmpS[-kk,]
			} else {
				if (dim(tmpS)[1]==length(kk)+1) {
					tmpS = rbind(tmpS[-kk,],tmpS[-kk,])
				}
			}
		}
		if (Rest1>0) {
			tt_end = length(which(tt[1,]<=dim(tmpB)[1]))
			if ((Rest1==1)&&(tt_end==2)) {
				tt_end = 1
			}
			TT = 1
			j = 26
			for (j in 1:tt_end) {
		  		if (trunc(tt[1,j])!=tt[1,j]) {
					kk = trunc(tt[1,j])
					if (tt[2,j]==1) {
	  					if (tmpS[TT,GAP::.GapEnv$c_len]<10) {
							if (!is.na(tmpS[TT,GAP::.GapEnv$c_lrr])) {
								if (kk<dim(tmpB)[1]) {
									pp = which.min(c(abs(tmpS[TT,GAP::.GapEnv$c_lrr]-tmpB[kk,GAP::.GapEnv$c_lrr]),abs(tmpS[TT,GAP::.GapEnv$c_lrr]-tmpB[kk+1,GAP::.GapEnv$c_lrr])))
								} else {
									pp = 1
								}
								tmpS[TT,GAP::.GapEnv$c_bafH] = -(kk+pp-1)
								tmpS[TT,GAP::.GapEnv$c_nS] = 0
								tt[2,j] = -tt[2,j]
							} else {
								tmpS[TT,GAP::.GapEnv$c_bafH] = (-kk)
								tmpS[TT,GAP::.GapEnv$c_nS] = 0
								tt[2,j] = (-tt[2,j])
							}
						} else {
	  						if (kk<dim(tmpB)[1]) {
								if ((abs(tmpS[TT,GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN])<1)&&(abs(tmpS[TT,GAP::.GapEnv$c_nC]-tmpB[kk+1,GAP::.GapEnv$c_CN])>=1)) {
									tmpS[TT,GAP::.GapEnv$c_bafH] = -kk
									tmpS[TT,GAP::.GapEnv$c_nS] = -1
								}
								if ((abs(tmpS[TT,GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN])>=1)&&(abs(tmpS[TT,GAP::.GapEnv$c_nC]-tmpB[kk+1,GAP::.GapEnv$c_CN])<1)) {
									tmpS[TT,GAP::.GapEnv$c_bafH] = (-kk-1)
									tmpS[TT,GAP::.GapEnv$c_nS] = -1
								}
								if ((abs(tmpS[TT,GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN])<1)&&(abs(tmpS[TT,GAP::.GapEnv$c_nC]-tmpB[kk+1,GAP::.GapEnv$c_CN])<1)) {
									pp = which.min(c(tmpS[TT,GAP::.GapEnv$c_lrr],tmpS[TT,GAP::.GapEnv$c_lrr])-tmpB[kk:(kk+1),GAP::.GapEnv$c_lrr])
									tmpS[TT,GAP::.GapEnv$c_bafH] = (-(kk+pp-1))
									tmpS[TT,GAP::.GapEnv$c_nS] = -1
								}
								if ((tmpS[TT,GAP::.GapEnv$c_bafH]<0)&&(abs(tmpS[TT,GAP::.GapEnv$c_bafH])==trunc(abs(tmpS[TT,GAP::.GapEnv$c_bafH])))&&(tmpS[TT,GAP::.GapEnv$c_len]<20)) {
									tmpS[TT,GAP::.GapEnv$c_nS] = 0
									tt[2,j] = -tt[2,j]
								}
							} else {
								if (abs(tmpS[TT,GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN])<1) {
									tmpS[TT,GAP::.GapEnv$c_bafH] = -kk
									tmpS[TT,GAP::.GapEnv$c_nS] = 0
									tt[2,j] = -tt[2,j]
								}
							}
						}
						kk = abs(tmpS[TT,GAP::.GapEnv$c_bafH])
						if ((trunc(kk)==kk)&&(kk!=0)) {
							tmpS[TT,GAP::.GapEnv$c_CN] = tmpB[kk,GAP::.GapEnv$c_CN]
							tmpS[TT,GAP::.GapEnv$c_BA] = tmpB[kk,GAP::.GapEnv$c_BA]
						}
					}
				}
				TT = TT+abs(tt[2,j])
			}
			k = 66
			if ((dim(tmpS)[1]==2)&&(tmpS[1,GAP::.GapEnv$c_ind]==tmpS[2,GAP::.GapEnv$c_ind])) {
				dim_tmpS = 1
			} else {
				dim_tmpS = dim(tmpS)[1]
			}
			for (k in 1:dim_tmpS) {
				if ((tmpS[k,GAP::.GapEnv$c_bafH]<0)&&(tmpS[k,GAP::.GapEnv$c_nS]==0)) {
					tmpS[k,GAP::.GapEnv$c_bafH] = abs(tmpS[k,GAP::.GapEnv$c_bafH])
					if (tmpS[k,GAP::.GapEnv$c_is]==tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_if]+1) {
						tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_if] = tmpS[k,GAP::.GapEnv$c_if]
						tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len] = tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len]+tmpS[k,GAP::.GapEnv$c_len]
					} else {
						if (tmpS[k,GAP::.GapEnv$c_if]==tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_is]-1) {
							tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_is] = tmpS[k,GAP::.GapEnv$c_is]
							tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len] = tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len]+tmpS[k,GAP::.GapEnv$c_len]
						} else {
							tmpS[k,GAP::.GapEnv$c_nS] = 1
							print("An unknown error occurred...!")
						}
					}
				}
			}
			Rest2 = 0
			kk = which(tt[2,]<0)
			if (length(kk)>0) {
				if (dim(tt)[2]>length(kk)+1) {
					tt = tt[,-kk]
					Rest2 = 2
				} else {
					if (dim(tt)[2]==length(kk)+1) {
						if (Rest1==2) {
							tt = cbind(tt[,-kk],tt[,-kk])
							Rest2 = 1
						}
					}
				}
			} else {
				Rest2 = Rest1
			}
			kk = which(tmpS[,GAP::.GapEnv$c_nS]==0)
			if (length(kk)>0) {
				if ((dim(tmpS)[1]==2)&&(tmpS[1,GAP::.GapEnv$c_ind]==tmpS[2,GAP::.GapEnv$c_ind])) {
					Rest2 = 0
				} else {
					if (dim(tmpS)[1]>length(kk)+1) {
						tmpS = tmpS[-kk,]
					} else {
						if (dim(tmpS)[1]==length(kk)+1) {
							tmpS = rbind(tmpS[-kk,],tmpS[-kk,])
						}
					}
				}
			}
			if (Rest2>0) {
				TT = 1
				tt_end = length(which(tt[1,]<=dim(tmpB)[1]))
				if ((Rest2==1)&&(tt_end==2)) {
					tt_end = 1
				}
				for (j in 1:tt_end) {
					if (trunc(tt[1,j])==tt[1,j]) {
						if (tt[2,j]>1) {
							if ((max(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_len])<10)&&(max(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_CN]),na.rm=TRUE)<1)) {
								tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nS] = 0
							} else {
								tt1 = which(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_len]>=10)
 								if (length(tt1)>0) {
 									tt1 = TT+tt1-1
									if ((max(abs(tmpS[tt1,GAP::.GapEnv$c_nC]-tmpS[tt1,GAP::.GapEnv$c_CN]))<1)&&(max(tmpS[tt1,GAP::.GapEnv$c_len])<20)) {
										tmpS[tt1,GAP::.GapEnv$c_nS] = 0
									} else {
										p = 1
										for (p in 1:length(tt1)) {
											if (tmpS[tt1[p],GAP::.GapEnv$c_len]>20) {
												if (((tmpS[tt1[p],GAP::.GapEnv$c_CN]==tmpS[tt1[p],GAP::.GapEnv$c_CN1])&&((tmpS[tt1[p],GAP::.GapEnv$c_BA1]==0)||(tmpS[tt1[p],GAP::.GapEnv$c_BA1]==tmpS[tt1[p],GAP::.GapEnv$c_BA])||(tmpS[tt1[p],GAP::.GapEnv$c_BA]==0)))
											   	   ||((tmpS[tt1[p],GAP::.GapEnv$c_CN]==tmpS[tt1[p],GAP::.GapEnv$c_CN1+3])&&((tmpS[tt1[p],GAP::.GapEnv$c_BA1+3]==0)||(tmpS[tt1[p],GAP::.GapEnv$c_BA1+3]==tmpS[tt1[p],GAP::.GapEnv$c_BA])||(tmpS[tt1[p],GAP::.GapEnv$c_BA1+3]==0)))) {
													tmpS[tt1[p],GAP::.GapEnv$c_nS] = 0
												}
											} else {
												if ((tmpS[tt1[p],GAP::.GapEnv$c_CN]==tmpS[tt1[p],GAP::.GapEnv$c_CN1])||(tmpS[tt1[p],GAP::.GapEnv$c_CN]==tmpS[tt1[p],GAP::.GapEnv$c_CN1+3])) {
													tmpS[tt1[p],GAP::.GapEnv$c_nS] = 0
												}
											}
										}
									}
								} else {
									tt1 = which(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_CN])>=1)
									if (length(tt1)<=1) {
										tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nS] = 0
									}
								}
							}
						}
  					} else {
  						kk = trunc(tt[1,j])
						if ((tt[2,j]>1)&&(kk<dim(tmpB)[1])) {
							tt1 = which(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_len]>=10)
							if (length(tt1)==0) {
								if ((max(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN]),na.rm=T)<1)&&(min(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpB[kk+1,GAP::.GapEnv$c_CN]),na.rm=TRUE)>=1)) {
									tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_bafH] = (-kk);tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nS] = 0
								}
								if ((max(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpB[kk+1,GAP::.GapEnv$c_CN]),na.rm=TRUE)<1)&&(min(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN]),na.rm=TRUE)>=1)) {
									tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_bafH] = (-kk-1);tmpB[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nS] = 0
								}
							} else {
								tt1 = TT + tt1-1
								if ((max(abs(tmpS[tt1,GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN]))<1)&&(min(abs(tmpS[tt1,GAP::.GapEnv$c_nC]-tmpB[kk+1,GAP::.GapEnv$c_CN]))>=1)) {
									tmpS[tt1,GAP::.GapEnv$c_bafH] = (-kk);tmpS[tt1,GAP::.GapEnv$c_nS] = 0
								}
								if ((min(abs(tmpS[tt1,GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN]))>=1)&&(max(abs(tmpS[tt1,GAP::.GapEnv$c_nC]-tmpB[kk+1,GAP::.GapEnv$c_CN]))<1)) {
									tmpS[tt1,GAP::.GapEnv$c_bafH] = (-kk-1);tmpS[tt1,GAP::.GapEnv$c_nS] = 0
								}
								if (sum(tmpS[tt1,GAP::.GapEnv$c_nS])!=0) {
									t_adH = tmpS[tt1[1],GAP::.GapEnv$c_bafH]
									p = 0
									while (((p+1)<=length(tt1))&&(abs(tmpS[tt1[p+1],GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN])<1)) {
										tmpS[tt1[p+1],GAP::.GapEnv$c_bafH] = (-kk)
										tmpS[tt1[p+1],GAP::.GapEnv$c_nS] = -1
										p = p + 1
									}
									if (p!=0) {
										p_up = seq(1,p)
									} else {
										p_up = numeric(0)
									}
									p = 0
									while (((p+1)<=length(tt1))&&(abs(tmpS[tt1[length(tt1)-p],GAP::.GapEnv$c_nC]-tmpB[kk+1,GAP::.GapEnv$c_CN])<1)) {
										tmpS[tt1[length(tt1)-p],GAP::.GapEnv$c_bafH] = (-kk-1)
										tmpS[tt1[length(tt1)-p],GAP::.GapEnv$c_nS] = -1
										p = p + 1
									}
									if (p!=0) {
										p_down = seq(length(tt1)-p+1,length(tt1))
									} else {
										p_down = numeric(0)
									}
									p_unkn = intersect(p_up,p_down)
									if (length(p_unkn)>0) {
										tmpS[tt1[p_unkn],GAP::.GapEnv$c_bafH] = (t_adH)
										tmpS[tt1[p_unkn],GAP::.GapEnv$c_nS] = 1
									} else {
										p_unkn = setdiff(seq(1,length(tt1)),union(p_up,p_down))
										if (length(p_unkn)>0) {
											tmpS[tt1[p_unkn],GAP::.GapEnv$c_bafH] = (t_adH)
											tmpS[tt1[p_unkn],GAP::.GapEnv$c_nS] = 1
										}
									}
								}	
							}
						} else {
							if (tt[2,j]>1) {
								tt1 = which(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_len]>=10)
								if (length(tt1)==0) {
									if (max(abs(tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nC]-tmpS[kk,GAP::.GapEnv$c_CN]),na.rm=TRUE)<1) {
										tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_bafH] = (-kk)
										tmpS[TT:(TT+tt[2,j]-1),GAP::.GapEnv$c_nS] = 0
									}
								} else {
									if (max(abs(tmpS[tt1,GAP::.GapEnv$c_nC]-tmpB[kk,GAP::.GapEnv$c_CN]))<1) {
										tmpS[tt1,GAP::.GapEnv$c_bafH] = (-kk)
										tmpS[tt1,GAP::.GapEnv$c_nS] = 0
									}
				  				}
							}
						}
					}
					TT = TT+abs(tt[2,j])
				}
				if ((dim(tmpS)[1]==2)&&(tmpS[1,GAP::.GapEnv$c_ind]==tmpS[2,GAP::.GapEnv$c_ind])) {
					dim_tmpS = 1
				} else{
					dim_tmpS = dim(tmpS)[1]
				}
				for (k in 1:dim_tmpS) {
					if ((tmpS[k,GAP::.GapEnv$c_bafH]<0)&&(tmpS[k,GAP::.GapEnv$c_nS]==0)) {
						tmpS[k,GAP::.GapEnv$c_bafH] = abs(tmpS[k,GAP::.GapEnv$c_bafH])
						if (tmpS[k,GAP::.GapEnv$c_is]>=tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_if]+1) {
							tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_if] = tmpS[k,GAP::.GapEnv$c_if]
							tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len] = tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len]+tmpS[k,GAP::.GapEnv$c_len]
						} else {
							if (tmpS[k,GAP::.GapEnv$c_if]<=tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_is]-1) {
								tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_is] = tmpS[k,GAP::.GapEnv$c_is]
								tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len] = tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len]+tmpS[k,GAP::.GapEnv$c_len]
							} else {
								if ((tmpS[k,GAP::.GapEnv$c_is]>=tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_is])&&(tmpS[k,GAP::.GapEnv$c_if]<=tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_if])) {
									tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len] = tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len]+tmpS[k,GAP::.GapEnv$c_len]
								} else {
									print("An unknown error occurred...!")
								}
							}
						}
					}
				}
				Rest3 = 2
				kk = which(tmpS[,GAP::.GapEnv$c_nS]==0)
				if (length(kk)>0) {
					if ((dim(tmpS)[1]==2)&&(tmpS[1,GAP::.GapEnv$c_ind]==tmpS[2,GAP::.GapEnv$c_ind])) {
						Rest3 = 0
					} else {
						if (dim(tmpS)[1]>length(kk)+1) {
							tmpS = tmpS[-kk,]
						} else {
							if (dim(tmpS)[1]==length(kk)+1) {
								tmpS = rbind(tmpS[-kk,],tmpS[-kk,])
								Rest3 = 1
							}
						}
					}
				} else {
					if ((dim(tmpS)[1]==2)&&(tmpS[1,GAP::.GapEnv$c_ind]==tmpS[2,GAP::.GapEnv$c_ind])) {
						Rest3 = 1
					}
				}
				if (Rest3>0) {
					tt = which(tmpS[,GAP::.GapEnv$c_nS]<0)
					if (length(tt)>0) {
						tmpS[tt,GAP::.GapEnv$c_nC] = tmpB[abs(tmpS[tt,GAP::.GapEnv$c_bafH]),GAP::.GapEnv$c_CN]
						if (length(tt)>1) {
							tmpO = tmpS[tt,]
						} else {
							tmpO = rbind(tmpS[tt,],tmpS[tt,])
						}
						tmpO = getGTwithRank(tmpO,Centr_CL,NETM)
						if (length(tt)>1) {
							tmpS[tt,GAP::.GapEnv$c_CN1:(GAP::.GapEnv$c_CN1+8)] = tmpO[,GAP::.GapEnv$c_CN1:(GAP::.GapEnv$c_CN1+8)]
							tmpS[tt,GAP::.GapEnv$c_CN] = tmpO[,GAP::.GapEnv$c_CN1]
							tmpS[tt,GAP::.GapEnv$c_BA] = tmpO[,GAP::.GapEnv$c_BA1]
						} else {
							tmpS[tt,GAP::.GapEnv$c_CN1:(GAP::.GapEnv$c_CN1+8)] = tmpO[1,GAP::.GapEnv$c_CN1:(GAP::.GapEnv$c_CN1+8)]
							tmpS[tt,GAP::.GapEnv$c_CN] = tmpO[1,GAP::.GapEnv$c_CN1]
							tmpS[tt,GAP::.GapEnv$c_BA] = tmpO[1,GAP::.GapEnv$c_BA1]
						}
						for (p in 1:length(tt)) {
							if (((tmpS[tt[p],GAP::.GapEnv$c_CN]==tmpS[tt[p],GAP::.GapEnv$c_CN1])&&((tmpS[tt[p],GAP::.GapEnv$c_BA1]==0)||(tmpS[tt[p],GAP::.GapEnv$c_BA1]==tmpS[tt[p],GAP::.GapEnv$c_BA])||(tmpS[tt[p],GAP::.GapEnv$c_BA]==0)))
							   ||((tmpS[tt[p],GAP::.GapEnv$c_CN]==tmpS[tt[p],GAP::.GapEnv$c_CN1+3])&&((tmpS[tt[p],GAP::.GapEnv$c_BA1+3]==0)||(tmpS[tt[p],GAP::.GapEnv$c_BA1+3]==tmpS[tt[p],GAP::.GapEnv$c_BA])||(tmpS[tt[p],GAP::.GapEnv$c_BA1+3]==0)))) {
								tmpS[tt[p],GAP::.GapEnv$c_nS] = 0
							}
						}
						for (k in 1:dim(tmpS)[1]) {
							if ((tmpS[k,GAP::.GapEnv$c_bafH]<0)&&(tmpS[k,GAP::.GapEnv$c_nS]==0)) {
								tmpS[k,GAP::.GapEnv$c_bafH] = abs(tmpS[k,GAP::.GapEnv$c_bafH])
								if (tmpS[k,GAP::.GapEnv$c_is]>=tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_if]+1) {
									tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_if] = tmpS[k,GAP::.GapEnv$c_if]
									tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len] = tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len]+tmpS[k,GAP::.GapEnv$c_len]
								} else {
									if (tmpS[k,GAP::.GapEnv$c_if]<=tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_is]-1) {
										tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_is] = tmpS[k,GAP::.GapEnv$c_is]
										tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len] = tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len]+tmpS[k,GAP::.GapEnv$c_len]
									} else {
										if ((tmpS[k,GAP::.GapEnv$c_is]>=tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_is])&&(tmpS[k,GAP::.GapEnv$c_if]<=tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_if])) {
											tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len] = tmpB[tmpS[k,GAP::.GapEnv$c_bafH],GAP::.GapEnv$c_len]+tmpS[k,GAP::.GapEnv$c_len]
										} else {
											print("An unknown error occurred...!")
										}
									}
								}
							}
						}
						kk = which(tmpS[,GAP::.GapEnv$c_nS]==0)
						if (length(kk)>0) {
							tmpS = tmpS[-kk,]
						}
					}
				}
			}
		}
		tmp[,GAP::.GapEnv$c_conf] = 0
		for (k in 1:dim(tmpB)[1]) {
			tt = which(tmp[,GAP::.GapEnv$c_is]>=tmpB[k,GAP::.GapEnv$c_is])
			if (length(tt)>0) {
				kk = which(tmp[tt,GAP::.GapEnv$c_if]<=tmpB[k,GAP::.GapEnv$c_if])
				if (length(kk)>0) {
					tmp[tt[kk],GAP::.GapEnv$c_CN] = tmpB[k,GAP::.GapEnv$c_CN]
					tmp[tt[kk],GAP::.GapEnv$c_BA] = tmpB[k,GAP::.GapEnv$c_BA]
					tmp[tt[kk],GAP::.GapEnv$c_conf] = k
				}
			}
		}
		if ((dim(tmpS)[1]==2)&&(tmpS[1,GAP::.GapEnv$c_ind]==tmpS[2,GAP::.GapEnv$c_ind])) {
			dim_tmpS = 1
		} else {
			dim_tmpS = dim(tmpS)[1]
		}
		for (k in 1:dim_tmpS) {
			tt = which(tmp[,GAP::.GapEnv$c_ind]==tmpS[k,GAP::.GapEnv$c_ind])
			if (length(tt)==1) {
				tmp[tt,GAP::.GapEnv$c_conf] = (-k)
				tmp[tt,GAP::.GapEnv$c_CN] = tmpS[k,GAP::.GapEnv$c_CN1]
				tmp[tt,GAP::.GapEnv$c_BA] = tmpS[k,GAP::.GapEnv$c_BA1]
			}
		}
	}
	return(invisible(tmp))
}
