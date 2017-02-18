'getGTwithRank' <- function(tmp, Centr_CL, NETM)
{
	Centr_CL = cbind(Centr_CL,Centr_CL[,which(Centr_CL[1,]-Centr_CL[2,]==0)])
	Centr_CL[2,26:34] = 0
	Centr_CL[3,26:34] = Centr_CL[1,26:34]+30
	Centr_CL[4,26:34] = 1
	Centr_CL = Centr_CL[,order(Centr_CL[3,],decreasing=FALSE)]
	tmp[,GAP::.GapEnv$c_nS] = 0
	k = 395
	for (k in 1:dim(tmp)[1]) {
		if (!is.na(tmp[k,GAP::.GapEnv$c_lrr]+tmp[k,GAP::.GapEnv$c_baf])) {
			if (trunc(tmp[k,GAP::.GapEnv$c_nC])==tmp[k,GAP::.GapEnv$c_nC]) {
				kk = c(tmp[k,GAP::.GapEnv$c_nC])
			} else {
				kk = c(trunc(tmp[k,GAP::.GapEnv$c_nC]),trunc(tmp[k,GAP::.GapEnv$c_nC])+1)
			}
			tt = which(kk<0)
			if (length(tt)>0) {
				kk = kk[-tt]
			}
			tt = which(kk>8)
			if (length(tt)>0) {
				kk = kk[-tt]
			}
		 	if (length(kk)==1) {
		 		kk = c(kk,kk)
		 	}
			kk = Centr_CL[3,which(Centr_CL[1,]%in%kk)]
			if (length(kk)==1) {
				kk = c(kk,kk)
			}
			tt = which(NETM[,3] %in% kk)
			if (length(tt)>1) {
				tt = NETM[tt,]
			} else {
				if (length(tt)==1) {
					tt = rbind(NETM[tt,],NETM[tt,])
				}
			}
			kk = which(NETM[,3]>=40)
			if (length(kk)>0) {
				tt = rbind(tt,NETM[kk,])
			}
			tt[,4] = sqrt((tt[,1]-tmp[k,GAP::.GapEnv$c_lrr])^2+4*(tt[,2]-tmp[k,GAP::.GapEnv$c_baf])^2)
			tt = tt[order(tt[,3],tt[,4],decreasing=FALSE),]
			hh = table(tt[,3])
			hh = c(0,hh)
			for (p in 2:length(hh)) {
				tt[(hh[p-1]+1),4] = sum(tt[(hh[p-1]+1):min((hh[p-1]+10),(hh[p-1]+hh[p])),4])/min(10,hh[p])
				if (1<hh[p]) {
					tt[(hh[p-1]+2):(hh[p-1]+hh[p]),3] = 0
				}
				hh[p] = hh[p-1]+hh[p]
			}
			ss = which(tt[,3]==0)
			if (length(ss)>0) {
				if (dim(tt)[1]-length(ss)>1) {
					tt = tt[-ss,]
				} else {
					tt = tt[-ss,];tt = rbind(tt,tt)
				}
			}
			tt = tt[order(tt[,4],decreasing=FALSE),]
			tt[,4] = signif(tt[,4],digits=1)
			if ((tt[1,3]>=40)&&(tmp[k,GAP::.GapEnv$c_len]>=50)) {
				tmp[k,GAP::.GapEnv$c_nS] = tt[1,3]
			}
			pp = which(tt[,3]>=40)
			if (length(pp)>0) {
				if (length(pp)<(dim(tt)[1]-2)) {
					tt = tt[-which(tt[,3]>=40),]
				} else {
					tt = tt[-which(tt[,3]>=40),];tt = rbind(tt,tt,tt)
				}
			}
			for (p in 1:min(3,dim(tt)[1])) {
				if (tt[p,3]>=30) {
					tmp[k,GAP::.GapEnv$c_CN1+3*(p-1)] = tt[p,3]-30
					tmp[k,GAP::.GapEnv$c_BA1+3*(p-1)] = 0
				} else {
					tmp[k,GAP::.GapEnv$c_CN1+3*(p-1)] = Centr_CL[1,tt[p,3]]
					tmp[k,GAP::.GapEnv$c_BA1+3*(p-1)] = Centr_CL[2,tt[p,3]]
				}
				tmp[k,GAP::.GapEnv$c_rnk1+3*(p-1)] = tt[p,4]
			}
		}
	}
	return(invisible(tmp))
}
