'updateMeans_FinalPolish' <- function(tmp, datA, datL, germHomozyg.mBAF.thr=0.97)
{
	tt = which(is.na(datA))
	if (length(tt)!=0) {
		datA[tt] = 1
	}
	tt = which((datA>1))
	if (length(tt)!=0) {
		datA[tt] = 1
	}
	for (k in 1:dim(tmp)[1]) {
		if (tmp[k,GAP::.GapEnv$c_if]>tmp[k,GAP::.GapEnv$c_is]+10) {
			tmp[k,GAP::.GapEnv$c_lrr] = median(datL[tmp[k,GAP::.GapEnv$c_is]:tmp[k,GAP::.GapEnv$c_if]],na.rm=TRUE)
		} else {
			if (tmp[k,GAP::.GapEnv$c_if]>tmp[k,GAP::.GapEnv$c_is]) {
				tmp[k,GAP::.GapEnv$c_lrr] = mean(datL[tmp[k,GAP::.GapEnv$c_is]:tmp[k,GAP::.GapEnv$c_if]],na.rm=TRUE)
			} else {
				tmp[k,GAP::.GapEnv$c_lrr] = datL[tmp[k,GAP::.GapEnv$c_is]]
			}
		}
	}
	ad = datA
	tt = which(datA>germHomozyg.mBAF.thr)
	if (length(tt)!=0) {
		ad[tt] = NA
	}
	for (k in 1:dim(tmp)[1]) {
		adv = ad[tmp[k,GAP::.GapEnv$c_is]:tmp[k,GAP::.GapEnv$c_if]]
	  	tt = which(!is.na(adv))
	  	if (length(tt)>length(adv)/8) {
			adv = adv[tt]
			if (length(adv)>30) {
		  		tmp[k,GAP::.GapEnv$c_baf] = median(adv,na.rm=TRUE)
		  		if (tmp[k,GAP::.GapEnv$c_baf]<0.6) {
					hh = hist(adv, breaks=seq(0.5,1,0.02),right=FALSE,plot=FALSE)
					if ((sum(hh$counts[1:3])/sum(hh$counts)>0.7)&&((hh$counts[1]+hh$counts[2])>2*hh$counts[3]*1.2)) {
						if (hh$counts[1]>hh$counts[2]) {
							tmp[k,GAP::.GapEnv$c_baf] = 0.5
						} else {
							tmp[k,GAP::.GapEnv$c_baf] = 0.52
						}
					} else {
						dd = density(adv,n=20)
						tmp[k,GAP::.GapEnv$c_baf] = dd$x[which.max(dd$y)]
					}
		  		}
			} else {
				if (length(adv)>0) {
					hh = hist(adv, breaks=seq(0.5,1,0.05),plot=FALSE)
					tmp[k,GAP::.GapEnv$c_baf] = hh$mids[which.max(hh$counts)]
				} else {
					tmp[k,GAP::.GapEnv$c_baf] = NA
				}
			}
		} else {
	  		tt = which(is.na(adv))
			adv = datA[tmp[k,GAP::.GapEnv$c_is]:tmp[k,GAP::.GapEnv$c_if]]
			if (length(tt)>10) {
				adv = adv[tt]
				if (length(adv)>30) {
					dd = density(adv,n=20)
					tmp[k,GAP::.GapEnv$c_baf] = dd$x[which.max(dd$y)]
				} else {
					hh = hist(adv, breaks=seq(0.5,1,0.05),right=FALSE,plot=FALSE)
					tmp[k,GAP::.GapEnv$c_baf] = hh$mids[which.max(hh$counts)]
				}
			} else {
				if (length(tt)>0) {
					adv = adv[tt]
					hh = hist(adv, breaks=seq(0.5,1,0.1),plot=FALSE)
					tmp[k,GAP::.GapEnv$c_baf] = hh$mids[which.max(hh$counts)]
				} else {
					tmp[k,GAP::.GapEnv$c_baf] = NA
				}
			}
		}
	}
	tmp = round(tmp,4)
	return(invisible(tmp))
}
