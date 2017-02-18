'getCenterClusterGAP' <- function(CN, pp, qq, Delta)
{
	Centr_CL = Centr_CL_annot(CN)
	SR_perf = (-4)
	for (cn in 1:8) {
		SR_perf = c(SR_perf,log2(1.1*cn/2))
	}
	SR_perf = SR_perf*qq
	if (SR_perf[6]-SR_perf[5]<0.1) {
		for (cn in 7:9) {
			SR_perf[cn:9] = SR_perf[cn:9]-(SR_perf[cn]-SR_perf[cn-1])/2+(SR_perf[6]-SR_perf[5])/2
		}
	}
	DD = SR_perf[Delta[1]+1]-Delta[2]
	SR_perf = round((SR_perf-DD),4)
	for (cn in 1:9) {
		Centr_CL[5,which(Centr_CL[1,]==(cn-1))] = SR_perf[cn]
	}
	Centr_CL[4,1] = 0.5
	nbn = 1
	for (cn in 1:8) {
		max_F = rep(0,cn+1)
		if (cn%%2==0) {
			max_F[cn/2+1] = 0.5
		}
		for (nbc in trunc(cn/2+1):cn) {
			baf = ((1-pp)*nbc+pp*nbn)/((1-pp)*cn+2*pp)
			max_F[nbc+1] = round(baf,4)
		}
		cn_CL = 1
		while (Centr_CL[1,cn_CL]<cn) {
			cn_CL = cn_CL+1
		}
		for (nbc in 1:trunc(cn/2+1)) {
			Centr_CL[4,cn_CL+nbc-1] = max(0.5,min(1,max_F[trunc((cn+1)/2)+nbc]))
		}
	}
	for (nbc in 1:trunc(cn/2+1)) {
		Centr_CL[4,cn_CL+nbc-1] = max(0.5,min(1,max_F[trunc((cn+1)/2)+nbc]))
	}
	Centr_CL[5,1] = Centr_CL[5,2]-1.5*(Centr_CL[5,4]-Centr_CL[5,2])
	if (Centr_CL[4,4]<0.98) {
		for (cn_CL in 5:25) {
			if (Centr_CL[2,cn_CL]==Centr_CL[1,cn_CL]) {
				Centr_CL[4,cn_CL] = min(0.97,Centr_CL[4,cn_CL])
			}
		}
	}
	return(invisible(Centr_CL))
}
