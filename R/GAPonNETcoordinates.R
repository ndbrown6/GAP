'GAPonNETcoordinates' <- function(Centr_CL, SF, SR)
{
	NF = length(SF)
	NR = length(SR)
	step_SF = SF[2]-SF[1]
	step_SR = SR[2]-SR[1]
	CN = max(Centr_CL[1,])
	NET_C = matrix(0,CN+1+1,NF)
	rownames(NET_C) = seq(0,CN+1,1)
	colnames(NET_C) = seq(1,NF,1)
	NET_C[,(dim(NET_C)[[2]]-3)] = 1
	cn_CL = 1
	for (cn_CL in 1:dim(Centr_CL)[[2]]) {
		max_R = min(abs(Centr_CL[5,cn_CL]-SR))
		if (max_R<step_SR) {
			max_R = which.min(abs(Centr_CL[5,cn_CL]-SR))
			if (Centr_CL[5,cn_CL]>=SR[max_R]) {
				max_R = min(NR,max_R+1)
			}
		} else {
			max_R = round(max_R/step_SR)
			if (Centr_CL[5,cn_CL]>SR[NR]) {
				max_R = NR+max_R
			}
			if (Centr_CL[5,cn_CL]<SR[1]) {
				max_R = (-max_R)
			}
		}
		rownames(NET_C)[Centr_CL[1,cn_CL]+1] = max_R
		max_F = which.min(abs(SF-Centr_CL[4,cn_CL]))
		NET_C[Centr_CL[1,cn_CL]+1,max_F] = cn_CL
		NET_C[Centr_CL[1,cn_CL]+1,(dim(NET_C)[[2]]-3)] = Centr_CL[1,cn_CL]+30
		NET_C[dim(NET_C)[[1]],max_F] = NET_C[dim(NET_C)[[1]],max_F]+1
	}
	tt = which(NET_C[dim(NET_C)[[1]],]==0)
	if (length(tt)>0) {
		NET_C = NET_C[,-tt]
	}
	NET_C = NET_C[-dim(NET_C)[[1]],]
	return(invisible(NET_C))
}
