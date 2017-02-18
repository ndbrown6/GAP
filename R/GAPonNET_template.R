'GAPonNET_template' <- function(Centr_CL, NET_C, SF)
{
	NF = length(SF)
	step_SF = SF[2]-SF[1]
	NET_tmp = as.numeric(NET_C)
	for (k in 1:4) {
		NET_tmp = cbind(NET_tmp,NET_tmp)
	}
	NET_tmp = NET_tmp[,1:9]
	colnames(NET_tmp) = c("N_LRR","N_BAF","Step_LRR_down","Step_LRR_up","Step_BAF_le","Step_BAF_ri","Copy","BAF","Clust_Num")
	tt = as.numeric(rownames(NET_C))
	NET_tmp[,1] = rep(tt,dim(NET_C)[2])
	NET_tmp[,2] = rep(as.numeric(colnames(NET_C)),each=dim(NET_C)[1])
	NET_tmp[,3] = rep(tt-1-c(tt[1]-3,tt[1:(length(tt)-1)]),dim(NET_C)[2])
	NET_tmp[,4] = rep(c(tt[2:length(tt)],tt[length(tt)]+3)-tt-1,dim(NET_C)[2])
	tt = which(NET_tmp[,9]==0)
	if (length(tt)>0) {
		NET_tmp = NET_tmp[-tt,]
	}
	tt = which(NET_tmp[,7]<30)
	if (length(tt)>0) {
		NET_tmp[tt,7] = Centr_CL[1,NET_tmp[tt,7]]
	}
	tt = which(NET_tmp[,7]>=30)
	if (length(tt)>0) {
		NET_tmp[tt,7] = NET_tmp[tt,7]-30
	}
	tt = which(NET_tmp[,8]<30)
	if (length(tt)>0) {
		NET_tmp[tt,8] = Centr_CL[2,NET_tmp[tt,8]]
	}
	tt = which(NET_tmp[,8]>=30)
	if (length(tt)>0) {
		NET_tmp[tt,8] = NET_tmp[tt,8] = 0
	}
	NET_tmp = NET_tmp[order(NET_tmp[,1],NET_tmp[,2],decreasing=TRUE),]
	NET_tmp[,3:4] = round((NET_tmp[,3:4]-2)/2.5,0)
	NET_tmp[1:(dim(NET_tmp)[[1]]-1),5] = NET_tmp[1:(dim(NET_tmp)[[1]]-1),2]-NET_tmp[2:dim(NET_tmp)[[1]],2]-1
	tt = which(NET_tmp[,5]<0)
	if (length(tt)>0) {
		NET_tmp[tt,5] = (4/3)*(NET_tmp[tt,2]-1)
	}
	NET_tmp[dim(NET_tmp)[[1]],5] = NET_tmp[dim(NET_tmp)[[1]],2]-1
	NET_tmp[2:dim(NET_tmp)[[1]],6] = NET_tmp[1:(dim(NET_tmp)[[1]]-1),2]-NET_tmp[2:dim(NET_tmp)[[1]],2]-1
	for (k in 1:dim(NET_tmp)[[1]]) {
		if ((NET_tmp[k,8]==0)&&(NET_tmp[k,5]>=7)) {
			NET_tmp[k,5] = 16
		}
		if ((NET_tmp[k,8]==NET_tmp[k,7])&&(abs(NET_tmp[k,6])<16)) {
			NET_tmp[k,6] = (97-NET_tmp[k,2])*4
		}
	}
	NET_tmp[,5:6] = round(NET_tmp[,5:6]/4,0)
	tt = which(NET_tmp[,6]<0)
	if (length(tt)>0) {
		NET_tmp[tt,6] = NF-NET_tmp[tt,2]
	}
	NET_tmp[1,6] = max(1,NF-NET_tmp[1,2])
	tt = which(NET_tmp[,5]>4)
	if (length(tt)>0) {
		NET_tmp[tt,5] = 4
	}
	for (k in 1:dim(NET_tmp)[1]) {
		if (NET_tmp[k,6]>4) {
			if (NET_tmp[k,2]==1) {
				NET_tmp[k,6] = min(7,NET_tmp[k,6])
			} else {
				NET_tmp[k,6] = 4
			}
		}
	}
	return(invisible(NET_tmp))
}
