'fittingGAP' <- function(NET)
{
	SR_perf = c(-1,0,0.58,1,1.32,2)
	pp_N = seq(0,0.6,0.05)
	qq_L = sort(c(seq(3,8,1)/32,seq(17,32,1)/64),decreasing=TRUE)
	
	NF = dim(NET)[2]
	NR = dim(NET)[1]
	SF = as.numeric(colnames(NET))
	SR = as.numeric(rownames(NET))
	step_SF = SF[2]-SF[1]
	step_SR = SR[2]-SR[1]
	step_L = 0.005
	Start = which.min(abs(SR+1))
	lim_CL = 3
	FF_all = matrix(0,24,2)
	rownames(FF_all) = c("Coverage","p_BAF","p_LRR","1","2","3","4","5", "copy_2","clust34","homo_clust","5_3","4_0","4_4","4_3","4_2","3_0","3_3","3_2","2_0","2_2","2_1","1_0","1_1")
	for (jj in 1:length(pp_N)) {
		NET_C = matrix(0,6,NF)
		rownames(NET_C) = SR_perf
		colnames(NET_C) = seq(1,NF,1)
		tt = 0
		for (cn in 1:5) {
			nbn = 1
			if (cn%%2==0) {
				NET_C[cn,1] = 1
			}
			for (nbc in trunc(cn/2+1):cn) {
				baf = ((1-pp_N[jj])*nbc+pp_N[jj]*nbn)/((1-pp_N[jj])*cn+2*pp_N[jj])
				if ((cn==2)&&(nbc==cn)&&(baf<=0.98)) {
					tt = 1
				}
				if ((cn>1)&&(nbc==cn)&&(tt!=0)) {
					NET_C[cn,min(95,which.min(abs(SF-baf)))] = 1
				} else {
					NET_C[cn,which.min(abs(SF-baf))] = 1
				}
			}
			nbn = 2
			nbc = cn
			baf = ((1-pp_N[jj])*nbc+pp_N[jj]*nbn)/((1-pp_N[jj])*cn+2*pp_N[jj])
			NET_C[cn,(dim(NET_C)[[2]]-3)] = 1

		}
		NET_C[6,] = NET_C[1,]+NET_C[2,]+NET_C[3,]+NET_C[4,]+NET_C[5,]
		NET_C = NET_C[,-which(NET_C[6,]==0)]
		NET_C[6,] = as.numeric(colnames(NET_C))
		NET_tmp = matrix(1,16,7)
		colnames(NET_tmp) = c("N_LRR","N_BAF","Step_LRR","Step_BAF_le","Step_BAF_ri","Copy","majorA")
		NET_tmp[,6] = c(5,5,5,5,4,4,4,4,3,3,3,2,2,2,1,1)
		NET_tmp[,7] = c(0,5,4,3,0,4,3,2,0,3,2,0,2,1,0,1)
		kk = 1
		for (kk in 1:length(qq_L)) {
			SR_C = SR_perf*qq_L[kk]
			for (k in 1:length(SR_C)) {
				SR_C[k] = round((SR_C[k]-SR[1])/step_L+1,0)
			}
			SR_C = SR_C-SR_C[2]+Start
			tt = as.numeric(NET_C[1:5,])
			tt = cbind(matrix(0,length(tt),6),tt)
			tt[,1] = rep(SR_C[1:5],dim(NET_C)[[2]])
			tt[,3] = rep(c(max(SR_C[2:5]-SR_C[1:4]),SR_C[2:5]-SR_C[1:4]),dim(NET_C)[[2]])
			tt[,2] = rep(NET_C[6,],each=5)
			tt = tt[-which(tt[,7]==0),]
			tt = tt[order(tt[,1],tt[,2],decreasing=TRUE),]
			NET_tmp[,1:5] = tt[,1:5]
			NET_tmp[,3] = round(NET_tmp[,3]/lim_CL,0)
			tt = which(NET_tmp[,3]>20)
			if (length(tt)>0) {
				NET_tmp[tt,3] = 20
			}
			tt = which(NET_tmp[,3]<6)
			if (length(tt)>0) {
				NET_tmp[tt,3] = round(NET_tmp[tt,3]*3/2-0.6)
			}
			NET_tmp[1:(dim(NET_tmp)[[1]]-1),4] = NET_tmp[1:(dim(NET_tmp)[[1]]-1),2]-NET_tmp[2:dim(NET_tmp)[[1]],2]-1
			tt = which(NET_tmp[,4]<0)
			if (length(tt)>0) {
				NET_tmp[tt,4] = NET_tmp[tt,2]-1
			}
			NET_tmp[dim(NET_tmp)[[1]],4] = max(0,NET_tmp[dim(NET_tmp)[[1]],2]-1)
			NET_tmp[2:dim(NET_tmp)[[1]],5] = NET_tmp[1:(dim(NET_tmp)[[1]]-1),2]-NET_tmp[2:dim(NET_tmp)[[1]],2]-1
			tt = which(NET_tmp[,5]<0)
			if (length(tt)>0) {
				NET_tmp[tt,5] = 2*(NF-NET_tmp[tt,2])
			}
			NET_tmp[1,5] = 2*max(1,NF-NET_tmp[1,2])
			NET_tmp[,4] = NET_tmp[,4]/2
			NET_tmp[,5] = NET_tmp[,5]/2
			NET_tmp[1,4] = min(NET_tmp[1,4],5)
			if (NET_tmp[2,2]>95) {
				NET_tmp[2,4] = min(NET_tmp[2,4],10)
			}
			NET_tmp[2,5] = min(max(min(NET_tmp[2,5],5),NET_tmp[2,5]/2),10)*2
			for (k in 2:dim(NET_tmp)[1]) {
				if (NET_tmp[k,2]==102) {
					NET_tmp[k,4] = min(NET_tmp[k,4],5)
					NET_tmp[k+1,5] = min(max(min(NET_tmp[k+1,5],5),NET_tmp[k+1,5]/2),10)*2
					NET_tmp[k-1,4] = min(NET_tmp[k-1,4]*2,10)
					if (NET_tmp[k+1,2]>95) {
						NET_tmp[k+1,4] = min(NET_tmp[k+1,4],10)
					}
				} else {
					if (NET_tmp[k,2]==1) {
						NET_tmp[k,5] = min(NET_tmp[k,5]*2/3,15)
					} else {
						NET_tmp[k,5] = min(NET_tmp[k,5]/2,10)
					}
					NET_tmp[k,4] = min(NET_tmp[k,4]/2,10)
				}
			}
			NET_tmp = round(NET_tmp+0.1,0)
			NET_tmp[,4] = NET_tmp[,2]-NET_tmp[,4]
			NET_tmp[,5] = NET_tmp[,2]+NET_tmp[,5]
			tt = which(NET_tmp[,4]<1)
			if (length(tt)>0) {
				NET_tmp[tt,4] = 1
			}
			tt = which(NET_tmp[,5]>NF)
			if (length(tt)>0) {
				NET_tmp[tt,5] = NF
			}
			pos = length(which(NET_tmp[,1]>0))
			outr = length(which(NET_tmp[,1]>NR))
			FF_max2 = rep(0,24)
			FF_max4 = rep(0,24)
			while (SR_C[2]<NR) {
				FF = rep(0,24)
				for (j in (outr+1):pos) {
					FF[8+j] = sum(NET[max(1,(NET_tmp[j,1]-NET_tmp[j,3])):min(NR,(NET_tmp[j,1]+NET_tmp[j,3])),NET_tmp[j,4]:NET_tmp[j,5]])
				}
				FF[1] = sum(FF)
				if (FF[1]>FF_max4[1]) {
					FF_max4 = c(FF[1],pp_N[jj],qq_L[kk],SR_C[1:5],FF[9:24])
				}
				FF[1] = sum(FF[17:24])
				if ((-FF[1])<FF_max2[1]) {
					FF_max2 = c(-FF[1],pp_N[jj],qq_L[kk],SR_C[1:5],FF[9:24])
				}
				NET_tmp[,1] = NET_tmp[,1]+1
				SR_C = SR_C+1
				if ((pos<dim(NET_tmp)[[1]])&&(NET_tmp[pos+1,1]>0)) {
					pos = length(which(NET_tmp[,1]>0))
				}
				if ((outr<dim(NET_tmp)[[1]])&&(NET_tmp[outr+1,1]>NR)) {
					outr = length(which(NET_tmp[,1]>NR))
				}
			}
			if (FF_max2[1]<0) {
				FF_all = cbind(FF_all,FF_max2)
			}
			if (FF_max4[1]>0) {
				FF_all = cbind(FF_all,FF_max4)
			}
		}
	}
	for (j in 4:8) {
		if (length(which(FF_all[j,]<0))>0) {
			FF_all[j,which(FF_all[j,]<0)] = 0
		}
		if (length(which(FF_all[j,]>NR))>0) {
			FF_all[j,which(FF_all[j,]>NR)] = 0
		}
	}
	FF_all[9,] = FF_all[20,]+FF_all[21,]+FF_all[22,]
	FF_all[10,] = FF_all[15,]+FF_all[19,]
	FF_all[11,] = FF_all[14,]+FF_all[18,]+FF_all[21,]+FF_all[24,]
	FF_all[12,] = sign(FF_all[14,]*FF_all[15,]*FF_all[18,]*FF_all[19,]*FF_all[2,])*500
	tt = which(FF_all[9,]+FF_all[12,]<500)
	if (length(tt)>0) {
		FF_all = FF_all[,-tt]
	}
	tt = which(FF_all[10,]>500)
	if (length(tt)>0) {
		FF_all[10,tt] = 1000
	}
	tt = which(FF_all[11,]-FF_all[10,]<0)
	ss = which(FF_all[19,]>2000)
	tt = setdiff(tt,ss)
	if (length(tt)>0) {
		FF_all = FF_all[,-tt]
	}
	FF_all = FF_all[1:11,]
	tt = which(FF_all[1,]<0)
	if (length(tt)>0) {
		FF_all2 = FF_all[,tt]
		if (dim(FF_all)[[2]]>length(tt)) {
			FF_all = FF_all[,-tt]
		} else {
			FF_all = FF_all2
		}
	} else {
		FF_all2 = NA
	}
	FF_all = abs(FF_all)
	FF_all = FF_all[,order(FF_all[1,],decreasing=TRUE)]
	FF_all[9,] = 0
	FF_all[10,] = round(FF_all[1,]/FF_all[1,1]*100,1)
	tt = which(FF_all[10,]<70)
	if (length(tt)>0) {
		FF_all = FF_all[,-tt]
	}
	for (k in 1:(dim(FF_all)[[2]]-1)) {
		if (FF_all[9,k]==0) {
		 	FF_all[9,k] = k
		 	kn = 1
			for (kk in (k+1):dim(FF_all)[[2]]) {
				if ((max(abs(FF_all[5:7,k]-FF_all[5:7,kk]))<10)&&(abs(FF_all[2,k]-FF_all[2,kk])<0.1)) {
					if ((FF_all[10,k]-FF_all[10,kk]<10)&&(kn<10)) {
						FF_all[9,kk] = k
    	                kn = kn+1
    	            } else {
    	            	FF_all[9,kk] = (-1)
    	            }
				}
			}
		}
	}
	k = dim(FF_all)[2]
	if (FF_all[9,k]==0) {
		FF_all[9,k] = k
	}
	tt = which(FF_all[9,]==(-1))
	if (length(tt)>1) {
		if (length(tt)>0) {
			FF_all = FF_all[,-tt]
		}
	} else {
		FF_all = cbind(FF_all,FF_all)
		tt = which(FF_all[9,]==(-1))
		if (length(tt)>0) {
			FF_all = FF_all[,-tt]
		}
	}
	FF_all = FF_all[,order(FF_all[9,],decreasing=FALSE)]
	FF_all = FF_all[1:10,]
	rownames(FF_all) = c("Coverage_SNP","p_BAF","q_LRR","1_copy","2_copy","3_copy","4_copy","5_copy","Set_pattern","Percent_coverage")
	return(invisible(FF_all))
}
