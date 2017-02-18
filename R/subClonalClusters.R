'subClonalClusters' <- function(NETM, Centr_CL)
{
	CN_CL = rbind(Centr_CL[1,],Centr_CL[5,])
	CN_CL = table(CN_CL[1,],CN_CL[2,])
	CN_CL = rbind(as.numeric(rownames(CN_CL)),as.numeric(colnames(CN_CL)))
	BAF_CL = rbind(Centr_CL[1,],Centr_CL[2,],Centr_CL[4,])
	BAF_CL = BAF_CL[,order(BAF_CL[3,],decreasing=FALSE)]
	SB_CL = numeric(0)
	if (length(which(NETM[,3]>=40))>1) {
		if (length(which(NETM[,3]<40))>0) {
			SB_CL = NETM[-which(NETM[,3]<40),]
			SB_CL = SB_CL[order(SB_CL[,3],decreasing=FALSE),]
			SB_CL[,1] = SB_CL[,1]*SB_CL[,4]
			SB_CL[,2] = SB_CL[,2]*SB_CL[,4]
			tt = SB_CL[1,]
			k = 4
			for (k in 2:dim(SB_CL)[1]) {
				if (SB_CL[k,3]==tt[3]) {
					SB_CL[k,1] = SB_CL[k,1]+SB_CL[k-1,1]
					SB_CL[k,2] = SB_CL[k,2]+SB_CL[k-1,2]
					SB_CL[k,4] = SB_CL[k,4]+SB_CL[k-1,4]
					SB_CL[k-1,2] = 0
				} else {
					SB_CL[k-1,1] = SB_CL[k-1,1]/SB_CL[k-1,4]
					SB_CL[k-1,2] = SB_CL[k-1,2]/SB_CL[k-1,4]
					tt = SB_CL[k,]
				}
			}
			SB_CL[k,1] = SB_CL[k,1]/SB_CL[k,4]
			SB_CL[k,2] = SB_CL[k,2]/SB_CL[k,4]
			if (length(which(SB_CL[,2]==0))>0) {
				if (length(which(SB_CL[,2]!=0))==1) {
					SB_CL = rbind(SB_CL,cbind(SB_CL[,1:2],SB_CL[,3:4]+1))
				}
				SB_CL = SB_CL[-which(SB_CL[,2]==0),]
				SB_CL = SB_CL[order(SB_CL[,3],decreasing=FALSE),]
			}
			SB_CL = cbind(SB_CL[,3],SB_CL[,2:4],SB_CL[,1:2])
			SB_CL[,2:4] = 0
			SB_CL = t(SB_CL)
			rownames(SB_CL) = c("Subclone","CN_LRR","B_All","CN_BAF","SB_LRR","SB_BAF")
			for (k in 1:dim(SB_CL)[2]) {
				p = 1
				while ((p<=dim(CN_CL)[2])&&(SB_CL[5,k]>CN_CL[2,p])) {
					p = p+1
				}
				if (p>dim(CN_CL)[2]) {
					SB_CL[2,k] = CN_CL[1,p-1]
				} else {
					if (p==1) {
						SB_CL[2,k] = CN_CL[1,p]
					} else {
						w = c((CN_CL[2,p]-SB_CL[5,k]),(SB_CL[5,k]-CN_CL[2,p-1]))
						w = w/(CN_CL[2,p]-CN_CL[2,p-1])
						SB_CL[2,k] = weighted.mean(c(CN_CL[1,p-1],CN_CL[1,p]),w)
					}
				}
				SB_CL[2,k] = round(SB_CL[2,k],1)
				if ((p>1)&&(SB_CL[2,k]==CN_CL[1,p-1])) {
					SB_CL[2,k] = SB_CL[2,k]+0.1
				}
				if ((p<=dim(CN_CL)[[2]])&&(SB_CL[2,k]==CN_CL[1,p])) {
					SB_CL[2,k] = SB_CL[2,k]-0.1
				}
				if (SB_CL[6,k]>=0.97) {
					SB_CL[3,k] = SB_CL[2,k]
					SB_CL[4,k] = SB_CL[2,k]
				} else {
					if ((p>1)&&(length(which(BAF_CL[1,]<CN_CL[1,p-1]))>0)) {
						tmp_BAF_CL = BAF_CL[,-which(BAF_CL[1,]<CN_CL[1,p-1])]
					} else {
						tmp_BAF_CL = BAF_CL
					}
					if (p<=dim(CN_CL)[2]) {
						tt = which(tmp_BAF_CL[1,]>CN_CL[1,p])
						if (length(tt)>0) {
							if (length(tt)==(dim(tmp_BAF_CL)[[2]]-1)) {
								tmp_BAF_CL = cbind(tmp_BAF_CL,tmp_BAF_CL)
							}
							tmp_BAF_CL = tmp_BAF_CL[,-which(tmp_BAF_CL[1,]>CN_CL[1,p])]
						}
					}
					tmp_BAF_CL = cbind(tmp_BAF_CL,tmp_BAF_CL[,dim(tmp_BAF_CL)[2]])
					tmp_BAF_CL[3,dim(tmp_BAF_CL)[2]] = 1
					p = 1
					while ((p<dim(tmp_BAF_CL)[2])&&(SB_CL[6,k]>tmp_BAF_CL[3,p])) {
						p = p + 1
					}
					if (p==1) {
						p = p + 1
					}
					w = c((tmp_BAF_CL[3,p]-SB_CL[6,k]),(SB_CL[6,k]-tmp_BAF_CL[3,p-1]))
					w = w/(tmp_BAF_CL[3,p]-tmp_BAF_CL[3,p-1])
					if (tmp_BAF_CL[1,p-1]!=0) {
						SB_CL[3,k] = weighted.mean(c(tmp_BAF_CL[2,p-1]/tmp_BAF_CL[1,p-1],tmp_BAF_CL[2,p]/tmp_BAF_CL[1,p]),w)
					} else {
						SB_CL[3,k] = weighted.mean(c(tmp_BAF_CL[2,p-1],tmp_BAF_CL[2,p]),w)
					}
					SB_CL[3,k] = round(SB_CL[3,k]*SB_CL[2,k],1)
					if (p!=dim(tmp_BAF_CL)[2]) {
						SB_CL[4,k] = weighted.mean(c(tmp_BAF_CL[1,p-1],tmp_BAF_CL[1,p]),w)
						SB_CL[4,k] = round(SB_CL[4,k],1)
						if (SB_CL[4,k]==tmp_BAF_CL[1,p-1]) {
							SB_CL[4,k] = SB_CL[4,k]+0.1
						}
						if (SB_CL[4,k]==tmp_BAF_CL[1,p]) {
							SB_CL[4,k] = SB_CL[4,k]-0.1
						}
					}
				}
			}
		}
	}
	return(invisible(SB_CL))
}
