'markNETtemplate' <- function(NET, NET_tmp)
{
	NF = dim(NET)[[2]]
	NR = dim(NET)[[1]]
	SF = as.numeric(colnames(NET))
	SR = as.numeric(rownames(NET))
	step_SF = SF[2]-SF[1]
	step_SR = SR[2]-SR[1]
	NET1 = NET
	NET1 = NET1*0
	NET2 = NET1*0
	NET3 = NET1*0
	p = 1
	for (p in 1:dim(NET_tmp)[1]) {
		m = max(1,(NET_tmp[p,1]-NET_tmp[p,3]))
		k = min((NET_tmp[p,1]+NET_tmp[p,4]),NR)
		s = (NET_tmp[p,2]-NET_tmp[p,5])
		n = (NET_tmp[p,2]+NET_tmp[p,6])
		NET1[m:k,s:n] = (-NET_tmp[p,9])
		NET3[m:k,s:n] = (-NET_tmp[p,9])
		j = 50	
		while (j>=1) {
			NET2[max(1,m-j):min(NR,k+j),max(1,s-j):min(NF,n+j)] = j
			j = j-1
		}
		for (j in max(1,m-50):min(NR,k+50)) {
			for (jj in max(1,s-50):min(NF,n+50)) {
				if ((NET1[j,jj]>0)&&(NET2[j,jj]<NET1[j,jj])||(NET1[j,jj]==0)) {
					NET1[j,jj] = NET2[j,jj]
					NET3[j,jj] = NET_tmp[p,9]
				}
			}
		}
	}
	for (j in 3:(dim(NET3)[1]-3)) {
		for (p in 3:(dim(NET3)[[2]]-3-8)) {
			tt = as.numeric(abs(NET3[j:(j+1),p:(p+1)]))
			if (length(which(tt==0))>0) {
				tt = tt[-which(tt==0)]
			}
			if ((length(tt)>0)&&(min(tt)!=max(tt))) {
				NET3[(j-2):(j+3),(p-2):(p+3)] = 0
			}
		}
	}
	for (j in 3:(dim(NET3)[1]-3)) {
		for (p in (dim(NET3)[2]-3-8):(dim(NET3)[2]-1)) {
			tt = as.numeric(abs(NET3[j:(j+1),p:(p+1)]))
			if (length(which(tt==0))>0) {
				tt = tt[-which(tt==0)]
			}
			if ((length(tt)>0)&&(min(tt)!=max(tt))) {
				NET3[(j-2):(j+3),p:(p+1)] = 0
			}
		}
	}
	for (j in 1:dim(NET3)[1]) {
		for (p in 1:dim(NET3)[2]) {
			if (NET3[j,p]>0) {
				tt = NET3[max(1,j-6):min(NR,j+6),max(1,p-6):min(NF,p+6)]+NET3[j,p]
				if (length(which(tt==0))==0) {
					NET3[j,p] = 0
				}
			}
		}
	}
	NET3[,(dim(NET3)[[2]]-2):dim(NET3)[[2]]] = 0
	for (p in 1:dim(NET1)[1]) {
		tt = which(NET1[p,]>0)
		if (length(tt)>0) {
			NET1[p,tt] = 0
		}
		tt = which(NET3[p,]>=0)
		if (length(tt)>0) {
			NET1[p,tt] = 0
		}
	}
	NET = NET*(1+sign(NET1))+NET1
	NETM = as.numeric(NET)
	NETM = cbind(NETM,NETM,NETM)
	NETM[,1] = rep(seq(1,length(SR)),length(SF))
	NETM[,2] = rep(seq(1,length(SF)),each=length(SR))
	NETM = NETM[-which(NETM[,3]==0),]
	NETM = cbind(NETM,NETM[,3])
	NETM = NETM[order(NETM[,3],decreasing=FALSE),]
	p = 1
	while ((p<=dim(NETM)[[1]])&&(NETM[p,3]<0)) {
		p = p+1
	}
	if (p<dim(NETM)[[1]]) {
		NETM[1:(p-1),4] = 4
	}
	for (k in 1:(p-1)) {
		if (p<=dim(NETM)[1]) {
			NETM[p:dim(NETM)[1],4] = abs(NETM[p:dim(NETM)[1],1]-NETM[k,1])
			tt = which(NETM[p:dim(NETM)[1],4]<4)
			if (length(tt)>0) {
				NETM[(p-1+tt),4] = abs(NETM[(p-1+tt),2]-NETM[k,2])
				tt = which(NETM[p:dim(NETM)[[1]],4]<4)
				if (length(tt)>0) {
					NETM = NETM[-(p-1+tt),]
				}
			}
		}
	}
	NETM[,4] = NETM[,3]
	NETM[,3] = NETM[,4]
	NETM = NETM[order(NETM[,3],decreasing=TRUE),]
	p = 1
	while (NETM[p,3]>0) {
		p = p + 1
	}
	if (p>1) {
		Dist = as.matrix(dist(NETM[1:(p-1),1:2],method="maximum",diag=TRUE,upper=TRUE))
		Dist[which(Dist>1)] = 0
		for (k in 1:dim(Dist)[1]) {
			Dist[k,k] = 1
		}
		cn_CL = -40
		for (k in 1:dim(Dist)[1]) {
			if (NETM[k,3]>0) {
				NETM[which(Dist[,k]!=0),3] = cn_CL
				if ((k+1)<=dim(Dist)[1]) {
					for (m in (k+1):dim(Dist)[1]) {
						if (sum(Dist[,k]*Dist[,m])>0) {
							Dist[,k] = Dist[,k]+Dist[,m]
							NETM[which(Dist[,k]!=0),3] = cn_CL
						}
					}
				}
				cn_CL = cn_CL-1
			}
		}
	}
	NET = NET*0
	for (k in 1:dim(NETM)[1]) {
		NET[NETM[k,1],NETM[k,2]] = NETM[k,3]
	}
	NET = abs(NET)
	NETM = abs(NETM)
	NETM[,1] = SR[NETM[,1]]
	NETM[,2] = SF[NETM[,2]]
	NETM = NETM[order(NETM[,3],decreasing=FALSE),]
	return(invisible(NETM))
}
