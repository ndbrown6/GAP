'CopyNumber_and_Genotype' <- function (seg_LRR, seg_BAF, seg_Len, showGAPTemplate=FALSE)
{
	tmp = matrix(0,length(seg_LRR),6)
	colnames(tmp) = c("Index", "LRR", "BAF", "Len", "sup_L", "sup_B")
	tmp[,1] = seq(1:length(seg_LRR))
	tmp[,2:4] = cbind(seg_LRR,seg_BAF,seg_Len)
	C_ind = 1
	C_lrr = 2
	C_baf = 3
	C_len = 4
	C_nL = 5
	C_nB = 6
	CN = 8
	R_min = round((min(tmp[which(tmp[,C_len]>10),C_lrr],na.rm=TRUE)-0.15)*10,0)/10
	R_max = round((max(tmp[which(tmp[,C_len]>10),C_lrr],na.rm=TRUE)+0.15)*10,0)/10
	step_SR = 0.005
	step_SF = 0.005
	NR = round(((R_max-R_min)/step_SR),0)
	SR = seq((R_min+step_SR),R_max,step_SR)
	NF = round(((1.02+step_SF-0.5)/step_SF),0)
	SF = seq(0.5+step_SF,1.02+step_SF,step_SF)
	NET = matrix(0,NR,NF)
	rownames(NET) = SR
	colnames(NET) = SF
	for (k in 1:dim(tmp)[1]) {
		tmp[k,C_nL] = which.min(abs(SR-tmp[k,C_lrr]))
		if (tmp[k,C_lrr]>=SR[tmp[k,C_nL]]) {
			tmp[k,C_nL] = min(tmp[k,C_nL]+1,NR)
		}
		tmp[k,C_nB] = which.min(abs(SF-tmp[k,C_baf]))
		if (tmp[k,C_baf]>=SF[tmp[k,C_nB]]) {
			tmp[k,C_nB] = min(tmp[k,C_nB]+1,NF)
		}
		if (tmp[k,C_len]>=50) {
			NET[tmp[k,C_nL],tmp[k,C_nB]] = NET[tmp[k,C_nL],tmp[k,C_nB]]+tmp[k,C_len]
		}
	}
	NET = smoothNET(NET,weightC=1)
	NET[which(NET<100)] = 0
	FF_all = fittingGAP(NET)
	kk = as.numeric(rownames(table(FF_all[9,])))[1]
	tt = which(FF_all[9,]!=kk)
	if (length(tt)>0) {
		if (length(tt)==(dim(FF_all)[[2]]-1)) {
			FF_final = cbind(FF_all,FF_all)
		} else {
			FF_final = FF_all
		}
		FF_final = FF_final[,-which(FF_final[9,]!=kk)]
	}
	p_BAF = mean(FF_final[2,])
	q_LRR = mean(FF_final[3,])
	for (k in 4:8) {
		if (sum(FF_final[k,])>0) {
			FF_final[k,1] = round(mean(FF_final[k,which(FF_final[k,]!=0)]),0)
		}
	}
	FF_final[4:8,1] = round(FF_final[4:8,1],0)
	Delta = which(FF_final[5:7,1]!=0)
	Delta = c(Delta[1]+1, SR[FF_final[Delta[1]+4,1]])
	Centr_CL = getCenterClusterGAP(CN,p_BAF,q_LRR,Delta)
	NETM = GetGAPtemplateZones(NET,Centr_CL)
	SB_CL = subClonalClusters(NETM,Centr_CL)
	for (k in 1:dim(Centr_CL)[2]) {
		tt = which(NETM[,3]==Centr_CL[3,k])
		if (length(tt)==0) {
			NETM = rbind(NETM,c(Centr_CL[5,k],Centr_CL[4,k],Centr_CL[3,k],Centr_CL[3,k]))
		}
	}
	NETM = NETM[order(NETM[,3],decreasing=FALSE),]
	if (showGAPTemplate) {
		tt = plotGAPdensity(NET)
		for (j in 1:dim(Centr_CL)[2]) {
			text(Centr_CL[4,j],Centr_CL[5,j],paste(as.character(Centr_CL[1,j]),"/",as.character(Centr_CL[2,j]),sep=""),cex=.9, col="grey40")
		}
		text(0.875,-0.7,paste("p_BAF = ",as.character(round(p_BAF,2)),sep=""),cex=.9, pos=4)
		text(0.875,-0.8,paste("q_LRR = ",as.character(round(q_LRR,2)),sep=""),cex=.9, pos=4)
		text(0.875,-0.9,paste("2 Copy LRR = ",as.character(round(Centr_CL[5,4],2)),sep=""),cex=.9, pos=4)
	}
	tmp = list(NETM=NETM, NET=NET, Centr_CL=Centr_CL, SB_CL=SB_CL, coeff=c(p_BAF,q_LRR))
	return(invisible(tmp))
}
