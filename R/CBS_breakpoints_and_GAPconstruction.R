'CBS_breakpoints_and_GAPconstruction' <- function(datL, datA, SNP_annot, germHomoThr=0.93)
{
	nSNP = dim(SNP_annot)[1]
	dataF = matrix(0,nSNP,4)
	colnames(dataF) = c("BP_CN", "BP_AD", "AD_hete", "Homo_segm")
	c_bpLRR = 1
	c_bpBAF = 2
	c_hetAD = 3
	c_homo = 4
	C_Ind = 1
	C_Chr = 2
	C_Pos = 4

	## segmentation of Log2 Ratio
	tt = which(abs(datL)>100)
	if (length(tt)>0) {
		datL[tt] = NA
	}
	notNAs = which(!is.na(datL))
	CNA.object = CNA(genomdat=datL[notNAs], chrom=round((SNP_annot[notNAs,C_Chr]+(SNP_annot[notNAs,C_Chr+1]-1)/2),1),
				 	 maploc=SNP_annot[notNAs,C_Ind], data.type="logratio", sampleid=NULL)
	segment.CNA.object = segment(CNA.object, alpha=0.1, nperm=100, p.method="perm", min.width=5, eta=0.05, trim=0.025,
								 undo.splits="sdundo", undo.SD=3, verbose=0)
									 
	tmp = segment.CNA.object$output
	dataF[tmp$loc.end,c_bpLRR] = 1

	## segmentation of B Allele frequency
	tt = which(abs(datA)>1)
	if (length(tt)>0) {
		datA[tt] = 1
	}
	tt = which(is.na(datA))
	if (length(tt)>0) {
		datA[tt] = 1
	}
	datA = abs(datA-0.5)+0.5

	homoA = datA*0
	tt = which(abs(datA)>germHomoThr)
	if (length(tt)>0) {
		homoA[tt] = 1
	}
	S_homo = homoA
	rh = S_homo*0
	lh = S_homo*0
	for (k in 1:9) {
		lh[10:length(lh)] = lh[10:length(lh)]+S_homo[(9-k+1):(length(S_homo)-k)]
		rh[1:(length(rh)-10+1)] = rh[1:(length(rh)-10+1)]+S_homo[(k+1):(length(S_homo)-9+k)]
	}
	S_hete = S_homo*0
	tt = which(rh+lh<15)
	if (length(tt)>0) {
		S_hete[tt] = 1
	}
	tt = which(rh<9)
	if (length(tt)>0) {
		S_hete[tt] = 1
	}
	tt = which(lh<9)
	if (length(tt)>0) {
		S_hete[tt] = 1
	}
	S_homo = (1-S_hete)
	rh = S_homo*0
	lh = S_homo*0

	for (k in 1:49) {
		lh[50:length(lh)]<-lh[50:length(lh)]+S_homo[(49-k+1):(length(S_homo)-k)]
		rh[1:(length(rh)-50+1)]<-rh[1:(length(rh)-50+1)]+S_homo[(k+1):(length(S_homo)-49+k)]
	}
	S_homo = S_homo*0
	tt = which(rh+lh>90)
	if (length(tt)>0) {
		S_homo[tt] = 1
	}
	tt = which(rh>45)
	if (length(tt)>0) {
		S_homo[tt] = 1
	}
	tt = which(lh>45)
	if (length(tt)>0) {
		S_homo[tt] = 1
	}
	S_homo = S_homo*(1-S_hete)
	seg = abs(S_homo[2:nSNP]-S_homo[1:(nSNP-1)])
	seg = c(1,which(seg!=0),nSNP)
	seg = as.numeric(rownames(table(seg)))
	for (k in 1:(length(seg)-1)) {
		if (S_homo[seg[k]+1]==0) {
			if (sum(1-homoA[(seg[k]+1):seg[k+1]])<5) {
				S_homo[(seg[k]+1):seg[k+1]] = 1
			}
		}
	}
	seg = abs(S_homo[2:nSNP]-S_homo[1:(nSNP-1)])+abs(SNP_annot[2:nSNP,C_Chr]-SNP_annot[1:(nSNP-1),C_Chr])+abs(SNP_annot[2:nSNP,C_Chr+1]-SNP_annot[1:(nSNP-1),C_Chr+1])
	seg = c(1,which(seg!=0),nSNP)
	seg = as.numeric(rownames(table(seg)))
	for (k in 2:length(seg)) {
		if (seg[k]-seg[k-1]+1<100) {
			S_homo[seg[k-1]:seg[k]] = 0
		}
	}
	dataF[,c_homo] = S_homo
	dataF[,c_hetAD] = 1-homoA
	tmp = cbind(datA,round((SNP_annot[,C_Chr]+(SNP_annot[,C_Chr+1]-1)/2),1),SNP_annot[,C_Ind])
	tmp = data.matrix(tmp)
	tt = which((1-dataF[,c_hetAD])+dataF[,c_homo]==0)
	if (length(tt)>0) {
		tmp = tmp[tt,]
	}
	CNA.object = CNA(genomdat=tmp[,1], chrom=tmp[,2], maploc=tmp[,3], data.type="logratio", sampleid=NULL)
	segment.CNA.object = segment(x=CNA.object, alpha=0.1, nperm=100, p.method="perm", min.width=5, eta=0.05,
							 	 trim=0.025, undo.splits="sdundo", undo.SD=2, verbose=0)
	tmp = segment.CNA.object$output
	dataF[tmp$loc.end,c_bpBAF] = 1
	seg = abs(SNP_annot[1:(nSNP-1),C_Chr]-SNP_annot[2:(nSNP-0),C_Chr])+abs(SNP_annot[1:(nSNP-1),C_Chr+1]-SNP_annot[2:(nSNP-0),C_Chr+1])+dataF[1:(nSNP-1),c_bpLRR]+dataF[1:(nSNP-1),c_bpBAF]
	seg = which(seg!=0)
	seg = c(seg,nSNP)
	seg = cbind(data.matrix(SNP_annot[seg,]),dataF[seg,])
	tmp = matrix(0,dim(seg)[1],10)
	colnames(tmp) = GAP::.GapEnv$ColNamesTMP[1:10]
	tmp[,GAP::.GapEnv$c_ind] = seq(1,dim(tmp)[1])
	tmp[,GAP::.GapEnv$c_lrrBP] = seg[,c_bpLRR+4]
	tmp[,GAP::.GapEnv$c_bafBP] = seg[,c_bpBAF+4]
	tmp[,GAP::.GapEnv$c_chr] = seg[,C_Chr]+(seg[,C_Chr+1]-1)/2
	tmp[,GAP::.GapEnv$c_if] = seg[,C_Ind]
	tmp[,GAP::.GapEnv$c_is] = c(1,seg[1:(dim(seg)[1]-1),C_Ind]+1)
	tmp[,GAP::.GapEnv$c_len] = tmp[,GAP::.GapEnv$c_if]-tmp[,GAP::.GapEnv$c_is]+1
	Flag = TRUE
	while (Flag) {
		Flag = FALSE
		for (k in 2:(dim(tmp)[1]-2)) {
			if (tmp[k,GAP::.GapEnv$c_chr]!=tmp[k+1,GAP::.GapEnv$c_chr]) {
				if ((tmp[k,GAP::.GapEnv$c_len]<5)&&(tmp[k,GAP::.GapEnv$c_chr]==tmp[k-1,GAP::.GapEnv$c_chr])) {
					tmp[k-1,GAP::.GapEnv$c_if] = tmp[k,GAP::.GapEnv$c_if]
					tmp[k-1,GAP::.GapEnv$c_len] = tmp[k-1,GAP::.GapEnv$c_len]+tmp[k,GAP::.GapEnv$c_len]
					tmp[k,GAP::.GapEnv$c_len] = 0
				}
				if ((tmp[k+1,GAP::.GapEnv$c_len]<5)&&(tmp[k+1,GAP::.GapEnv$c_chr]==tmp[k+2,GAP::.GapEnv$c_chr])) {
					tmp[k+2,GAP::.GapEnv$c_is] = tmp[k+1,GAP::.GapEnv$c_is]
					tmp[k+2,GAP::.GapEnv$c_len] = tmp[k+1,GAP::.GapEnv$c_len]+tmp[k+2,GAP::.GapEnv$c_len]
					tmp[k+1,GAP::.GapEnv$c_len] = 0
				}
			}
		}
		if ((tmp[1,GAP::.GapEnv$c_len]<5)&&(tmp[1,GAP::.GapEnv$c_chr]==tmp[2,GAP::.GapEnv$c_chr])) {
			tmp[2,GAP::.GapEnv$c_is] = tmp[1,GAP::.GapEnv$c_is]
			tmp[2,GAP::.GapEnv$c_len] = tmp[1,GAP::.GapEnv$c_len]+tmp[2,GAP::.GapEnv$c_len]
			tmp[1,GAP::.GapEnv$c_len] = 0
		}
		tt = dim(tmp)[1]
		if ((tmp[tt,GAP::.GapEnv$c_len]<5)&&(tmp[tt,GAP::.GapEnv$c_chr]==tmp[tt-1,GAP::.GapEnv$c_chr])) {
			tmp[tt-1,GAP::.GapEnv$c_if] = tmp[tt,GAP::.GapEnv$c_if]
			tmp[tt-1,GAP::.GapEnv$c_len] = tmp[tt-1,GAP::.GapEnv$c_len]+tmp[tt,GAP::.GapEnv$c_len]
			tmp[tt,GAP::.GapEnv$c_len] = 0
		}
		tt = which(tmp[,GAP::.GapEnv$c_len]==0)
		if (length(tt)!=0) {
			tmp = tmp[-tt,]
			Flag = TRUE
		}
	}
	for (k in 1:dim(tmp)[1]) {
		if (tmp[k,GAP::.GapEnv$c_if]>tmp[k,GAP::.GapEnv$c_is]) {
			tmp[k,GAP::.GapEnv$c_lrr] = median(datL[tmp[k,GAP::.GapEnv$c_is]:tmp[k,GAP::.GapEnv$c_if]], na.rm=TRUE)
		} else {
			tmp[k,GAP::.GapEnv$c_lrr] = datL[tmp[k,GAP::.GapEnv$c_is]]
		}
	}	
	tmp = getBAFforRecognition(tmp, datA, germHomoThr)
	tmp = round(tmp, 4)
	return(invisible(tmp))
}
