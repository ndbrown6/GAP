#==================================================
# David Brown
# dbrown@ulb.ac.be
# http://homepages.ulb.ac.be/~dbrown/
#==================================================
rm(list=ls(all=TRUE))
source('config.R')
GAP_Params = read.csv(file="../GAP_Params/GAP_Params.csv", header=TRUE, sep=",", row.names=1, stringsAsFactors=FALSE)
sampleNames = rownames(GAP_Params)
pb = txtProgressBar(min=1, max=length(sampleNames), style=3)
for (i in 1:length(sampleNames)) {
	setTxtProgressBar(pb, i)
	if (GAP_Params[sampleNames[i],"Array_QC"]==1) {
		load(paste("../Raw/", sampleNames[i], ".RData", sep=""))
		load(paste("../Segmented/", sampleNames[i], ".RData", sep=""))
		load(paste("../Results/", sampleNames[i], "/CopyNumber_and_Genotype.RData", sep=""))
		tmp2 = cbind(tmp, matrix(0,dim(tmp)[1],16))
		colnames(tmp2)[11:26] = GAP::.GapEnv$ColNamesTMP[11:26]
		tmp2[,GAP::.GapEnv$c_nC] = getCNwithPrecision(tmp2[,GAP::.GapEnv$c_lrr], resWOfit$Centr_CL,Precision=1/8)
		tmp2 = getGTwithRank(tmp2,resWOfit$Centr_CL,resWOfit$NETM)
		tmp2[,GAP::.GapEnv$c_CN] = tmp2[,GAP::.GapEnv$c_CN1]
		tmp2[,GAP::.GapEnv$c_BA] = tmp2[,GAP::.GapEnv$c_BA1]
		tmp3 = getConfidence(tmp2)
		tmp3 = shrinkReprTMP(tmp3)
		tmp3 = updateMeans_FinalPolish(tmp3, BAf_S, LogR_S, germHomozyg.mBAF.thr=0.99)
		PosSNPstart = SNP_annot[tmp3[,2],4]
		PosSNPend = SNP_annot[tmp3[,3],4]
		tmp4 = cbind(tmp3[,1:3], PosSNPstart, PosSNPend, tmp3[,4:ncol(tmp3)])
		save(tmp4, file=paste("../Results/", sampleNames[i], "/GAP_Segmented_Results.RData", sep=""))
		tmp4[,"Chr"] = floor(tmp4[,"Chr"])
		tmp4 = tmp4[,c("Chr", "PosSNPstart", "PosSNPend", "LRR", "MajorBAF", "Copy", "Mallele"),drop=FALSE]
		write.table(x=tmp4, file=paste("../Results/", sampleNames[i], "/GAP_Segmented_Results.csv", sep=""), append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE, sep=";")
	}
}
close(pb)
