#==================================================
# David Brown
# dbrown@ulb.ac.be
# http://homepages.ulb.ac.be/~dbrown/
#==================================================
source(paste0(.libPaths(), "/GAP/etc/config.R"))
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
	
		## summary stats
		info = summaryStats(tmp3, SNP_annot=SNP_annot)
		GAP_Params[sampleNames[i], 4:7] = as.numeric(info)
		write.csv(GAP_Params, file="../GAP_Params/GAP_Params.csv")
	}
}
close(pb)
