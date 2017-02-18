#==================================================
# David Brown
# dbrown@ulb.ac.be
# http://homepages.ulb.ac.be/~dbrown/
#==================================================
rm(list=ls(all=TRUE))
source('config.R')
wD = getwd()
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
	
		## summary per chromosome
		if (!file.exists(paste("../Results/", sampleNames[i], "/Summary_per_Chrom/", sep=""))) {
			dir.create(path=paste("../Results/", sampleNames[i], "/Summary_per_Chrom/", sep=""))
		}
		setwd(paste("../Results/", sampleNames[i], "/Summary_per_Chrom/", sep=""))
		summaryPlots_perChrom(tmp3, sampleName=sampleNames[i], SNP_annot, BAf_R, LogR_R, resWOfit$Centr_CL)
		setwd(wD)
	
		## updated GAP pattern recognition
		pdf(file=paste("../Results/", sampleNames[i], "/GAP_Pattern_Final.pdf", sep=""))
		GAPplot(tmp3, sampleName=sampleNames[i], SNP_annot, resWOfit$Centr_CL, resWOfit$NETM, resWOfit$coeff[1], resWOfit$coeff[2])
		dev.off()
	
		## summary genome plot and GAP
		jpeg(filename=paste("../Results/", sampleNames[i], "/GAP_Summary.jpeg", sep=""), width=1920, height=2400, pointsize=24)
		info = summaryPlots(tmp3, sampleName=sampleNames[i], SNP_annot=SNP_annot, datA=BAf_R, datL=LogR_R, resWOfit$Centr_CL, resWOfit$NETM, p_BAF=as.numeric(GAP_Params[sampleNames[i],"p_BAF"]), q_LRR=as.numeric(GAP_Params[sampleNames[i],"q_LRR"]))
		dev.off()
	}
}
close(pb)
