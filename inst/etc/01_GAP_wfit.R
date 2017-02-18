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
		if (!file.exists("../Results/")) {
			dir.create(path="../Results/")
		}
		if (!file.exists(paste("../Results/", sampleNames[i], "/", sep=""))) {
			dir.create(path=paste("../Results/", sampleNames[i], "/", sep=""))
		}
		load(paste("../Segmented/", sampleNames[i], ".RData", sep=""))
	
		## copy number and b allele frequency segmentation
		tmp = CBS_breakpoints_and_GAPconstruction(datL=LogR_S, datA=BAf_S, SNP_annot=SNP_annot, germHomoThr=0.99)
		notNAs = sort(intersect(which(!is.na(tmp[,GAP::.GapEnv$c_lrr]+tmp[,GAP::.GapEnv$c_baf])),which(tmp[,GAP::.GapEnv$c_chr]<23)))
	
		## fit GAP
		resWfit = CopyNumber_and_Genotype(tmp[notNAs,GAP::.GapEnv$c_lrr],tmp[notNAs,GAP::.GapEnv$c_baf],tmp[notNAs,GAP::.GapEnv$c_len],showGAPTemplate=FALSE)
		
		## write parameters to file
		GAP_Params[sampleNames[i],"p_BAF"] = resWfit$coeff[1]
		GAP_Params[sampleNames[i],"q_LRR"] = resWfit$coeff[2]
		GAP_Params[sampleNames[i],"X2copy_LRR"] = resWfit$Centr_CL[5,3]
		write.csv(GAP_Params, file="../GAP_Params/GAP_Params.csv")

		## plot preliminary GAP pattern recognition
		pdf(file=paste("../Results/", sampleNames[i], "/GAP_Preliminary_pattern_wfit.pdf", sep=""))
		GAPplot_preliminary(tmp, resWfit$Centr_CL, resWfit$coeff[1], resWfit$coeff[2], resWfit$NETM, sampleName=sampleNames[i])
		dev.off()

		## save data
		save(tmp, notNAs, resWfit, file=paste("../Results/", sampleNames[i], "/CopyNumber_and_Genotype.RData", sep=""))
	}
}
close(pb)
