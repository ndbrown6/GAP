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
		load(paste("../Results/", sampleNames[i], "/CopyNumber_and_Genotype.RData", sep=""))

		# p_BAF, q_LRR, Delta
		p_BAF = as.numeric(GAP_Params[sampleNames[i],"p_BAF"])
		q_LRR = as.numeric(GAP_Params[sampleNames[i],"q_LRR"])
		Delta = c(2, as.numeric(GAP_Params[sampleNames[i],"X2copy_LRR"]))
	
		## fit GAP manually
		resWOfit = CopyNumber_and_GenotypeWOfitting(tmp[notNAs,GAP::.GapEnv$c_lrr], tmp[notNAs,GAP::.GapEnv$c_baf], tmp[notNAs,GAP::.GapEnv$c_len], showGAPTemplate=FALSE, p_BAF, q_LRR, Delta)

		## plot preliminary GAP pattern recognition
		pdf(file=paste("../Results/", sampleNames[i], "/GAP_Preliminary_pattern_wofit.pdf", sep=""))
		GAPplot_preliminary(tmp, resWOfit$Centr_CL, resWOfit$coeff[1], resWOfit$coeff[2], resWOfit$NETM, sampleName=sampleNames[i])
		dev.off()

		## save data
		save(tmp, notNAs, resWfit, resWOfit, file=paste("../Results/", sampleNames[i], "/CopyNumber_and_Genotype.RData", sep=""))
	}
}
close(pb)
