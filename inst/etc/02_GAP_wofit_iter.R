library('GAP')
GAP_Params = read.csv(file="../GAP_Params/GAP_Params.csv", header=TRUE, sep=",", row.names=1, stringsAsFactors=FALSE)
Sample_Names = rownames(GAP_Params)

i = 1

print(Sample_Names[i])
load(paste("../Res/", Sample_Names[i], "/CopyNumber_and_Genotype.RData", sep=""))

# p_BAF, q_LRR, Delta
p_BAF = as.numeric(GAP_Params[Sample_Names[i],"p_BAF"])
q_LRR = as.numeric(GAP_Params[Sample_Names[i],"q_LRR"])
Delta = c(2, as.numeric(GAP_Params[Sample_Names[i],"X2copy_LRR"]))
	
## fit GAP manually
resWOfit = CopyNumber_and_GenotypeWOfitting(tmp[notNAs,GAP::.GapEnv$c_lrr], tmp[notNAs,GAP::.GapEnv$c_baf], tmp[notNAs,GAP::.GapEnv$c_len], showGAPTemplate=FALSE, p_BAF, q_LRR, Delta)

## plot preliminary GAP pattern recognition
pdf(file=paste("../Res/", Sample_Names[i], "/GAP_Preliminary_pattern_wofit.pdf", sep=""))
GAPplot_preliminary(tmp, resWOfit$Centr_CL, resWOfit$coeff[1], resWOfit$coeff[2], resWOfit$NETM, sampleName=Sample_Names[i])
dev.off()

## save data
save(tmp, notNAs, resWfit, resWOfit, file=paste("../Res/", Sample_Names[i], "/CopyNumber_and_Genotype.RData", sep=""))
