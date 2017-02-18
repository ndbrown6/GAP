'getRecCHR' <- function(SNP_annot)
{
	C_Chr = 2
	tt = table(SNP_annot[,C_Chr],SNP_annot[,C_Chr+1])
	CHR = c(0)
	tCHR = c(0)
	for (k in 1:dim(tt)[1]) {
		CHR = c(CHR,CHR[2*k-1]+tt[k,1])
		CHR = c(CHR,CHR[2*k]+tt[k,2])
		tCHR = c(tCHR,as.numeric(rownames(tt))[k])
		tCHR = c(tCHR,as.numeric(rownames(tt))[k]+0.5)
	}
	names(CHR) = tCHR
	return(invisible(CHR))
}
