'Centr_CL_annot' <- function(CN)
{
	CC = matrix(0,2,2)
	k = 2
	for (cn in 2:CN) {
		CC = cbind(CC,matrix(0,2,k))
		if (cn%%2==1) {
			k = k+1
		}
	}
	for (cn in 1:(CN+1)) {
		if (cn-1<=2) {
			cn_CL = cn-1
		}
		if ((cn-1)>2)	{
			cn_CL = 2
			for (k in 2:(cn-1)) {
				for (nbc in trunc(k/2+1):k) {
					cn_CL = cn_CL+1
				}
			}
			cn_CL = cn_CL-1
		}
		for (nbc in 1:trunc((cn-1)/2+1)) {
			CC[1,cn_CL+nbc] = cn-1
			CC[2,cn_CL+nbc] = trunc((cn-1+1)/2)+nbc-1
		}
	}
	CC = rbind(CC,seq(1,dim(CC)[2]))
	CC = rbind(CC,matrix(0,2,dim(CC)[2]))
	rownames(CC) = c("Copy","Major allele","Index","BAF","LRR")
	return(invisible(CC))
}
