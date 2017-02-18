'getCNwithPrecision' <- function(vectCN, Centr_CL, Precision)
{
	ss = which(!is.na(vectCN))
	vectCN_init = vectCN
	vectCN = vectCN_init[ss]
	if (length(vectCN)>0) {
		c_CL = table(Centr_CL[5,])
		c_CL = rbind(seq(0,8),as.numeric(rownames(c_CL)))
		c_CL = rbind(c_CL,c_CL,seq(0,8))
		c_CL[3,2:dim(c_CL)[[2]]] = c_CL[2,2:dim(c_CL)[[2]]]-c_CL[2,1:(dim(c_CL)[[2]]-1)]
		c_CL[3,1] = 1
		c_CL[3,] = round(c_CL[3,]*Precision,4)
		c_CL[5,] = 1
		Centr_CL = cbind(Centr_CL,Centr_CL[,which(Centr_CL[1,]-Centr_CL[2,]==0)])
		Centr_CL[2,26:34] = 0
		Centr_CL[3,26:34] = Centr_CL[1,26:34]+30
		Centr_CL[4,26:34] = 1
		Centr_CL = Centr_CL[,order(Centr_CL[3,],decreasing=FALSE)]
		vectRang = vectCN*0
		for (k in 1:length(vectCN)) {
			c_CL[5,] = 1
			c_CL[5,] = getCNrank(vectCN[k],c_CL)
			tt = which(c_CL[5,]==10)
			tt1 = which(c_CL[5,]==9)
			if (length(tt1)==1) {
				kk = round(abs(c_CL[4,tt]-vectCN[k])/abs(c_CL[4,tt]-c_CL[4,tt1]),1)
				vectRang[k] = c_CL[1,tt]*(1-kk)+c_CL[1,tt1]*kk
			} else {
				vectRang[k] = c_CL[1,tt]
			}
		}
		vectCN_init[ss] = vectRang
	}
	return(invisible(vectCN_init))
}
