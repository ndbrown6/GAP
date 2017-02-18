'getCNrank' <- function(cnV, c_CL)
{
	c_CL[4,] = c_CL[2,]-cnV
	kk = which.min(abs(c_CL[4,]))
	c_CL[5,kk] = 10
	if (kk==1) {
		if ((c_CL[4,kk]<0)&&(c_CL[3,kk+1]<abs(c_CL[4,kk]))) {
			c_CL[5,kk+1] = 9
		} else {
			c_CL[5,kk+1] = 8
		}
	} else {
		if (kk==9) {
			if ((c_CL[4,kk]>0)&&(c_CL[3,kk]<c_CL[4,kk])) {
				c_CL[5,kk-1] = 9
			} else {
				c_CL[5,kk-1] = 8
			}
		} else {
			if (c_CL[4,kk]>=0) {
				if (c_CL[3,kk]<c_CL[4,kk]) {
					c_CL[5,kk-1] = 9
					c_CL[5,kk+1] = 8
				} else {
					c_CL[5,kk-1] = 8
					c_CL[5,kk+1] = 8
				}
			} else {
				if (c_CL[3,kk+1]<abs(c_CL[4,kk])) {
					c_CL[5,kk+1] = 9
					c_CL[5,kk-1] = 8
				} else {
					c_CL[5,kk+1] = 8
					c_CL[5,kk-1] = 8
				}
			}
		}
	}
	kk = which(c_CL[5,]==1)
	c_CL[5,kk] = order(abs(c_CL[4,kk]),decreasing=TRUE)+1
	tmp = c_CL[5,]
	return(invisible(tmp))
}
