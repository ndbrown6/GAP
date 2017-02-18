'smoothNET' <- function(NET, weightC)
{
	NR = dim(NET)[1]
	NF = dim(NET)[2]
	for (k in 1:2) {
		NET1 = NET
		for (j in 2:(NF-1)) {
			for (p in 2:(NR-1)) {
				NET1[p,j] = (sum(NET[(p-1):(p+1),(j-1):(j+1)])+weightC*NET[p,j])/(9+weightC)
				if ((NET[p,j]>50)&&(NET1[p,j]<50)) {
					NET1[p,j] = NET[p,j]
				}
			}
			NET1[1,j] = sum(NET[1:2,(j-1):(j+1)])/6
			NET1[NR,j] = sum(NET[(NR-1):NR,(j-1):(j+1)])/6
		}
		for (p in 2:(NR-1)) {
			NET1[p,1] = sum(NET[(p-1):(p+1),1:2])/6
			NET1[p,NF] = sum(NET[(p-1):(p+1),(NF-1):NF])/6
		}
		NET1[1,1] = (NET[1,1]+NET[2,1]+NET[1,2]+NET[2,2])/4
		NET1[NR,NF] = (NET[NR,NF]+NET[NR-1,NF]+NET[NR,NF-1]+NET[NR-1,NF-1])/4
		NET = NET1
	}
	return(invisible(NET))
}
