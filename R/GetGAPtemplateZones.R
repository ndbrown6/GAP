'GetGAPtemplateZones' <- function(NET,Centr_CL)
{
	NF = dim(NET)[2]
	NR = dim(NET)[1]
	SF = as.numeric(colnames(NET))
	SR = as.numeric(rownames(NET))
	step_SF = SF[2]-SF[1]
	step_SR = SR[2]-SR[1]
	NET_C = GAPonNETcoordinates(Centr_CL,SF,SR)
	NET_tmp = GAPonNET_template(Centr_CL,NET_C,SF)
	p = 1
	while ((p<=dim(NET_tmp)[1])&&(NET_tmp[p,1]>dim(NET)[1])) {
		p = p+1
	}
	if ((p>1)&&(p<=dim(NET_tmp)[1])) {
		tt = matrix(0,NET_tmp[p-1,1]-NET_tmp[p,1]+1,dim(NET)[2])
		rownames(tt) = seq(SR[NR]+step_SR,SR[NR]+dim(tt)[1]*step_SR,step_SR)
		NET = rbind(NET,tt)
	}
	while ((p<=dim(NET_tmp)[1])&&(NET_tmp[p,1]>0)) {
		p = p + 1
	}
	if (p<dim(NET_tmp)[1]) {
		tt = matrix(0,(-NET_tmp[p,1]+1),dim(NET)[2])
		rownames(tt) = seq(SR[1]-step_SR*dim(tt)[1],SR[1]-step_SR,step_SR)
		NET = rbind(tt,NET)
		NET_tmp[,1] = NET_tmp[,1]-NET_tmp[p,1]+1
	}
	tt = which(NET_tmp[,1]<0)
	if (length(tt)>0) {
		NET_tmp = NET_tmp[-tt,]
	}
	tt = which(NET_tmp[,1]>dim(NET)[1])
	if (length(tt)>0) {
		NET_tmp = NET_tmp[-tt,]
	}
	SR = as.numeric(rownames(NET))
	NR = length(SR)
	tmp = markNETtemplate(NET,NET_tmp)
	return(invisible(tmp))
}
