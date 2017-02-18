.GapEnv <- new.env()
vN = c("c_ind", "c_is", "c_if", "c_chr", "c_len",
	   "c_lrrBP", "c_bafBP", "c_lrr", "c_baf", "c_bafH",
	   "c_nC", "c_nA", "c_nS", "c_CN", "c_BA",
	   "c_conf", "c_sb", "c_CN1", "c_BA1", "c_rnk1")
for (i in 1:length(vN)) {
	assign(vN[i], i, envir=.GapEnv)
}
v2N = c("c2_ind", "c2_is", "c2_if", "c2_posS",
		"c2_posF", "c2_chr", "c2_len", "c2_cnBP",
		"c2_adBP", "c2_cn", "c2_adL", "c2_adH",
		"c2_nC", "c2_nA", "c2_nS", "c2_CN",
		"c2_BA", "c2_conf", "c2_sb", "c2_CN1",
		"c2_BA1", "c2_rnk1")
for (i in 1:length(v2N)) {
	assign(v2N[i], i, envir=.GapEnv)
}
v3N = c("c3_ind", "c3_is", "c3_if",
		"c3_posS", "c3_posF", "c3_isCN",
		"c3_ifCN", "c3_chr", "c3_len")
for (i in 1:length(v3N)) {
	assign(v3N[i], i, envir=.GapEnv)
}
cN = c("Index", "IndxSnpStart", "IndxSnpEnd",
 	   "Chr", "LenSnp", "LRR_BP",
 	   "BAF_BP", "LRR", "MajorBAF",
 	   "HomoSegm", "nCN", "nAD",
 	   "nSB", "Copy", "Mallele",
 	   "Confidence", "Subclone", "CN1",
 	   "BA1", "Rank1", "CN2",
 	   "BA2", "Rank2", "CN3",
 	   "BA3", "Rank3")
assign("ColNamesTMP", cN, envir=.GapEnv)
c2N  = c("Index", "IndxSnpStart", "IndxSnpEnd",
		 "PosSNPStart", "PosSNPEnd", "IndxCnStart",
		 "IndxCnEnd", "Chr", "LenSnp")
assign("ColNamesTMP2", c2N, envir=.GapEnv)
rm(vN, v2N, cN, c2N, i)
