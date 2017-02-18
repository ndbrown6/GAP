#==================================================
# David Brown
# dbrown@ulb.ac.be
# http://homepages.ulb.ac.be/~dbrown/
#==================================================
library('colorspace')
library('copynumber')
library('affy2sv')
library('doMC')
library('foreach')
library('mclust')
library('GAP')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
library('TxDb.Hsapiens.UCSC.hg18.knownGene')
library('org.Hs.eg.db')
library('annotate')
library('cluster')
library('stringr')
library('gdata')
library('car')
data('CytoBand')

ncores = 24
registerDoMC(ncores)
