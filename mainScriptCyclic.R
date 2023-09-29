require(limma)
require(parallel)
require(reshape2)
require(d3heatmap)
require(umap)
require(plot3D)
require(gplots)require(rgl)
require(plot3D)
require(corrplot)
require(factoextra)
require(FactoMineR)
require(stringdist)
require(qualV)
require(igraph)
require(pbapply)
require(stringi)
require(matrixStats)


require(multimode)
require(multcomp)
require(Biostrings)
library(lsmeans)
require(boot)
require(Rfast)
require(sva)
require(mixtools)
require(clusterCrit)
require(e1071)
require(abind)


cpl=colorRampPalette(c("#0010303F","#0050AA9F","#10AA109F","#FFFF009F","#FFA0009F","#B500009F"), alpha=T)

#
# Collect the data from the .gpr files in a folder /gpr containing also a key.csv file describing the slides
#

alldata=chpr("gpr//", sh=T)
Coor=alldata[[1]] # 1-chip/2-zone/3-diag/4-channel/5-row/6-column
iD=alldata[[2]]   # ID of spots
Res=alldata[[3]]  # Results - locally normalized data (the background is reconstituted using the duplicate diff. and subtracted using normexp)
rownames(Res)=iD  # 
Fl=alldata[[4]]   # Flags data
Nam=alldata[[5]]  # Names data (not used)

coix=as.data.frame(t(sapply(alldata[[6]], function(x){as.character(x)})), stringsAsFactors=FALSE) # 1.Chip;2.Zone;3.Channel;4.Patient ID _ Diagnosis;5.Positive Column; 6.Negative Column
FLdf=as.data.frame(Fl[1:ncol(Res)])
FL=rowSums(FLdf)==0
Cop=Coor[as.double(coix[,5])]

WD=Res[FL,]
WDa=aggregate.data.frame(WD, by=list(rownames(WD)), FUN=mean)

WDam=data.frame(log10(WDa[,2:ncol(WDa)]), stringsAsFactors = FALSE, row.names = WDa[,1])
# normalize for amino acid composition dependent binding (e.g. non-specific stickiness of charged amino acids)
D=pepnorm(WDam)
# normalize arrays 
Dn=normalizeCyclicLoess(D, method="affy", span=.15, iterations=3) 
pp=rownames(Dn)
l=length(pp)

fdn=rep(1,5)
fdn[1:2]=2
fdn[3:4]=3
fdn=factor(fdn)
calls=limfit_(fdn,Dn, p=0.05)