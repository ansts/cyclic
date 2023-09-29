require(igraph)
require(matrixStats)
require(reshape2)
require(parallel)
require(stringdist)


load("allcallscliques")
load("allcallscliq_grouped")
load("allpeps")

aclAL=adjL(allcalls, t=7)
Gcls7=adjL2G(aclAL)
Gcls7=simplify(Gcls7)
grps=rep(0,length(allcalls))
names(grps)=allcalls
grps[ppcqsigADhi]="1"
grps[ppcqsigADlo]="2"
grps[ppcqsigFTDhi]="3"
grps[ppcqsigFTDlo]="4"

vGcls7=names(V(Gcls7))

Gcls7=set_vertex_attr(Gcls7, name="group", value=grps[vGcls7] )
write.graph(Gcls7, format="graphml",file="Gcls7.graphml")

appAL=adjL(pp, t=7)
appAM=as.matrix(stringdistmatrix(pp, method = "lcs"))
appAM=appAM*(appAM<7)
appAM[appAM>0]=1/appAM[appAM>0]
rownames(appAM)=pp
colnames(appAM)=pp
#Gpp7=adjL2G(appAL)
Gpp7=graph_from_adjacency_matrix(appAM, mode="undirected",weighted = T)
Gpp7=simplify(Gpp7)  
grpps=rep(0,length(pp))
names(grpps)=pp
grpps[ppcqsigADhi]="1"
grpps[ppcqsigADlo]="2"
grpps[ppcqsigFTDhi]="3"
grpps[ppcqsigFTDlo]="4"

vGpp7=names(V(Gpp7))

Gpp7=set_vertex_attr(Gpp7, name="group", value=grpps[vGpp7])
write.graph(Gpp7, format="graphml",file="Gpp7.graphml")
  
assortativity(Gpp7, grpps==4)


# Peptides to print on chips

ADhisel=sample(ppcqsigADhi, 43)
ADlosel=sample(ppcqsigADlo, 18)
FTDhisel=sample(ppcqsigFTDhi, 49)
FTDlosel=sample(ppcqsigFTDlo, 20)

pepsel=c(ADhisel,ADlosel,FTDhisel,FTDlosel)
write.csv(pepsel, file="peptidesPashov_linear_cyclic.csv")
