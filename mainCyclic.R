# Packages -----

require(limma)
require(parallel)
require(reshape2)
require(d3heatmap)
require(uwot)
require(plot3D)
require(gplots)
require(rgl)
require(plot3D)
require(corrplot)
require(factoextra)
require(FactoMineR)
require(stringdist)
require(qualV)
require(pbapply)
require(stringi)
require(matrixStats)
require(Biostrings)
require(clusterCrit)
require(e1071)
require(abind)
require(igraph)
require(eulerr)
require(future.apply)
require(chisq.posthoc.test)
require(alluvial)

# ... from a previous work, using a joint reactivity graph based on 4 sera from
# patients with AD and 4 sera from patients with FTD compared to pooled human
# IgM, an analysis of the cliques determined cliques with significant over/under
# expression of reactivity in AD/FTD. The actual number of peptides are smaller
# to get bellow the limit for the PEPPerPRINT 5 field arrays - 130 peptides

load("allcallscliq_grouped")

## Peptides to print on chips ----

ADhisel=sample(ppcqsigADhi, 43)
ADlosel=sample(ppcqsigADlo, 18)
FTDhisel=sample(ppcqsigFTDhi, 49)
FTDlosel=sample(ppcqsigFTDlo, 20)

pepsel=c(ADhisel,ADlosel,FTDhisel,FTDlosel)
pepgrp=(pepsel %in% ADhisel)*1+(pepsel %in% ADlosel)*2+(pepsel %in% FTDhisel)*3+(pepsel %in% FTDlosel)*4
names(pepgrp)=pepsel
gr=c("ADhisel","ADlosel","FTDhisel","FTDlosel")
names(gr)=1:4
pepgrp=gr[pepgrp]
names(pepgrp)=pepsel

write.csv(pepsel, file="peptidesPashov_linear_cyclic.csv")




# Analysis of the IgG and IgM Binding ----

cpl=colorRampPalette(c("#000000","#0050AA9F","#10AA109F","#FFFF009F","#FFA0009F","#B50000"), alpha=T)

## Collect the data from the .gpr files in a folder /gpr containing also a key.csv file describing the slides ----

pth="gprfiles//"
dgn=cbind(1:16,c("A","K","F","X","K","A","X","F","A","K","F","X","K","A","X","F"))

lf=list.files(path = p)
lf=lf[lf!="key.csv"]
fnms=t(apply(sapply(strsplit(lf,split=list("_|-")), unlist),2,function(l){
  x1=paste(l[4:5], collapse="",sep="")
  x2=stri_extract(l[6], regex="\\d+")
  x3=paste("p",x2,"_",substr(l[1],1,1),sep="")
  if (l[4]=="sera")  d=dgn[dgn[,1]==x2,2] else  d="B"
  return(c(x3,1,d))
}))

fnms=rbind(fnms,fnms)
key=cbind(c(lf,lf),fnms,c(rep("G",length(lf)),rep("R",length(lf))))
colnames(key)=c("file","pat","block","diag","ch")
write.csv(key, file=paste(p,"key.csv", sep=""))

# chpr1 extracts the data from the gpr files and organizes it subtracting the background
alldata=chpr1(pth, sh=T)
Coor=alldata[[1]] # 1-chip/2-zone/3-diag/4-channel/5-row/6-column
iD=alldata[[2]]   # ID of spots
Res=alldata[[3]]  # Results - locally normalized data (the background is reconstituted using the duplicate diff. and subtracted using normexp)
rownames(Res)=iD  # 
Fl=alldata[[4]]   # Flags data

coix=as.data.frame(t(sapply(alldata[[5]], function(x){as.character(x)})), stringsAsFactors=FALSE) # 1.Chip;2.Channel;3.Patient ID _ Diagnosis
FLdf=as.data.frame(Fl[1:ncol(Res)])
FL=rowSums(FLdf)==0

# WD - data filtered for flags from the densitormetry then aggregated to average
# the duplicate peptide spots
WD=Res[FL,]
WDa=aggregate.data.frame(WD, by=list(rownames(WD)), FUN=mean)
WDam=data.frame(log10(WDa[,2:ncol(WDa)]), stringsAsFactors = FALSE, row.names = WDa[,1])

#### normalize for aminoacid composition dependent binding (e.g. non-specific stickiness of charged amino acids) ----
D=pepnorm(WDam)

#### normalize between arrays within topology/isotype groups ----
Dn=array(0, dim=dim(WDam), dimnames=dimnames(WDam))
prms=as.data.frame(t(sapply(colnames(WDam), function(l) unlist(strsplit(l, split="_")))))
colnames(prms)=c("Sample","Topology","Diagnosis","Channel")
fld=3
Dn[,prms$Channel=="G"&prms$Topology=="l"]=normalizeCyclicLoess(D[,prms$Channel=="G"&prms$Topology=="l"], method="affy", iterations = fld) 
Dn[,prms$Channel=="R"&prms$Topology=="l"]=normalizeCyclicLoess(D[,prms$Channel=="R"&prms$Topology=="l"], method="affy", iterations = fld)
Dn[,prms$Channel=="G"&prms$Topology=="c"]=normalizeCyclicLoess(D[,prms$Channel=="G"&prms$Topology=="c"], method="affy", iterations = fld) 
Dn[,prms$Channel=="R"&prms$Topology=="c"]=normalizeCyclicLoess(D[,prms$Channel=="R"&prms$Topology=="c"], method="affy", iterations = fld)
pp=rownames(Dn)
l=length(pp)





# ReaGraph ----

mn=rowMeans(Dn)
std=rowSds(Dn)

### Graphs by Topology, Isotype and Diagnosis ----

allgr=expand.grid(Topology=c("l","c"), Diagnosis=c("A","F","X","K"), Channel=c("G","R")) # grouping categories
grnames=apply(allgr,1,paste, collapse="")
rownames(allgr)=grnames

# the actual graphs are in the list Gi
Gi=lapply(grnames, function(gr){
  g=ReGr(Dn[,prms$Channel  ==as.character(allgr[grnames==gr,"Channel"]) &
             prms$Topology ==as.character(allgr[grnames==gr,"Topology"]) &
             prms$Diagnosis==as.character(allgr[grnames==gr,"Diagnosis"]) ])
  g=set_edge_attr(g, name="Channel", value=as.character(allgr[grnames==gr,"Channel"]) )
  g=set_edge_attr(g, name="Topology", value=as.character(allgr[grnames==gr,"Topology"]))
  g=set_edge_attr(g, name="Diagnosis", value=as.character(allgr[grnames==gr,"Diagnosis"]))
  g=set_edge_attr(g, name="SourceGraph", value=gr)
  g=set_vertex_attr(g, name="Group", value=pepgrp)
  g=set_vertex_attr(g, name="Mean", value=mn)
  g=set_vertex_attr(g, name="SD", value=std)
  g=set_graph_attr(g,name="name", value=gr)
  write.graph(g, format = "graphml", file=paste(gr,".graphml", collapse="",sep=""))
  return(g)
})
names(Gi)=sapply(Gi, function(g) graph_attr(g)$name)

### Graphs order, size, density, etc. params ----
gf=apply(allgr[,c(1,3)],1,function(l) paste(l,collapse=""))                   
gf=sub("G", "IgM", gf)
gf=sub("R", "IgG", gf)
gf=sub("l", "linear ", gf)
gf=sub("c", "cyclic ", gf)
gf=as.factor(gf)
sapply(Gi, components)

ordGi=vcount(Gi[[1]])
sizGi=sapply(Gi, function(g) ecount(g))
kruskal.test(sizGi, gf)

#### Figure 4 ----
dnsGi=sapply(Gi,function(g) graph.density(g))
kruskal.test(dnsGi, gf)
gff=gf
gff=sub("linear","l",gff)
gff=sub("cyclic","c",gff)
gff=sub("IgM","M",gff)
gff=sub("IgG","G",gff)
dunn.test::dunn.test(dnsGi, gff)
boxplot(dnsGi~gf, xlab="Topology/Isotype", ylab="Graph density",  par(bty="n"))
legend("bottomleft", legend="Dunn test p<0.05", bty = "n")
lines(c(1,3),c(0.382,0.382))
lines(c(3,4),c(0.381,0.381))

meanGi=t(aggregate(t(Dn), by=list(prms$Topology,prms$Channel), mean)[,-(1:2)])
colnames(meanGi)=c("cyclic_IgM","linear_IgM","cyclic_IgG","linear_IgG")
x=melt(meanGi)
dunn.test::dunn.test(x$value, g=x$Var2)
boxplot(x$value~x$Var2, xlab="Topology/Isotype", ylab="Mean intensity",  par(bty="n"))
legend("topleft", legend="Dunn test p<0.05/ all differences significant", bty = "n")

aggregate(dnsGi, by=list(gf), "c")

## Graph reunited ----

Gall=do.call(igraph::union,Gi)

j=grepl("weight",edge_attr_names(Gall))  
x=do.call(cbind,edge_attr(Gall)[j])
x=apply(x,1,function(y) sum(y[!is.na(y)]))

Gall=set_edge_attr(Gall, name="weight", value=x)

j=grepl("Channel",edge_attr_names(Gall))  
x=do.call(cbind,edge_attr(Gall)[j])
x=t(apply(x,1,function(y) {
  t0=rep(0,2)
  names(t0)=c("R","G")
  tx=table(y[!is.na(y)])
  t0[names(tx)]=tx
  return(t0)
}))

Gall=set_edge_attr(Gall, name="ChR", value=x[,1])
Gall=set_edge_attr(Gall, name="ChG", value=x[,2])

j=grepl("Topology",edge_attr_names(Gall))  
x=do.call(cbind,edge_attr(Gall)[j])
x=t(apply(x,1,function(y) {
  t0=rep(0,2)
  names(t0)=c("l","c")
  tx=table(y[!is.na(y)])
  t0[names(tx)]=tx
  return(t0)
}))

Gall=set_edge_attr(Gall, name="TopL", value=x[,1])
Gall=set_edge_attr(Gall, name="TopC", value=x[,2])

j=grepl("Diagnosis",edge_attr_names(Gall))  
x=do.call(cbind,edge_attr(Gall)[j])
x=t(apply(x,1,function(y) {
  t0=rep(0,4)
  names(t0)=c("A","F","X","K")
  tx=table(y[!is.na(y)])
  t0[names(tx)]=tx
  return(t0)
}))

Gall=set_edge_attr(Gall, name="DiaA", value=x[,1])
Gall=set_edge_attr(Gall, name="DiaF", value=x[,2])
Gall=set_edge_attr(Gall, name="DiaX", value=x[,3])
Gall=set_edge_attr(Gall, name="DiaK", value=x[,4])

j=grepl("SourceGraph",edge_attr_names(Gall))  
x=do.call(cbind,edge_attr(Gall)[j])
x=t(apply(x,1,function(y) {
  t0=rep(0,16)
  names(t0)=grnames
  tx=table(y[!is.na(y)])
  t0[names(tx)]=tx
  return(t0)
}))

Gall=set_vertex_attr(Gall, name="Group", value=vertex_attr(Gall)$Group_1)
Gall=set_vertex_attr(Gall, name="Mean", value=vertex_attr(Gall)$Mean_1)
Gall=set_vertex_attr(Gall, name="SD", value=vertex_attr(Gall)$SD_1)

for (i in 1:16) {
  Gall=delete_vertex_attr(Gall,name = paste("Group",i,sep="_"))
  Gall=delete_vertex_attr(Gall,name = paste("Mean",i,sep="_"))
  Gall=delete_vertex_attr(Gall,name = paste("SD",i,sep="_"))
  nmw=paste("weight",i,sep="_")
  att=edge_attr(Gall)[edge_attr_names(Gall)==nmw][[1]]
  att[is.na(att)]=0
  Gall=set_edge_attr(Gall, name=paste("w",graph_attr(Gall)[[1]], sep="_"), value=att)
  Gall=delete_edge_attr(Gall, name=nmw)
  Gall=delete_edge_attr(Gall,name = paste("Topology",i,sep="_"))
  Gall=delete_edge_attr(Gall,name = paste("Diagnosis",i,sep="_"))
  Gall=delete_edge_attr(Gall,name = paste("SourceGraph",i,sep="_"))
  Gall=delete_edge_attr(Gall,name = paste("Channel",i,sep="_"))
  Gall=delete_graph_attr(Gall, name=paste("name",i,sep="_"))
}
Gall=subgraph.edges(Gall, E(Gall)[edge_attr(Gall)$weight>0])

#### Figure 3 ----
# the graphml file is passed to Gephi for visualization with x and y
# attributes calculated here based on the Laplasian eigenvectors
Gallsparse=subgraph.edges(Gall, E(Gall)[edge_attr(Gall)$weight>17.5])
LGall=embed_laplacian_matrix(Gallsparse, no=129)
XY=umap(LGall$X[,114:129], n_neighbors = 80)
Gallsparse=set_vertex_attr(Gallsparse, name="x", value=800*XY[,1])
Gallsparse=set_vertex_attr(Gallsparse, name="y", value=800*XY[,2])
write.graph(Gallsparse, format = "graphml", file="Gallsparse.graphml")

# Fig 3 - modularities
modGsparse=modularity(Gallsparse, as.numeric(as.factor(vertex_attr(Gallsparse)$Group)))
BSmodsparse=sapply(1:1000, function(i){
  x=sample(as.numeric(as.factor(vertex_attr(Gallsparse)$Group)))
  modularity(Gallsparse,x)
})
fBSmodsparse=ecdf(BSmodsparse)
fBSmodsparse(modGsparse)

qntlmodbyGr=pbsapply(unique(vertex_attr(Gall)$Group), function(gr){
  mo=modularity(Gallsparse, as.numeric(as.factor(vertex_attr(Gallsparse)$Group==gr)))
  BS=sapply(1:1000, function(i){
    x=sample(as.numeric(as.factor(vertex_attr(Gallsparse)$Group==gr)))
    modularity(Gallsparse,x)
  })
  fBS=ecdf(BS)
  fBS(mo)
})

### Distribution of common edges in structure graphs ----

vgrep=Vectorize(grep,vectorize.args = "pattern")

allgruStr=unique(allgr[,c(1,3)])

ya=do.call(cbind,edge_attr(Gall)[vgrep(grnames, edge_attr_names(Gall))])
ya=apply(allgruStr, 1, function(l) 
  rowSums(ya[,which(as.character(allgr[,1])==l[1] & as.character(allgr[,3])==l[2])]))
n=colnames(ya)
n1=paste(substr(n,1,1),substr(n,3,3), sep="")    
n1=sub("G","_IgM",n1)
n1=sub("R","_IgG",n1)
colnames(ya)=n1
names(colnames(ya))=n

Ya=sapply(c(2,5,8),function(j) {
  yy=ya*(ya>=j & ya<(j+2))
  x=aggregate(yy, by=list(apply(yy>0,1,function(l) paste((colnames(ya)[l]), collapse="_"))), "sum")
  xn=rowSums(x[,2:5])
  x[,1][x[,1]==""]="none"
  names(xn)=x[,1]
  return(xn)
})
nYa=unique(unlist(lapply(Ya,names)))
mYa=array(0,dim=c(length(nYa),3))
rownames(mYa)=nYa
for (i in 1:3) {
  mYa[names(Ya[[i]]),i]=Ya[[i]]
}
mYa=mYa[rownames(mYa)!="none",]

#### bootstrap the overlap distributions ----
plan("multisession", workers=10)
BSa=future_lapply(1:1000, function(co){
  g=Gall
  for (ea in edge_attr_names(g)[vgrep(grnames, edge_attr_names(g))]){
    edge_attr(g)[[ea]]=sample(edge_attr(g)[[ea]])
  }
  y=do.call(cbind,edge_attr(g)[vgrep(grnames, edge_attr_names(g))])
  y=apply(allgruStr, 1, function(l) 
    rowSums(y[,which(as.character(allgr[,1])==l[1] & as.character(allgr[,3])==l[2])]))
  n=colnames(y)
  n1=paste(substr(n,1,1),substr(n,3,3), sep="")    #c("IgM-l","IgM-c","IgG-l","IgG-c")
  n1=sub("G","_IgM",n1)
  n1=sub("R","_IgG",n1)
  colnames(y)=n1
  names(colnames(y))=n
  Y=sapply(c(2,5,8),function(j) {
    yy=y*(y>j & y<(j+2))
    x=aggregate(yy, by=list(apply(yy>0,1,function(l) paste((colnames(y)[l]), collapse="_"))), "sum")
    xn=rowSums(x[,2:5])
    x[,1][x[,1]==""]="none"
    names(xn)=x[,1]
    return(xn)
  })
  nY=unique(unlist(lapply(Y,names)))
  mY=array(0,dim=c(length(nY),3))
  rownames(mY)=nY
  for (i in 1:3) {
    mY[names(Y[[i]]),i]=Y[[i]]
  }
  mY=mY[rownames(mY)!="none",]
  return(mY)
})
closeAllConnections()

x=BSa
x=unlist(x)
x=array(x, dim=c(dim(BSa[[1]]),1000), dimnames = list(dimnames(BSa[[1]])[[1]],(1:3)*2,1:1000))
xmax=apply(x,c(1,2), quantile, 0.95)
xmin=apply(x,c(1,2), quantile, 0.05)
xmn=apply(x,c(1,2), quantile, 0.5)

#### Figure 5 ----
par0=par()
par(mai=par()$mai+c(1.2,0,0,0))
for (i in 1:4){
  barplot(log10(mYa[,i]+1),ylim=c(0,5), las=2, xlab="", yaxt="n", bty="n", 
        main=(paste("Weight threshold = ",seq(2,8,2)[i],collapse = "",sep="")))
  par(new=T)
  plot((0:14)+0.5,log10(xmax[rownames(mYa),i]+1),ylim=c(0,5), ty="l",col=2, 
        xlab="", ylab="",xaxt="n", xlim=c(0,16), bty="n")
  par(new=T)
  plot((0:14)+0.5,log10(xmin[rownames(mYa),i]+1),ylim=c(0,5), ty="l",col=3, 
       main=colnames(BSa[[1]][i]), ylab="Number of edges", xlab="", 
       xaxt="n", xlim=c(0,16), bty="n")
  par(new=F)
}
par(mai=par()$mai-c(1.2,0,0,0))

par0=par()
par(mai=par()$mai+c(1.2,0,0,0))
psc=min(xmn[xmn>0])
Y=(mYa+psc)/(xmn+psc)
cl=1*(mYa>xmax)+
   2*(mYa<xmin)+1
for (i in 1:3){
  j=order(Y[,i],decreasing = T)
  barplot(Y[j,i],las=2, xlab="",col=cl[j,i])
}
par(mai=par()$mai-c(1.2,0,0,0))

Yluv=Y[order(nchar(rownames(Y))),]
colnames(Yluv)=c("Low","Medium","High")
Yluv=melt(Yluv)
colnames(Yluv)=c("Subgraph","Crossreactivity","Ttl(E)Weight")
ylord=apply(Yluv,1,function(l) paste(l[1],l[2],sep="_"))
cluv=cl
colnames(cluv)=c("Low","Medium","High")
cluv=melt(cluv)
clord=apply(cluv,1,function(l) paste(l[1],l[2],sep="_"))
rownames(cluv)=clord
cluv=cluv[ylord,]
alluvial(Yluv[,1:2], freq = Yluv[,3], col=cluv[,3], border = cluv[,3], cex=0.9)





# Cyclization and biomarker relevance ----
## Laplacian embedding ----
eigLapEmb=lapply(seq_along(Gi), function(i){
  g=Gi[[i]]
  eL=embed_laplacian_matrix(g, no=vcount(g)-1)
  l=length(eL$D[eL$D>0])
  n=dim_select(rev(max(eL$D[1:l])-eL$D[1:l]))-1
  Gi[[i]]=set_vertex_attr(Gi[[i]], name="x", value=eL$X[,l])
  Gi[[i]]=set_vertex_attr(Gi[[i]], name="y", value=eL$X[,l-1])
  write.graph(g, format = "graphml", file=paste(graph_attr(g)$name,".graphml", collapse="",sep=""))
  eL$X[,(l-n):l]
})

# splitting the data
grtop=paste(prms$Topology,prms$Channel, sep="_")
grtop=sub("G","IgM",grtop)
grtop=sub("R","IgG",grtop)
grtonms=unique(grtop)
Mcs=lapply(grtonms, function(ti){
  m=Dn[,grtop==ti]
})
names(Mcs)=grtonms
diags=sapply(Mcs,function(m){
  x=colnames(m)
  unlist(stri_extract_all(x, regex="(?<=[c,l]_)\\w"))
})
diags=as.factor(diags[,1])

# RFE 
clssf_16=lapply(names(Mcs), function(n){
  m=Mcs[[n]]
  x=rfe(m,diags)
  plotMDS(x[[1]], col=as.numeric(diags), main=n, xlab="D1", ylab="D2")
  dg=as.character(unique(diags))
  ccmx=sapply(dg, function(d1){
    sapply(dg, function(d2){
      if (d1!=d2) clucri(x[[1]][,diags %in% c(d1,d2)], diags[diags %in% c(d1,d2)]) else 0
    })
  })
  return(list(max(x[[3]]), nrow(x[[1]]), rownames(x[[1]]),ccmx))
})


clucris=sapply(1:4, function(i){
  cc=clssf_16[[i]][[4]]
  cc[lower.tri(cc)]
})
colnames(clucris)=names(Mcs)
X=melt(clucris)
kruskal.test(X[,3]~X[,2])
dunn.test::dunn.test(X[,3], g=X[,2])

lmbig=lm(data=X, value~Var2)
almbig=aov(lmbig)
TukeyHSD(almbig, method="bh")

ctop=grepl("c_",X$Var2)
wilcox.test(X$value~ctop)

isot=grepl("M",X$Var2)
wilcox.test(X$value~isot)

lmciso=lm(X$value~ctop*isot)
anova(lmciso)

wilcox.test(X$value~ctop, subset=!isot)

boxplot(X$value~X$Var2, xlab="Subgraph", ylab="Sum of clustering criterion", bty="n")
lines(c(1,4), c(26,26))
lines(c(2,4), c(25,25))
legend("bottomleft", "Kruskal-Wallis, p=0.03, Dunn p<0.05", bty="n")

# Using Kruskal-Wallis and Wilcoxon tests, the sum of the clustering criteria
# for the best classifiers depend on the mostly on
# isotyope and not on topology. ANOVA does not show interaction between isotype
# and topology. Although the more sensitive IgG assay has lower performance on
# cyclic peptides, the difference is not significant (p=0.24)


x=unlist(lapply(clssf_16,function(l) l[[3]]))
xu=unique(x)
x1=lapply(clssf_16,function(l) l[[3]])
X=cbind(xu %in% x1[[1]], xu %in% x1[[2]],xu %in% x1[[3]],xu %in% x1[[4]])
colnames(X)=names(Mcs)
euX=euler(X, , shape="ellipse", control=list(extraopt=T))
plot(euX, quantities=T, cex=2)

# Among the selected features, those based on c_IgM and l_IgM had the least
# overlap (8) while those of l_IgG and c_IgG - the most (17) (Jaccard distance 
# confirms it but the differentces are not significant). The sum of the
# clustering criterion is highest for c_IgM and lowest for l_IgG.

clssf_16_DM=sapply(clssf_16,function(l1) {           # Jaccard distance
          sapply(clssf_16, function(l2) {
              if (all(l1[[3]]==l2[[3]])) {
                return(0)
            } else {return((length(union(l1[[3]],l2[[3]])))/(length(intersect(l1[[3]],l2[[3]]))))}
          })
 })

ij=combn(4,2)[,1:3]

clssf_16_ps=apply(ij,2,function(i1) {           # Jaccard distance
            i2= (1:4)[!(1:4) %in% i1]
                print(i2)
            n=c(length(union(clssf_16[[i1[1]]][[3]],clssf_16[[i1[2]]][[3]])),
                length(union(clssf_16[[i2[1]]][[3]],clssf_16[[i2[2]]][[3]])))
            print(n)
            x=c(length(intersect(clssf_16[[i1[1]]][[3]],clssf_16[[i1[2]]][[3]])),
               length(intersect(clssf_16[[i2[1]]][[3]],clssf_16[[i2[2]]][[3]])))
              p=prop.test(x,n,alternative="two.sided", correct = T)$p.value
              gr1=paste(names(Mcs)[i1], collapse="_")
              gr2=paste(names(Mcs)[i2], collapse="_")
            return(paste(gr1,"vs",gr2,"p=",p, collapse="_",sep="_"))
  })

plot(cmdscale(clssf_16_DM), cex=0, xlab="D1", ylab="D2")
text(cmdscale(clssf_16_DM), labels = names(Mcs))
