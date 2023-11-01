
# A function which iteratively scans the set of features and finds at each cycle the feature  
# the removal of which imrpoves most the separation between the classes of interest c, then removes it 
# leaving one less feature for the next cycle. The separation is measured using a composite clustering criterion 
# consisting of the sum of Dunn's, BHgamma criteria as well as negative Connectivity.
# The criteria are first transformed to compensate for their dependence on the dimensionality 
# and to restrict them in the range between -1 and 1. The conversion functions dunnfix, connfix and bhgamfix 
# are vectorized for faster calculations. 


rfe=function(m0,c){
  require(clValid)
  require(parallel)
  require(cluster)
  require(Rfast)
  require(reshape2)
  require(matrixStats)

  if (is.logical(c)) cc=as.integer(c)+1 else cc=as.integer(c)
  #pt=proc.time()
  D0=c()
  Ds0=c()
  bad=list()
  bads=list()
  m=m0
  N0=nrow(m0)
  nm=rownames(m0)
  nc=unique(c)
  vconf=Vectorize(connfix)
  vbhgf=Vectorize(bhgamfix)
  vdunf=Vectorize(dunnfix)
  repeat {
      D=c()
      N=nrow(m)
      if (N<3) break
      D=t(sapply(1:N, function(i) {
            d1=Dist(t(m[-i,]))
            x=c(BHgamma(d1,cc), dunn(as.dist(d1),cc),connectivity(as.dist(d1),cc))
            return(x)
        }))
      NN=rep(N-1,N)
      D=cbind(vbhgf(D[,1]),vdunf(D[,2],NN),-vconf(D[,3]))
      D=rowsums(D)
      D=c(which.max(D), max(D))                       
      
      md=D[[1]]                    
      dmin=D[[2]]
      bad=c(bad,unlist(rownames(m)[md]))
      m=m[-md,]
      D0=rbind(D0,c(N-1,dmin)) 
      print(c(N,D[[2]]))  
      fnm="scanlog.txt"
      cat(c(N,D[[2]]), file=fnm, sep="\n", append=T)
  }
  
  bad=c(unlist(bad),unlist(rownames(m)))
  yr=range(D0[,2])
  plot(nrow(D0):1,D0[,2])


  if (length(D0)>0) crit=D0[which.max(D0[,2]),1] else return(list(NA,NA,NA))
  if (crit<2) crit=2
  m=m0[!(nm %in% bad[1:(length(bad)-crit)]),]
  
  return(list(m,bad,D0))
}
