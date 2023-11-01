require(parallel)
require(pbapply)
require(future.apply)
require(stringdist)
require(Biostrings)
require(igraph)

aa=AA_STANDARD

BSDvr=sapply(1:20, function(i0){
    print(i0)
    proct=proc.time()
    plan("multisession", workers=16)
    ABREP0=future_sapply(5:10e6, function(i){
      paste(sample(aa, 20, replace=T), sep = "", collapse="")
    })
    closeAllConnections()
    
    probes=pbsapply(1:10e4, function(i){
      paste(sample(aa, 7, replace = T), sep = "", collapse="")
    })
    assay=sample(probes,1000)
    
    DvrsRep=sapply(1:10,function(i) sample(ABREP0,500000))
    SprsRep=sapply(1:10,function(i) sample(sample(ABREP0,50000),500000,replace = T))
    
    plan("multisession", workers=16)
    testDvr=future_apply(DvrsRep,2,function(rp){
      x=stringdistmatrix(assay,rp, method="lcs")
      rowSums(x==13)
    })
    closeAllConnections()
    rownames(testDvr)=assay
    testDvr=testDvr/1000
    
    plan("multisession", workers=16)
    testSpr=future_apply(SprsRep,2,function(rp){
      x=stringdistmatrix(assay,rp, method="lcs")
      rowSums(x==13)
    })
    closeAllConnections()
    rownames(testSpr)=assay
    testSpr=testSpr/1000
    
    ALdvr=pblapply(seq_along(testDvr[,1]),function(i1){
            x=sapply(seq_along(testDvr[,1]),function(i2){
              vCCV(testDvr[c(i1,i2),])
            })
            assay[x>2.4]
    })
    Gdvr=adjL2G(ALdvr)
    
    ALspr=pblapply(seq_along(testSpr[,1]),function(i1){
      x=sapply(seq_along(testSpr[,1]),function(i2){
        vCCV(testSpr[c(i1,i2),])
      })
      assay[x>2.4]
    })
    Gspr=adjL2G(ALspr)
    print(proc.time()-proct)
    c(MeanDvr=10^mean(log10(testDvr+1e-4)), MeanSpr=10^mean(log10(testSpr+1e-4)), DensityDvr=graph.density(Gdvr),DensitySpr=graph.density(Gspr))

})

boxplot(t(BSDvr[1:2,]), names=c("Diverse", "Sparse"), ylab="Geometric mean of intensity", ylim=c(0,0.017))
boxplot(t(BSDvr[3:4,]), names=c("Diverse", "Sparse"), ylab="Graph density", ylim=c(0,0.022))
