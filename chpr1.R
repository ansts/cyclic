#
# Read the densitometric data from microarrays. Needs a key.csv file 
# describing the files. bgcode is the code used for prescanned arrays (background)
#

chpr1=function(p, sh=FALSE, bgcode="ps", cntr=c("G","YPYDVPDYAG","KEVPALTAVETGAT")){
  require(parallel)
  require(future.apply)
  require(stringi)
  require(limma)
  require(reshape2)
  require(matrixStats)
  
  source("bgmy.R")
  cpl=colorRampPalette(c("#001030","#0050AA","#10AA10","#FFFF00","#FFA000","#B50000"))
  lf=list.files(path = p)
  lf=lf[lf!="key.csv"]
  key=read.csv(paste(p,"key.csv", sep=""))
  N=nrow(key)
  
  # Get the files ----
  fls=lapply(lf,function(f){
    finp=paste(p,f,sep = "")
    fn=f
    fcon=file(finp)
    f2l=readLines(con=fcon, n=2)
    nlns=as.double(stri_extract(f2l[2], regex="[1-9]+"))
    fhead=readLines(con = fcon, n=nlns+2)
    x=read.delim(finp, skip=nlns+2, header=T, check.names = F, stringsAsFactors = FALSE)
    list(fn,x)
  })
  closeAllConnections()
  flns=lapply(fls,function(f){f[[1]]})
  fls=lapply(fls,function(f){f[[2]]})
  names(fls)=flns
  
  # Rearrange and summarize data ----
  bl=lapply(seq_along(key$file), function(i){
    f=fls[[key$file[i]]] 
    if (key$ch[i]=="G") 
        f=f[f$Block==key$block[i],c("Column","Row","ID","F532 Median","Flags")] 
    else 
        f=f[f$Block==key$block[i],c("Column","Row","ID","F635 Median","Flags")]
    colnames(f)=c("C", "R","ID","V","Fl")
    #key$file=stri_extract_first(key$file, regex=".+(?<=_)\\d+")                                    
    return(list(key$pat[i],key$diag[i],key$ch[i],f,key$file[i])) #1-patient/2-diag/3-channel/4-[1-Col/2-Row/3-ID/4-F635 Med/5-Flags]/5-file
    })
  
  coor=lapply(bl, function(x){
    list("File"=unlist(x[5]), "Diagnosis"=unlist(x[2]),"Channel"=unlist(x[3]),"Row"=x[4][[1]][,2], "Column"=x[4][[1]][,1])
    })                                                           #1-file/2-diag/3-channel/4-row/5-column
  if (all(unlist(lapply(coor, function(x){unlist(coor[[1]][4:5])==unlist(x[4:5])})))) print("Ordered, homogeneous Sets") else return("Error: Nonhomogeneous Sets!")
  fl=lapply(bl,function(i){i[[4]][,5]})
  ID0=lapply(bl,function(i){i[[4]][,3]})
  if( all(rowAlls(sapply(ID0, function(x) ID0[[1]]==x)))) ID0=ID0[[1]]

  pt=lapply(bl,function(i){as.character(i[[1]])})
  dg=lapply(bl,function(i){as.character(i[[2]])})
  ch=lapply(bl,function(i){as.character(i[[3]])})
  cn=paste(pt,dg,ch, sep="_")
  
  dgch=t(as.data.frame(lapply(coor, function(x){c(as.character(x$Diagnosis),as.character(x$Channel))})))
  colnames(dgch)=c("D","C")
  Res0=as.data.frame(lapply(bl, function(x){x[[4]]$V}), col.names=cn)
  
  ij=aggregate(seq_along(fl), by=list(unlist(pt),unlist(ch)), function(x) as.numeric(x))
  R0=apply(ij,1,function(j){
    i=as.numeric(j[3:4])
    x=array(i, dimnames = list(unlist(dg)[i]))
    x=Res0[,x[names(x)!="B"]]-Res0[,x[names(x)=="B"]]
    return(x)
  })
  cnm=apply(ij,1,function(j){
    ij=as.numeric(j[3:4])
    x=array(ij, dimnames = list(unlist(dg)[ij]))
    paste(j[1],names(x)[names(x)!="B"], j[2], sep="_")
  })
  colnames(R0)=cnm
  # Reusing the name fls
  fls=unlist(lapply(coor,function(x){x$File}))
  inf0 =list()
  co=1
  for (t in seq_along(fls)) {
    ci=fls[[t]]
    inf0[[co]]=list(File=ci,Channel=coor[[t]][[3]],Patient_Diag=colnames(Res0)[t]) # 1.File;2.Channel;3.Diagnosis
    co=co+1
  }
  
  # R0=apply(R0,2,function(x){x=x-min(x)+1})          
  #rownames(R0)=ID
  #plan ("multisession", workers=16)
  B0=apply(R0, 2, function(clm){
      IM=data.frame(R=coor[[1]]$Row, C=coor[[1]]$Column, V=clm)
      ID=data.frame(R=coor[[1]]$Row, C=coor[[1]]$Column, V=ID0, stringsAsFactors = FALSE)
      ID[,1]=as.double(ID[,1])
      ID[,2]=as.double(ID[,2])
      x=bgmy(IM,ID)
      x=x[order(x$R,x$C),]
      y=x[,1:2]
      x=x[,3]
      return(list(x,y))
  })                                                           
  closeAllConnections()
  
  B0xy=B0[[1]][[2]]
  B1=sapply(B0, function(l){l[[1]]})
  
  R0=R0[!(ID0 %in% cntr),]
  B1=B1[!(ID0 %in% cntr),]
  fl=lapply(fl,function(y) y[!(ID0 %in% cntr)])
  ID0=ID0[!(ID0 %in% cntr)]
  R1=sapply(seq_along(R0[1,]),function(i){
      S=R0[,i]; B=B1[,i]
      abZ=lm(S~B)
      j=cut(B,10,labels=F)
      bi=aggregate(B, by=list(j), "mean")
      z=abZ$residuals
      sdZ=sapply(1:10, function(i) x=sd(z[j==i]))
      zz=lm(log10(sdZ)~log10(bi$x))
      a=zz$coefficients[2]
      b=zz$coefficients[1]
      sdf=10^b*B^a
      msdf=10^(mean(log10(sdf)))
      res=msdf*z/sdf
      return(res-min(res)+1)
  })
  dimnames(R1)=dimnames(R0)
  # res=backgroundCorrect.matrix(R0,B1, method = "normexp", normexp.method = "mle")
  # res=sapply(seq_along(R0[1,]),function(i){
  #   S=R0[,i]; B=B1[,i]
  #   Sr=S-B
  #   abZ=lm(Sr~B)
  #   Sr=(Sr-abZ$coefficients[2]*B)
  #   S=Sr+B
  #   x=backgroundCorrect.matrix(S, t(B), method="normexp", normexp.method="mle")
  #     # 
  #     # IM=data.frame(R=coor[[1]]$Row, C=coor[[1]]$Column, V=R0[,i])
  #     # y=acast(data=IM, R~C)
  #     # par(mfrow=c(2,1))
  #     # image(y, col=cpl(100), main=colnames(R0)[i])
  #     # y=acast(data=data.frame(IM[,1:2],V=x), R~C)
  #     # image(y, col=cpl(100), main=paste(colnames(R0)[i],"correctred",sep=" "))
  #     # par(mfrow=c(1,1))
  # 
  #   return(x)
  # })

  return(list(coor,ID0,R1,fl, inf0))
}