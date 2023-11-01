vCCV=function(m){
  
  Cor=cor(t(m))[1,2]+1
  CV=sd(m)/mean(m)
  CVf=(1-CV/(0.3+CV))
  x=Cor+CVf    
  return(x)

}