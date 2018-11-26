BASS.standard.generator<-function(m,p,q,tstart=1,n=50){
  data=data.frame(tstart,0,0)
  names(data)=c("time","sales","cumulative sales")
  i=1
  a=q/p
  b=p+q
  x=exp(-b*i)
  y=1+(a*x)
  adoption=m*((b^2)*x)/(p*(y^2))
  for( i in 1:n){
    data[i+1,1]=tstart+i
    data[i+1,2]=adoption
    data[i+1,3]=data[i,3]+data[i+1,2]
    i=i+1
    a=q/p
    b=p+q
    x=exp(-b*i)
    y=1+(a*x)
    adoption=m*((b^2)*x)/(p*(y^2))
  }
  return(data)
}


