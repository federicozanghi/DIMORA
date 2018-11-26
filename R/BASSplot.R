BASS.plot<-function(data){
  par(mfrow=c(1,2))
	plot(data[,1],data[,2],type="o",lwd=2,col="dark blue",xlab="[t]",ylab="Sales")
	plot(data[,1],data[,3],type="o",lwd=2,col="dark green",xlab="[t]",ylab="Cumulative sales")
	par(mfrow=c(1,1))
}
