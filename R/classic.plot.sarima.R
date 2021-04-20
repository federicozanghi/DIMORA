classic.plot.sarima <- function(model,typeog=0,...){
  if(typeog==0){
    plot(residuals(model),type="b",main="Residuals vs Time",ylab="residuals")
    abline(h=0,col="grey",lty=2)
    invisible(readline(prompt="Press <enter> to continue: "))
    acf(residuals(model),main="ACF Residuals")
    invisible(readline(prompt="Press <enter> to continue: "))
    plot(model$x,main="Fitted values",ylab="Instantaneous")
    lines(fitted(model),col=2)
  }
  else if(typeog==1){
    plot(residuals(model),main="Residuals vs Time",ylab="residuals",...)
    abline(h=0,col="grey",lty=2)
  }
  else if(typeog==2){
    acf(residuals(model),main="ACF Residuals")
  }
  else if(typeog==3){
    plot(model$x,main="Fitted values",ylab="Instantaneous",...)
    lines(fitted(model),col=2)
  }

}
