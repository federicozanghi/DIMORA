#
# plot.Dimora <- function(fit,mode=c("i","c"),typeog=1,oos=c("y","n"),xlim=NULL,ltycurve=1,lwdcurve=1.5,typecurve="l",colcurve=2,...){
#   if(is.null(xlim)) xmax <- length(fit$data)
#   else xmax <- xlim[2]
#   maxx <- max(length(fit$data),xmax)
#   if(oos=="y") maxx <- maxx
#   else maxx <- length(fit$data)
#   if(mode=="i"){
#     dati <- fit$data
#     stime <- predict(fit,1:maxx)
#   }
#   else{
#     dati <- cumsum(fit$data)
#     stime <- cumsum(predict(fit,1:maxx))
#   }
#   if(typeog == 0){
#     plot(dati,xlim = xlim,main="Fitted values",ylab="Fitted values",xlab="t",...)
#     lines(stime,lty=ltycurve,type = typecurve,lwd=lwdcurve,col=colcurve)
#
#     invisible(readline(prompt="Press <enter> to continue: "))
#
#     plot(residuals(fit)~fitted(fit),main="Residuals vs Fitted",ylab="Residuals",xlab="Fitted values",...)
#     abline(h=0,lty=2,col="gray")
#
#     invisible(readline(prompt="Press <enter> to continue: "))
#
#     plot(residuals(fit)~I(1:length(residuals(fit))),main="Residuals vs Time",ylab="Residuals",xlab="t",...)
#     abline(h=0,lty=2,col="gray")
#   }
#   else if(typeog==1){
#     plot(dati,xlim = xlim,main="Fitted values",ylab="Fitted values",xlab="t",...)
#     lines(stime,lty=ltycurve,type = typecurve,lwd=lwdcurve,col=colcurve)
#   }
#   else if(typeog==2){
#     plot(residuals(fit)~fitted(fit),main="Residuals vs Fitted",ylab="Residuals",xlab="Fitted values",...)
#     abline(h=0,lty=2,col="gray")
#   }
#   else if(typeog==3){
#     plot(residuals(fit)~I(1:length(residuals(fit))),main="Residuals vs Time",ylab="Residuals",xlab="t",...)
#     abline(h=0,lty=2,col="gray")
#   }
#
# }




plot.Dimora <- function(x,..., product=1, mode=c("i","c"),typeog=1,oos=c("y","n"), xlim=NULL,ltycurve=1,lwdcurve=1.5,typecurve="l",colcurve=2){

  if(product==1){
    if(is.null(xlim)) xmax <- length(x$data)
    else xmax <- xlim[2]
    maxx <- max(length(x$data),xmax)
    if(oos=="y") maxx <- maxx
    else maxx <- length(x$data)
    if(mode=="i"){
      dati <- x$data
      stime <- predict(x,newx=1:maxx)
      laby <- "Instantaneous"
    }
    else{
      dati <- cumsum(x$data)
      stime <- cumsum(predict(x,newx=1:maxx))
      laby <- "Cumulative"
    }
    if(typeog == 0){
      plot(dati,xlim = xlim,main="Fitted values",ylab=laby,xlab="t",...)
      lines(stime,lty=ltycurve,type = typecurve,lwd=lwdcurve,col=colcurve)

      invisible(readline(prompt="Press <enter> to continue: "))

      plot(residuals(x)~fitted(x),main="Residuals vs Fitted",ylab="Residuals",xlab="Fitted values",...)
      abline(h=0,lty=2,col="gray")

      invisible(readline(prompt="Press <enter> to continue: "))

      plot(residuals(x)~I(1:length(residuals(x))),main="Residuals vs Time",ylab="Residuals",xlab="t",...)
      abline(h=0,lty=2,col="gray")
    }
    else if(typeog==1){
      plot(dati,xlim = xlim,main="Fitted values",ylab=laby,xlab="t",...)
      lines(stime,lty=ltycurve,type = typecurve,lwd=lwdcurve,col=colcurve)
    }
    else if(typeog==2){
      plot(residuals(x)~fitted(x),main="Residuals vs Fitted",ylab="Residuals",xlab="Fitted values",...)
      abline(h=0,lty=2,col="gray")
    }
    else if(typeog==3){
      plot(residuals(x)~I(1:length(residuals(x))),main="Residuals vs Time",ylab="Residuals",xlab="t",...)
      abline(h=0,lty=2,col="gray")
    }
  }

  if(product==2){
    if(mode=="i"){
      dati1 <- x$data[[1]]
      dati2 <- x$data[[2]]
      stime1 <- x$fitted[[1]]
      stime2 <- x$fitted[[2]]
      laby <- "Instantaneous"
    }
    else{
      dati1 <- cumsum(x$data[[1]])
      stime1 <- cumsum(x$fitted[[1]])
      dati2 <- cumsum(x$data[[2]])
      stime2 <- cumsum(x$fitted[[2]])
      laby <- "Cumulative"
    }
    if(typeog == 0){
      par(mfrow=c(1,2))
      plot(dati1,main="Fitted values Product 1",ylab=laby,xlab="t",...)
      lines(stime1,lty=ltycurve,type = typecurve,lwd=lwdcurve,col=colcurve)
      plot(dati2,main="Fitted values Product 2",ylab=laby,xlab="t",...)
      lines(stime2,lty=ltycurve,type = typecurve,lwd=lwdcurve,col=colcurve)

      invisible(readline(prompt="Press <enter> to continue: "))

      plot(x$residuals[[1]]~x$fitted[[1]],main="Residuals 1 vs Fitted 1",ylab="Residuals",xlab="Fitted values",...)
      abline(h=0,lty=2,col="gray")
      plot(x$residuals[[2]]~x$fitted[[2]],main="Residuals 2 vs Fitted 2",ylab="Residuals",xlab="Fitted values",...)
      abline(h=0,lty=2,col="gray")

      invisible(readline(prompt="Press <enter> to continue: "))

      plot(x$residuals[[1]]~I(1:length(x$residuals[[1]])),main="Residuals 1 vs Time",ylab="Residuals",xlab="t",...)
      abline(h=0,lty=2,col="gray")
      plot(x$residuals[[2]]~I(1:length(x$residuals[[2]])),main="Residuals 2 vs Time",ylab="Residuals",xlab="t",...)
      abline(h=0,lty=2,col="gray")
      par(mfrow=c(1,1))
    }
    else if(typeog==1){
      par(mfrow=c(1,2))
      plot(dati1,main="Fitted values Product 1",ylab=laby,xlab="t",...)
      lines(stime1,lty=ltycurve,type = typecurve,lwd=lwdcurve,col=colcurve)
      plot(dati2,main="Fitted values Product 2",ylab=laby,xlab="t",...)
      lines(stime2,lty=ltycurve,type = typecurve,lwd=lwdcurve,col=colcurve)
      par(mfrow=c(1,1))

    }
    else if(typeog==2){
      par(mfrow=c(1,2))
      plot(x$residuals[[1]]~x$fitted[[1]],main="Residuals 1 vs Fitted 1",ylab="Residuals",xlab="Fitted values",...)
      abline(h=0,lty=2,col="gray")
      plot(x$residuals[[2]]~x$fitted[[2]],main="Residuals 2 vs Fitted 2",ylab="Residuals",xlab="Fitted values",...)
      abline(h=0,lty=2,col="gray")
      par(mfrow=c(1,1))
    }
    else if(typeog==3){
      par(mfrow=c(1,2))
      plot(x$residuals[[1]]~I(1:length(x$residuals[[1]])),main="Residuals 1 vs Time",ylab="Residuals",xlab="t",...)
      abline(h=0,lty=2,col="gray")
      plot(x$residuals[[2]]~I(1:length(x$residuals[[2]])),main="Residuals 2 vs Time",ylab="Residuals",xlab="t",...)
      abline(h=0,lty=2,col="gray")
      par(mfrow=c(1,1))
    }
  }

}
