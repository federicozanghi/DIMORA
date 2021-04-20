SARMAX.refinement <- function(fit,arima_order,seasonal_order,period=1,display=F){
  #require(forecast)
  xreg <- fitted(fit)
  model <- Arima(fit$data, order = arima_order,seasonal = list(order=seasonal_order,period=period),
                 xreg = xreg)
  model$prediction <- predict(model,newxreg=1:100)
  pred_sarmax <- xreg
  model$xreg <- xreg
  res <- residuals(model)
  if(display==F){
    return(model)
  }
  else if(display==T){
    plot(residuals(model),main="Residuals")
    abline(h=0,col="grey",lty=2)
    invisible(readline(prompt="Press <enter> to continue: "))
    acf(residuals(model),main="ACF Residuals")
    invisible(readline(prompt="Press <enter> to continue: "))
    plot(fit$data,main="Fitted values",ylab="Instantaneous")
    lines(pred_sarmax,col=2)
    return(model)
  }
}
