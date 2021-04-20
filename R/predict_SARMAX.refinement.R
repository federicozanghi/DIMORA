predict_SARMAX.refinement <- function(fit,newx){
  xreg <- fit$xreg
  prediction <- predict(fit,newxreg=newx)
  return(prediction)
}
