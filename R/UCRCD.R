
UCRCD <- function(sales1, sales2, c2, display=T, alpha=0.05,
                  delta=0.01, gamma=0.01, par="double",
                  m1 =BASS.standard(sales1,display = F)$Estimate[1,1],
                  m2 =BASS.standard(sales2,display = F)$Estimate[1,1],
                  p1c=BASS.standard(sales1,display = F)$Estimate[2,1],
                  q1c=BASS.standard(sales1,display = F)$Estimate[3,1],
                  p2c=BASS.standard(sales2,display = F)$Estimate[2,1],
                  q2c=BASS.standard(sales2,display = F)$Estimate[3,1])
{
  tot<-c(sales1,sales2)
  data1<-cumsum(sales1)
  data2<-cumsum(sales2)
  end <-length(sales1)

  if (c2 > 0){
    s1 <- sales1[1:c2]
    sales1 <- sales1[(c2+1):length(sales1)]
    t <- c(1:c2)
    s2 <- c(rep(0,c2))
    Z1 <- cumsum(s1)
    Z2 <- cumsum(s2)

    BMs1<-BASS.standard(s1,display = F)
    model1<-BMs1$model
    m1 <-BMs1$Estimate[1,1]
    p1a<-BMs1$Estimate[2,1]
    q1a<-BMs1$Estimate[3,1]

    m1_ <-BMs1$Estimate[1,]
    p1a_<-BMs1$Estimate[2,]
    q1a_<-BMs1$Estimate[3,]

    pred_BM1<-fitted(BMs1)

    o_bass<-matrix()
    o_bass <- cbind(t, s1, s2, Z1, Z2)
    p_bass <- matrix()
    p_bass <- cbind(t, pred_BM1, s2)
    o_bass<-as.data.frame(o_bass)
    p_bass<-as.data.frame(p_bass)
    colnames(o_bass)<-c("t","s1","s2", "Z1", "Z2")
    colnames(p_bass)<-c("t","pred_1","pred_2")
  }



  #######################################################
  if(par=="unique"){
    parms<-list(mc=(m1+m2)*2,
                p1c=p1c,
                p2c=p2c,
                q1c=q1c,
                q2c=q2c,
                delta=delta
    )

    t <- seq(c2+1, end, by = 1)

    Z1 <- cumsum(sales1)
    Z2 <- cumsum(sales2)
    data<-matrix()
    data<-cbind(t, sales1, sales2, Z1, Z2)
    data<-as.data.frame(data)
    colnames(data)<-c("t","s1","s2", "Z1", "Z2")
    expdf=melt(data[1:3],id.var="t",variable.name="product",value.name="response")

    model<-function(t,parms){
      Z=Z1+Z2
      mc=parms$mc
      p1c=parms$p1c
      p2c=parms$p2c
      q1c=parms$q1c
      q2c=parms$q2c
      delta=parms$delta
      z1 = mc*( p1c+(q1c+delta)*Z1/mc + q1c*Z2/mc ) * (1-Z/mc)
      z2 = mc*( p2c+(q2c-delta)*Z1/mc + q2c*Z2/mc ) * (1-Z/mc)
      return(list(z1=z1,z2=z2,t=t))
    }

    res.model<-function(parms) {
      stime=model(t=t,parms)
      data<-data.frame(t=stime$t, z1.s=stime$z1, z2.s=stime$z2)
      preddf=melt(data,id.var="t",variable.name="product",value.name="response")
      residui<-preddf$response-expdf$response
      return (residui)
    }

    fitval1=nls.lm(par=parms,fn=res.model, control=nls.lm.control(maxiter=1000,maxfev=10000))
    summary <- summary(fitval1)
    sssss<-signif(summary$coefficients, digits = 3)
    aa <- data.frame(summary$coefficients[, c(1, 2)], 0, 0,0)
    for (i in 1:NROW(aa)) {
      aa[i, c(3, 4)] <- aa[i, 1] + c(-1, 1) * qnorm(1 -
                                                      alpha/2) * aa[i, 2]
    }
    aa[, 5] <- sssss[, 4]
    names(aa) <- c("Estimate", "Std.Error", "Lower", "Upper",
                   "p-value")

    parest=as.list(coef(fitval1))
    mc=parest$mc
    p1c=parest$p1c
    q1c=parest$q1c
    p2c=parest$p2c
    q2c=parest$q2c
    delta=parest$delta


    Estimate<-c(m1,p1a,q1a,mc,p1c,q1c+delta,q1c,p2c,q2c,q2c-delta)
    names(Estimate)<-c("m1","p1a","q1a","mc","p1c","q1c+delta","q1c","p2c",
                       "q2c","q2c-delta")

    Estimate1<-rbind(m1_, p1a_, q1a_, aa)
    stime=model(t=t,parest)
    z_prime<-data.frame(t=stime$t, pred_1=stime$z1, pred_2=stime$z2)

    if(c2 > 0){
      data$t <- c((c2+1):end)
      data_o<-matrix()
      data_o <- rbind(o_bass, data)
      data <- data_o
      z_prime$t <- c((c2+1):end)
      data_p <- matrix()
      data_p <- rbind(p_bass, z_prime)
      z_prime <- data_p
    }

    oss=melt(data[,c(1:3)],id.var=c("t"),variable.name="product",value.name="consumption")
    pred=melt(z_prime,id.var=c("t"),variable.name="product",value.name="consumption")
    res=oss$consumption-pred$consumption
  }

  ##################################################################################

  if(par=="double"){
    parms<-list(mc=(m1+m2)*2,
                p1c=p1c,
                p2c=p2c,
                q1c=q1c,
                q2c=q2c,
                delta=delta,
                gamma=gamma)

    t <- seq(c2+1, end, by = 1)

    Z1 <- cumsum(sales1)
    Z2 <- cumsum(sales2)
    data<-matrix()
    data<-cbind(t, sales1, sales2, Z1, Z2)
    data<-as.data.frame(data)
    colnames(data)<-c("t","s1","s2", "Z1", "Z2")
    expdf=melt(data[1:3],id.var="t",variable.name="product",value.name="response")

    model<-function(t,parms){
      Z=Z1+Z2
      mc=parms$mc
      p1c=parms$p1c
      p2c=parms$p2c
      q1c=parms$q1c
      q2c=parms$q2c
      delta=parms$delta
      gamma=parms$gamma
      z1 = mc*( p1c+(q1c+delta)*Z1/mc + q1c*Z2/mc ) * (1-Z/mc)
      z2 = mc*( p2c+(q2c-gamma)*Z1/mc + q2c*Z2/mc ) * (1-Z/mc)
      return(list(z1=z1,z2=z2,t=t))
    }

    res.model<-function(parms) {
      stime=model(t=t,parms)
      data<-data.frame(t=stime$t, z1.s=stime$z1, z2.s=stime$z2)
      preddf=melt(data,id.var="t",variable.name="product",value.name="response")
      residui<-preddf$response-expdf$response
      return (residui)
    }

    fitval1=nls.lm(par=parms,fn=res.model, control=nls.lm.control(maxiter=1000,maxfev=10000))
    summary <- summary(fitval1)
    sssss<-signif(summary$coefficients, digits = 3)
    aa <- data.frame(summary$coefficients[, c(1, 2)], 0, 0,0)
    for (i in 1:NROW(aa)) {
      aa[i, c(3, 4)] <- aa[i, 1] + c(-1, 1) * qnorm(1 -
                                                      alpha/2) * aa[i, 2]
    }
    aa[, 5] <- sssss[, 4]
    names(aa) <- c("Estimate", "Std.Error", "Lower", "Upper",
                   "p-value")

    parest=as.list(coef(fitval1))
    mc=parest$mc
    p1c=parest$p1c
    q1c=parest$q1c
    p2c=parest$p2c
    q2c=parest$q2c
    delta=parest$delta
    gamma=parest$gamma

    Estimate<-c(m1,p1a,q1a,mc,p1c,q1c+delta,q1c,p2c,q2c,q2c-gamma)
    names(Estimate)<-c("m1","p1a","q1a","mc","p1c","q1c+delta","q1c","p2c",
                       "q2c","q2c-gamma")

    Estimate1<-rbind(m1_, p1a_, q1a_, aa)
    stime=model(t=t,parest)
    z_prime<-data.frame(t=stime$t, pred_1=stime$z1, pred_2=stime$z2)

    if(c2 > 0){
      data$t <- c((c2+1):end)
      data_o<-matrix()
      data_o <- rbind(o_bass, data)
      data <- data_o
      z_prime$t <- c((c2+1):end)
      data_p <- matrix()
      data_p <- rbind(p_bass, z_prime)
      z_prime <- data_p
    }

    oss=melt(data[,c(1:3)],id.var=c("t"),variable.name="product",value.name="consumption")
    pred=melt(z_prime,id.var=c("t"),variable.name="product",value.name="consumption")
    res=oss$consumption-pred$consumption
  }

  ################################################

  data<-list(oss$consumption[1:end],oss$consumption[(end+1+c2):(2*end)] )
  fitted<-list(pred$consumption[1:end],pred$consumption[(end+1+c2):(2*end)] )
  residu<-list(oss$consumption[1:end]-pred$consumption[1:end],
               oss$consumption[(end+1+c2):(2*end)]-pred$consumption[(end+1+c2):(2*end)] )

  Rrres <- c(oss$consumption[1:end]-pred$consumption[1:end],
             oss$consumption[(end+1+c2):(2*end)]-pred$consumption[(end+1+c2):(2*end)] )

  Description <- c("Market Potential 1 stand-alone","Innovation coefficient 1 stand-alone",
                   "Imitation coefficient 1 stand-alone","Market Potential 1 competition",
                   "Innovation coefficient 1 competition","Within-product wom 1 competition",
                   "Cross-product wom 1 competition","Innovation coefficient 2 competition",
                   "Within-product wom 2 competition", "Cross-product wom 2 competition")
  parameters<-matrix()
  parameters<-cbind(Estimate, Description)
  parameters<-as.data.frame(parameters)

  cum<-c(data1,data2)
  pred_c1<-cumsum(pred$consumption[1:end])
  pred_c2<-cumsum(pred$consumption[(end+1+c2):(2*end)])
  cum_pred<-c(pred_c1,pred_c2)
  tss <- sum((cum - mean(cum))^2)
  rss <- sum((cum - cum_pred)^2)
  r.squared <- 1 - rss/tss
  r.squared.adj <- 1 - ((1 - r.squared) * ((2*end-c2) - 1))/((2*end-c2)-1-NROW(aa))

  one<-end+1
  two<-one+c2-1

  ss1 <- oss$consumption[1:end]
  ss2 <- oss$consumption[(end+1+c2):(2*end)]
  cc1 <- cumsum(ss1)
  cc2 <- cumsum(ss2)

  pp1 <- pred$consumption[1:end]
  pp2 <- pred$consumption[(end+1+c2):(2*end)]
  gg1 <- cumsum(pp1)
  gg2 <- cumsum(pp2)

  t <- c(1:end)
  t2 <- c((c2+1):end)


  if(display==T){
    plot(t,ss1,main = "", pch = 1, ylab = "Instantaneous")
    lines(t2, pp2,col = 2)
    lines(t, pp1, col = 2)
    lines(t2, ss2, pch = 1, type = "p", col=3)
  }

  cl <- match.call()
  ao <- list(model=model, type="UCRCD Model",Estimate = Estimate1, Rsquared = r.squared,
             RsquaredAdj = r.squared.adj,
             residuals = residu, fitted=fitted,RSS=rss, Res_tot=Rrres,
             coefficients=parameters,data=data, call=cl)
  class(ao) <- "Dimora"
  invisible(ao)
}

