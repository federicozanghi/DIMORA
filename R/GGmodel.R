GG.model <- function(sales,prelimestimates=NULL,mt='base',alpha=0.05,...){
  x <- NULL
  t <- seq(1,length(sales),by=1)
  s <- sales
  c <- cumsum(sales)


  if(is.function(mt)==T){
    ff <- function(t,k,ps,qs){
      k*mt(t) * (1-exp(-(ps+qs)*t)) / (1+(qs/ps)*exp(-(ps+qs)*t))
    }
    ff1 <- function(t,par) {c - ff(t,par[1],par[2],par[3])}
    ff2 <- function(t,par) {ff(t,par[1],par[2],par[3])}
    #cat(is.null(prelimestimates)==TRUE)
    if(is.null(prelimestimates)==TRUE){
      prelimestimates <- BASS.standard(sales = s,display = F)$Estimate[,1]
    }
    ########### ottengo stime e std.Error e IC

    stime <- nls.lm(par=prelimestimates, fn=ff1, t=t,control = nls.lm.control(maxiter = 100))$par
    aa <- data.frame(summary(nls.lm(par=prelimestimates, fn=ff1,t=t,control = nls.lm.control(maxiter = 100)))$coefficients[,c(1,2)],0,0,0)
    names(aa) <- c('Estimate','Std.Error','Lower','Upper','P-value')
    row.names(aa) <- c('k :','ps :','qs :')
    for(i in 1:NROW(aa)){
      aa[i,c(3,4)] <- aa[i,1]+c(-1,1)*qnorm(1-alpha/2)*aa[i,2]
    }
    sssss <- signif(summary(nls.lm(par=prelimestimates, fn=ff1,t=t,control = nls.lm.control(maxiter = 100)))$coefficients,digits=3)
    aa[,5] <- sssss[,4]
  }
  else{
    ### mercato potenziale variabile

    mt <- function(t,k,pc,qc){
      k*sqrt( (1-exp(-(pc+qc)*t) )/( 1+(qc/pc)*exp( -(pc+qc)*t )) )
    }

    ### modello GUSEO GUIDOLIN cumulato

    ff0 <- function(t,k,pc,qc,ps,qs){
      mt(t,k,pc,qc) * (1-exp(-(ps+qs)*t)) / (1+(qs/ps)*exp(-(ps+qs)*t))
    }

    ########## stimo e ritorno i risultati

    ff1 <- function(t,par) {c - ff0(t,par[1],par[2],par[3],par[4],par[5])}
    ff2 <- function(t,par) {ff0(t,par[1],par[2],par[3],par[4],par[5])}

    if(is.null(prelimestimates)==TRUE){
      prelimestimates <- c(BASS.standard(sales = s,display=F)$Estimate[,1],1.00000e-03,1.00000e-01)
    }

    ########### ottengo stime e std.Error e IC

    stime <- nls.lm(par=prelimestimates, fn=ff1, t=t,control = nls.lm.control(maxiter = 100))$par
    aa <- data.frame(summary(nls.lm(par=prelimestimates, fn=ff1,t=t,control = nls.lm.control(maxiter = 100)))$coefficients[,c(1,2)],0,0,0)
    names(aa) <- c('Estimate','Std.Error','Lower','Upper','P-value')
    row.names(aa) <- c('k :','pc :','qc :','ps :','qs :')
    for(i in 1:NROW(aa)){
      aa[i,c(3,4)] <- aa[i,1]+c(-1,1)*qnorm(1-alpha/2)*aa[i,2]
    }
    sssss <- signif(summary(nls.lm(par=prelimestimates, fn=ff1,t=t,control = nls.lm.control(maxiter = 100)))$coefficients,digits=3)
    aa[,5] <- sssss[,4]

  }
  ########### plotto grafico dati e curva teorica

  par(mfrow=c(1,2))

  plot(t,c,main='Cumulative Sales of GGM',xlim=c(0,max(t)+200),ylim=c(0,ff2(max(t)+200,stime)+100),ylab = 'Cumulative')
  curve(ff2(x,stime),add=T,xlim=c(0,max(t)+200),...)

  plot(t,sales,xlim=c(0,max(t)+200),main='Instantaneous Sales of GGM',ylab='Instantaneous')#,ylim=c(0,50000),)
  curve(grad(function(t) ff2(t,stime),x,method = 'simple'),add=T,...)

  par(mfrow=c(1,1))

  ########### mi calcolo R^2

  s.hat <- ff2(t,stime)
  tss <- sum( (c - mean(s))^2 )
  rss <- sum( (c - s.hat)^2 )
  r.squared <- 1-rss/tss
  r.squared.adj <- 1 - ( (1- r.squared)*(length(s)-1) )/( length(s) - 1 - NROW(aa))

  ########### restituisco i valori

  res <- nls.lm(par=prelimestimates, fn=ff1, t=t,control = nls.lm.control(maxiter = 100))$fvec
  list(Estimate=aa , Rsquared=r.squared , RsquaredAdj=r.squared.adj, residuals=res )
}


