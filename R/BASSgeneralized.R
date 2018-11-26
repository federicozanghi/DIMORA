BASS.generalized <- function(sales,shock=c('exp','rett','mixed','armonic'),nshock,prelimestimates,alpha=0.05){
  a <- b <- c <- a1 <- b1 <- c1 <- a2 <- b2 <- c2 <- a3 <- b3 <- c3 <- x <-  NULL
  t <- seq(1,length(sales),by=1)
  s <- sales
  c <- cumsum(sales)
  if(nshock==1){
      if(shock=='exp'){
        ########### definisco la funzione: integrale di x(t)
        intx <- function(t,a1,b1,c1) {
          ( t + c1*(1/b1)*(exp(b1*(t-a1)) - 1 )*(t>=a1) )
        }
        cat('################## Esponential shock ################## \n')

        ##
        xt <- function(t,a1,b1,c1){  1 + (c1* exp(b1*(t-a1)))*(t>=a1)  }
        ##

      }
      else if (shock=='rett') {
        ########### definisco la funzione: integrale di x(t)
        intx <- function(t,a1,b1,c1) {
          ( t + c1 * (t-a1) *(a1<=t)*(t<=b1) + c1*(b1-a1)*(b1<t) )
        }
        cat('################## Rectangular shock  ################## \n')

        ##
        xt <- function(t,a1,b1,c1){ 1 + c1*(t>=a1)*(t<=b1)}
        ##

      }
      else if (shock=='mixed'){
        ########### definisco la funzione: integrale di x(t)
        intx <- function(t,a1,b1,c1) {
          t + c1*(1/b1)*(exp(b1*(t-a1)) - 1 )*(t>=a1)
        }
        cat('################## Mixed shock ################## \n')

        ##
        xt <- function(t,a1,b1,c1){  1 + (c1* exp(b1*(t-a1)))*(t>=a1)  }
        ##

      }
      else if (shock=='armonic'){
        ########### definisco la funzione: integrale di x(t)

        intx <- function(t,a1,b1,c1) {
          t + c1*( (b1-a1)/(2*pi) )*cos( 2*pi*( (t-a1)/(b1-a1) ) )*(t>=a1)*(b<=b1)
        }

        cat('################## Armonic shock  ################## \n')

        ##
        xt <- function(t,a1,b1,c1){
          1 + c1*cos( 2*pi*( (t-a1)/(b1-a1) ) )*(t>=a1)*(b<=b1)
        }
        ##
      }
      ########## stimo e ritorno i risultati
      #### !!!! cambia s in c quando passo i dati istantanei
      ff <- function(t,m,p,q,a1,b1,c1) { m*(1 - exp(- (p+q)* intx(t,a1,b1,c1))) / (1 + (q/p) * exp(- (p+q) * intx(t,a1,b1,c1))) }
      ff1 <- function(t,par) {c - ff(t,par[1],par[2],par[3],par[4],par[5],par[6])}
      ff2 <- function(t,par) {ff(t,par[1],par[2],par[3],par[4],par[5],par[6])}
      zprimo <- function(t,m,p,q,a1,b1,c1) {
        m*(p+q*(ff(t,m,p,q,a1,b1,c1)/m))*(1-(ff(t,m,p,q,a1,b1,c1)/m))*xt(t,a1,b1,c1)
      }
      ########### ottengo stime e std.Error e IC

      stime <- nls.lm(par=prelimestimates, fn=ff1, t=t)$par
      aa <- data.frame(summary(nls.lm(par=prelimestimates, fn=ff1,t=t))$coefficients[,c(1,2)],0,0,0)
      names(aa) <- c('Estimate','Std.Error','Lower','Upper','P-value')
      row.names(aa) <- c('m :','p :','q :','a1 :','b1 :','c1 :')
      for(i in 1:NROW(aa)){
        aa[i,c(3,4)] <- aa[i,1]+c(-1,1)*qnorm(1-alpha/2)*aa[i,2]
      }
      sssss <- signif(summary(nls.lm(par=prelimestimates, fn=ff1,t=t))$coefficients,digits=3)
      aa[,5] <- sssss[,4]
      resid <- nls.lm(par=prelimestimates, fn=ff1, t=t)$fvec

      ########### plotto grafico dati e curva teorica

      par(mfrow=c(1,2))
      plot(t,c,xlim=c(0,length(t)+200),ylim = c(0,sum(s)+sum(s)*0.4) ,main='Cumulative Sales of GBM',ylab = 'Cumulative')
      curve(ff2(x,stime),add=T,col=2)
      plot(t,sales,main='Instantaneous Sales of GBM',xlim=c(0,length(t)+200),ylab = 'Instantaneous')
      curve(zprimo(x,stime[1],stime[2],stime[3],stime[4],stime[5],stime[6]), col=2,add=T)
      par(mfrow=c(1,1))

      ########### mi calcolo R^2

      s.hat <- ff2(t,stime)
      tss <- sum( (c - mean(s))^2 )
      rss <- sum( (c - s.hat)^2 )
      r.squared <- 1-rss/tss
      r.squared.adj <- 1 - ( (1- r.squared)*(length(s)-1) )/( length(s) - 1 - NROW(aa))

      ########### restituisco i valori


      list(Estimate=aa , Rsquared=r.squared , RsquaredAdj=r.squared.adj, RSS=rss, residuals=resid)
    }
  else if (nshock==2){
    if(shock=='exp'){
      ########### definisco la funzione: integrale di x(t)
      intx1 <- function(t,a1,b1,c1,a2,b2,c2) {
        t + c1*(1/b1)*(exp(b1*(t-a1)) - 1 )*(t>=a1) + c2*(1/b2)*(exp(b2*(t-a2)) - 1 )*(t>=a2)
      }
      cat('################## Esponential shock ################## \n')

      ##
      xt1 <- function(t,a1,b1,c1,a2,b2,c2){  1 + (c1* exp(b1*(t-a1)))*(t>=a1) + (c2* exp(b2*(t-a2)))*(t>=a2) }
      ##
    }
    else if (shock=='rett') {
      ########### definisco la funzione: integrale di x(t)
      intx1 <- function(t,a1,b1,c1,a2,b2,c2) {
        t + c1 * (t-a1) *(a1<=t)*(t<=b1) + c1*(b1-a1)*(b1<t) + c2 * (t-a2) *(a2<=t)*(t<=b2) + c2*(b2-a2)*(b2<t)
      }
      cat('################## Rectangular shock  ################## \n')

      ##
      xt1 <- function(t,a1,b1,c1,a2,b2,c2){  1 + c1*(t>=a1)*(t<=b1) + c2*(t>=a2)*(t<=b2) }
      ##

    }
    else if (shock=='mixed'){
      ########### definisco la funzione: integrale di x(t)
      intx1 <- function(t,a1,b1,c1,a2,b2,c2) {
        t + (c1/b1)*(exp(b1*(t-a1)) - 1 )*(a1<=t) + c2*(t-a2)*(a2<=t)*(t<=b2) + c2*(b2-a2)*(b2<t)
      }
      cat('################## Mixed shock ################## \n')

      ##
      xt1 <- function(t,a1,b1,c1,a2,b2,c2){  1 + (c1* exp(b1*(t-a1)))*(t>=a1) + c2*(t>=a2)*(t<=b2) }
      ##

    }
    else if (shock=='armonic'){
      ########### definisco la funzione: integrale di x(t)

      xt1 <- function(t,a1,b1,c1,a2,b2,c2){
        1 + c1*cos( 2*pi*( (t-a1)/(b1-a1) ) )*(t>=a1)*(b<=b1) + c2*cos( 2*pi*( (t-a2)/(b2-a2) ) )*(t>=a2)*(b<=b2)
      }
      intx1 <- function(t,a1,b1,c1,a2,b2,c2) {
        t + c1*( (b1-a1)/(2*pi) )*cos( 2*pi*( (t-a1)/(b1-a1) ) )*(t>=a1)*(b<=b1) +
          c2*( (b2-a2)/(2*pi) )*cos( 2*pi*( (t-a2)/(b2-a2) ) )*(t>=a2)*(b<=b2)
      }

      cat('################## Armonic shock  ################## \n')

    }
    ########## stimo e ritorno i risultati
    # ricorda cambiare s con c
    ff0 <- function(t,m,p,q,a1,b1,c1,a2,b2,c2) { m*(1 - exp(- (p+q)* intx1(t,a1,b1,c1,a2,b2,c2))) / (1 + (q/p) * exp(- (p+q) * intx1(t,a1,b1,c1,a2,b2,c2))) }
    ff1 <- function(t,par) {c - ff0(t,par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8],par[9])}
    ff2 <- function(t,par) {ff0(t,par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8],par[9])}
    #zprimo <- function(t,m,p,q,a1,b1,c1,a2,b2,c2) { m*(( p* ((p+q)^2) ) * xt(t,a1,b1,c1,a2,b2,c2) * exp(- (p+q)* intx1(t,a1,b1,c1,a2,b2,c2)))/( (p*exp(- (p+q) * intx1(t,a1,b1,c1,a2,b2,c2)) + q)^2 ) }
    zprimo1 <- function(t,m,p,q,a1,b1,c1,a2,b2,c2) {
      m*(p+q*(ff0(t,m,p,q,a1,b1,c1,a2,b2,c2)/m))*(1-(ff0(t,m,p,q,a1,b1,c1,a2,b2,c2)/m))*xt1(t,a1,b1,c1,a2,b2,c2)
    }

    ########### ottengo stime e std.Error e IC

    stime <- nls.lm(par=prelimestimates, fn=ff1, t=t)$par
    aa <- data.frame(summary(nls.lm(par=prelimestimates, fn=ff1,t=t))$coefficients[,c(1,2)],0,0,0)
    names(aa) <- c('Estimate','Std.Error','Lower','Upper','P-value')
    row.names(aa) <- c('m :','p :','q :','a1 :','b1 :','c1 :','a2 :','b2 :','c2 :')
    for(i in 1:NROW(aa)){
      aa[i,c(3,4)] <- aa[i,1]+c(-1,1)*qnorm(1-alpha/2)*aa[i,2]
    }
    sssss <- signif(summary(nls.lm(par=prelimestimates, fn=ff1,t=t))$coefficients,digits=3)
    aa[,5] <- sssss[,4]
    resid <- nls.lm(par=prelimestimates, fn=ff1, t=t)$fvec

    ########### plotto grafico dati e curva teorica
    par(mfrow=c(1,2))
    plot(t,c,xlim=c(0,length(t)+100),ylim = c(0,max(s)+5000),main='Cumulative Sales of GBM')
    curve(ff2(x,stime),type='l',add=T,col=2)
    plot(t,sales,main='Instantaneous Sales of GBM')
    curve(zprimo1(x,stime[1],stime[2],stime[3],stime[4],stime[5],stime[6],stime[7],stime[8],stime[9]),col=2,add=T)
    par(mfrow=c(1,1))
    ########### mi calcolo R^2

    s.hat <- ff2(t,stime)
    tss <- sum( (c - mean(s))^2 )
    rss <- sum( (c - s.hat)^2 )
    r.squared <- 1-rss/tss
    r.squared.adj <- 1 - ( (1- r.squared)*(length(s)-1) )/( length(s) - 1 - NROW(aa))

    ########### restituisco i valori


    list(Estimate=aa , Rsquared=r.squared , RsquaredAdj=r.squared.adj, RSS=rss, residuals=resid)
  }
  else if (nshock==3){
    if(shock=='exp' | shock=='rett' | shock=='armonic'){
      if(shock=='exp'){
        ########### definisco la funzione: integrale di x(t)
        intx2 <- function(t,a1,b1,c1,a2,b2,c2,a3,b3,c3) {
          t + c1*(1/b1)*(exp(b1*(t-a1)) - 1 )*(t>=a1) + c2*(1/b2)*(exp(b2*(t-a2)) - 1 )*(t>=a2) + c3*(1/b3)*(exp(b3*(t-a3)) - 1 )*(t>=a3)
        }
        cat('################## Esponential shock ################## \n')

        ##
        xt2 <- function(t,a1,b1,c1,a2,b2,c2,a3,b3,c3){  1 + (c1* exp(b1*(t-a1)))*(t>=a1) + (c2* exp(b2*(t-a2)))*(t>=a2) + (c3* exp(b3*(t-a3)))*(t>=a3) }
        ##

      }
      else if (shock=='rett') {
        ########### definisco la funzione: integrale di x(t)
        intx2 <- function(t,a1,b1,c1,a2,b2,c2,a3,b3,c3) {
          ( t + c1 * (t-a1) *(a1<=t)*(t<=b1) + c1*(b1-a1)*(b1<t) + c2 * (t-a2) *(a2<=t)*(t<=b2) + c2*(b2-a2)*(b2<t) + c3 * (t-a3) *(a3<=t)*(t<=b3) + c3*(b3-a3)*(b3<t) )
        }
        cat('################## Rectangular shock  ################## \n')

        ##
        xt2 <- function(t,a1,b1,c1,a2,b2,c2,a3,b3,c3){  1 + c1*(t>=a1)*(t<=b1) + c2*(t>=a2)*(t<=b2) + c3*(t>=a3)*(t<=b3) }
        ##

      }
      else if (shock=='armonic'){
        ########### definisco la funzione: integrale di x(t)

        xt2 <- function(t,a1,b1,c1,a2,b2,c2,a3,b3,c3){
          1 + c1*cos( 2*pi*( (t-a1)/(b1-a1) ) )*(t>=a1)*(b<=b1) + c2*cos( 2*pi*( (t-a2)/(b2-a2) ) )*(t>=a2)*(b<=b2) + c3*cos( 2*pi*( (t-a3)/(b3-a3) ) )*(t>=a3)*(b<=b3)
        }
        intx2 <- function(t,a1,b1,c1,a2,b2,c2,a3,b3,c3) {
          t + c1*( (b1-a1)/(2*pi) )*cos( 2*pi*( (t-a1)/(b1-a1) ) )*(t>=a1)*(b<=b1) +
            c2*( (b2-a2)/(2*pi) )*cos( 2*pi*( (t-a2)/(b2-a2) ) )*(t>=a2)*(b<=b2) +
            c3*( (b3-a3)/(2*pi) )*cos( 2*pi*( (t-a3)/(b3-a3) ) )*(t>=a3)*(b<=b3)
        }

        cat('################## Armonic shock  ################## \n')

      }

      ########## stimo e ritorno i risultati

      ff00 <- function(t,m,p,q,a1,b1,c1,a2,b2,c2,a3,b3,c3) { m*(1 - exp(- (p+q)* intx2(t,a1,b1,c1,a2,b2,c2,a3,b3,c3))) / (1 + (q/p) * exp(- (p+q) * intx2(t,a1,b1,c1,a2,b2,c2,a3,b3,c3))) }
      ff1 <- function(t,par) {c - ff00(t,par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8],par[9],par[10],par[11],par[12])}
      ff2 <- function(t,par) {ff00(t,par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8],par[9],par[10],par[11],par[12])}
      zprimo2 <- function(t,m,p,q,a1,b1,c1,a2,b2,c2,a3,b3,c3) {
        m*(p+q*(ff00(t,m,p,q,a1,b1,c1,a2,b2,c2,a3,b3,c3)/m))*(1-(ff00(t,m,p,q,a1,b1,c1,a2,b2,c2,a3,b3,c3)/m))*xt2(t,a1,b1,c1,a2,b2,c2,a3,b3,c3)
      }

      ########### ottengo stime e std.Error e IC


      stime <- nls.lm(par=prelimestimates, fn=ff1, t=t)$par
      aa <- data.frame(summary(nls.lm(par=prelimestimates, fn=ff1,t=t))$coefficients[,c(1,2)],0,0,0)
      names(aa) <- c('Estimate','Std.Error','Lower','Upper','P-value')
      row.names(aa) <- c('m :','p :','q :','a1 :','b1 :','c1 :','a2 :','b2 :','c2 :','a3 :','b3 :','c3 :')
      for(i in 1:NROW(aa)){
        aa[i,c(3,4)] <- aa[i,1]+c(-1,1)*qnorm(1-alpha/2)*aa[i,2]
      }
      sssss <- signif(summary(nls.lm(par=prelimestimates, fn=ff1,t=t))$coefficients,digits=3)
      aa[,5] <- sssss[,4]
      resid <- nls.lm(par=prelimestimates, fn=ff1, t=t)$fvec

      ########### plotto grafico dati e curva teorica

      par(mfrow=c(1,2))
      plot(t,c,xlim=c(0,length(t)+100),ylim = c(0,max(s)+5000) ,main='Cumulative Sales of GBM')
      curve(ff2(x,stime),type='l',add=T,col=2)
      plot(t,sales,main='Instantaneous Sales of GBM')
      curve(zprimo2(x,stime[1],stime[2],stime[3],stime[4],stime[5],stime[6],stime[7],stime[8],stime[9],stime[10],stime[11],stime[12]),col=2,add=T)
      par(mfrow=c(1,1))

      ########### mi calcolo R^2

      s.hat <- ff2(t,stime)
      tss <- sum( (c - mean(s))^2 )
      rss <- sum( (c - s.hat)^2 )
      r.squared <- (tss-rss)/tss
      r.squared.adj <- 1 - ( (1- r.squared)*(length(s)-1) )/( length(s) - 1 - NROW(aa))

      ########### restituisco i valori


      list(Estimate=aa , Rsquared=r.squared , RsquaredAdj=r.squared.adj, RSS=rss, residuals=resid )
    }
    else if (shock=='mixed'){
      cat("--- Sorry you don't have 3 shocks in this case ---")
    }
  }
  else{
    cat("I'm sorry but we don't have implemented yet a model with more than 3 shocks.")
  }
}
