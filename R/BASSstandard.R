BASS.standard <- function(sales,method='nls',prelimestimates=c(sum(sales)+100, 0.01, 0.1),alpha=0.05,display=T){
  #alpha=0.05
  x <- NULL
  t <- seq(1,length(sales),by=1)
  s <- sales
  c <- cumsum(s)

  ff <- function(t,m,p,q) { m* (1 - exp(- (p+q)*t )) / (1 + q/p * exp(- (p+q)*t )) }
  ff1 <- function(t,par) { c - ff(t,par[1],par[2],par[3]) }
  ff2 <- function(t,par) { ff(t,par[1],par[2],par[3]) }
  zprimo <- function(t,m,p,q){ m * ( p+q*(ff(t,m,p,q)/m) )*( 1 - (ff(t,m,p,q)/m) ) }

  if(method=='nls'){
    stime <- nls.lm(par=prelimestimates, fn=ff1, t=t)$par
    sssss <- signif(summary(nls.lm(par=prelimestimates, fn=ff1,t=t))$coefficients,digits=3)
    #print(sssss)
    aa <- data.frame(summary(nls.lm(par=prelimestimates, fn=ff1,t=t))$coefficients[,c(1,2)],0,0,0)
    names(aa) <- c('Estimate','Std.Error','Lower','Upper','p-value')
    row.names(aa) <- c('m :','p :','q :')
    for(i in 1:NROW(aa)){
      aa[i,c(3,4)] <- aa[i,1]+c(-1,1)*qnorm(1-alpha/2)*aa[i,2]
    }
    aa[,5] <- sssss[,4]
    res <- nls.lm(par=prelimestimates, fn=ff1, t=t)$fvec
  }
  else if(method=='optim'){
    mass <- sum(s)+1000
    ff3 <- function(par) { ff(t,par[1],par[2],par[3]) }
    max <- sum(s)+10000
    stima_optim <- function(c){
      f <- function(p) { sum( (c - ff3(p) )^2 )}
      optim(par = prelimestimates, fn = f ,method = "L-BFGS-B",lower = c(0.0000000001,0.0000000001,0.0000000001) , upper = c(mass,1,1))$par
    }
    stime  <- stima_optim(c)
    aa <- stime
    res <- c - ff2(t,aa)
  }


  if(display==T){
    par(mfrow=c(1,2),mar=c(5, 4, 4, 2))
    plot(t,c,main='Cumulative Sales of BM',xlim = c(0,max(t)+100),ylim = c(0,sum(s)+sum(s)*50/100),ylab = 'Cumulative')
    curve(ff2(x,stime),add=T,col=2,xlim = c(0,max(t)+100))
    plot(t,sales,main='Instantaneous Sales of BM',xlim = c(0,max(t)+100),ylab = 'Instantaneous')
    curve(zprimo(x,stime[1],stime[2],stime[3]), col=2,add=T,xlim = c(0,max(t)+100))
    par(mfrow=c(1,1),mar=c(5, 4, 4, 2))
  }

  s.hat <- ff2(t,stime)
  tss <- sum( (c - mean(s))^2 )
  rss <- sum( (c - s.hat)^2 )
  r.squared <- 1-rss/tss
  r.squared.adj <- 1 - ( (1- r.squared)*(length(s)-1) )/( length(s) - 1 - NROW(aa))


  ########### restituisco i valori
  list(Estimate=aa , Rsquared=r.squared , RsquaredAdj=r.squared.adj,residuals=res)
}

