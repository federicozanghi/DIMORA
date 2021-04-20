
GG.model<-
  function (sales, prelimestimates = NULL, mt = "base", alpha = 0.05, ous=100, display=T, max.iter=100,
            ...)
  {
    x <- NULL
    t <- seq(1, length(sales), by = 1)
    s <- sales
    c <- cumsum(sales)
    if (is.function(mt) == T) {
      ff <- function(t, k, ps, qs) {
        k * mt(t) * (1 - exp(-(ps + qs) * t))/(1 + (qs/ps) *
                                                 exp(-(ps + qs) * t))
      }
      ff1 <- function(t, par) {
        c - ff(t, par[1], par[2], par[3])
      }
      ff2 <- function(t, par) {
        ff(t, par[1], par[2], par[3])
      }
      if (is.null(prelimestimates) == TRUE) {
        prelimestimates <- BASS.standard(sales = s, display = F)$Estimate[,
                                                                          1]
      }
      stime <- nls.lm(par = prelimestimates, fn = ff1, t = t,
                      control = nls.lm.control(maxiter = max.iter))$par
      aa <- data.frame(summary(nls.lm(par = prelimestimates,
                                      fn = ff1, t = t, control = nls.lm.control(maxiter = max.iter)))$coefficients[,
                                                                                                              c(1, 2)], 0, 0, 0)
      names(aa) <- c("Estimate", "Std.Error", "Lower", "Upper",
                     "p-value")
      row.names(aa) <- c("k :", "ps :", "qs :")
      for (i in 1:NROW(aa)) {
        aa[i, c(3, 4)] <- aa[i, 1] + c(-1, 1) * qnorm(1 -
                                                        alpha/2) * aa[i, 2]
      }
      sssss <- signif(summary(nls.lm(par = prelimestimates,
                                     fn = ff1, t = t, control = nls.lm.control(maxiter = max.iter)))$coefficients,
                      digits = 3)
      aa[, 5] <- sssss[, 4]
    }
    else {
      mt <- function(t, k, pc, qc) {
        k * sqrt((1 - exp(-(pc + qc) * t))/(1 + (qc/pc) *
                                              exp(-(pc + qc) * t)))
      }
      ff0 <- function(t, k, pc, qc, ps, qs) {
        mt(t, k, pc, qc) * (1 - exp(-(ps + qs) * t))/(1 +
                                                        (qs/ps) * exp(-(ps + qs) * t))
      }
      ff1 <- function(t, par) {
        c - ff0(t, par[1], par[2], par[3], par[4], par[5])
      }
      ff2 <- function(t, par) {
        ff0(t, par[1], par[2], par[3], par[4], par[5])
      }
      ff2.return <- function(t,par){
        grad(func = function(t){ff2(t,par)},t,method = "simple")
      }

      if (is.null(prelimestimates) == TRUE) {
        prelimestimates <- c(BASS.standard(sales = s, display = F)$Estimate[,
                                                                            1], 0.001, 0.1)
      }
      stime <- nls.lm(par = prelimestimates, fn = ff1, t = t,
                      control = nls.lm.control(maxiter = max.iter))$par
      aa <- data.frame(summary(nls.lm(par = prelimestimates,
                                      fn = ff1, t = t, control = nls.lm.control(maxiter = max.iter)))$coefficients[,
                                                                                                              c(1, 2)], 0, 0, 0)
      names(aa) <- c("Estimate", "Std.Error", "Lower", "Upper",
                     "p-value")
      row.names(aa) <- c("k  ", "pc  ", "qc  ", "ps  ", "qs  ")
      for (i in 1:NROW(aa)) {
        aa[i, c(3, 4)] <- aa[i, 1] + c(-1, 1) * qnorm(1 -
                                                        alpha/2) * aa[i, 2]
      }
      sssss <- signif(summary(nls.lm(par = prelimestimates,
                                     fn = ff1, t = t, control = nls.lm.control(maxiter = max.iter)))$coefficients,
                      digits = 3)
      aa[, 5] <- sssss[, 4]
    }
    if(display == TRUE){
      par(mfrow = c(1, 2))
      plot(t, c, main = "", xlim = c(0,
                                               max(t) + ous), ylim = c(0, ff2(max(t) + ous, stime)), ylab = "Cumulative")
      curve(ff2(x, stime), add = T, xlim = c(0, max(t) + ous),col=2,
            ...)
      plot(t, sales, xlim = c(0, max(t) + ous), main = "",
           ylab = "Instantaneous")
      curve(grad(function(t) ff2(t, stime), x, method = "simple"),
            add = T,col=2, ...)
      par(mfrow = c(1, 1))
    }
    s.hat <- ff2(t, stime)
    tss <- sum((c - mean(c))^2)
    rss <- sum((c - s.hat)^2)
    r.squared <- 1 - rss/tss
    r.squared.adj <- 1 - ((1 - r.squared) * (length(s) - 1))/(length(s) -
                                                                1 - NROW(aa))
    res <- nls.lm(par = prelimestimates, fn = ff1, t = t, control = nls.lm.control(maxiter = max.iter))$fvec
    coefi <- aa$Estimate
    names(coefi) <- rownames(aa)
    s.hat <- make.instantaneous(s.hat)
    cl <- match.call()
    ao <- list(model=ff2.return,type="Guseo Guidolin Model",Estimate = aa, Rsquared = r.squared, RsquaredAdj = r.squared.adj,
         residuals = res, fitted=s.hat,RSS=rss,coefficients=coefi,data=sales,call=cl)
    class(ao) <- "Dimora"
    invisible(ao)
  }

