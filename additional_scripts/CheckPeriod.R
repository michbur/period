#library(qpcR)
#library(lmtest)
#library(lawstat)
#library(gplots)

CheckPeriod <- function(x, cycle = 20, loess.span = 0.1, peak.window = 5, plot = TRUE, ...) {
  
  ## grab data, if 'modlist', else use vector
  if (class(x)[1] == "modlist") {
    DATA <- sapply(x, function(a) a$DATA[, 2])
    SEL <- DATA[cycle, ]
  } else SEL <- x
  
  ## make quadratic model
  nSample <- 1L:length(SEL)
  LM <- lm(SEL ~ I(nSample) + I(nSample^2))
  FITTED <- fitted(LM)
  RESID <- residuals(LM)
  COEFS <- as.numeric(summary(LM)$coefficients[2:3, 4])
 
  ## plot setup
  if (plot) {
    par(mfrow = c(3, 1))
    par(mar = c(1, 3, 3, 1))
  }
  
  ## plot data with fitted line from quadratic model 
  ## and show coefficients
  COL <- rainbow(length(FITTED))
  if (plot) {
    plot(SEL, col = COL, pch = 16)
    lines(FITTED, col = 1, lwd = 2)
    title(main = paste("lin:", signif(COEFS[1], 3), "      quad:", signif(COEFS[2], 3)), line = -1.5, cex.main = 1.5)
  }
    
  ## plot residuals with Loess fit
  ## and do Runs test & Ljung-Box Q-test
  if (plot) plot(RESID, col = COL, pch = 16)
  LOESS <- loess(RESID ~ nSample, span = loess.span)
  if (plot) lines(fitted(LOESS), lwd = 2)
  RUNS <- runs.test(RESID)
  BOX <- Box.test(RESID, lag = 1, type = "Ljung-Box")
  if (plot) title(main = paste("Runs test:", signif(RUNS$p.value, 3), "      Ljung-Box test:", signif(BOX$p.value, 3)),
        line = -1.5, cex.main = 1.5)
  
  ## plot autocorrelation 
  ACF <- acf(RESID, lag = length(RESID), plot = FALSE)
  if (plot) {
    BP <- barplot(as.numeric(ACF$acf), col = "darkgrey", space = 2)
    xPos <- BP[, 1]
  }
  FP <- findpeaks(as.numeric(ACF$acf), bw = peak.window)
  if (plot) points(xPos[FP$xmax], 1.2 * FP$ymax, pch = 25,  col = "darkred")
  PERIOD <- diff(FP$xmax)
  if (plot) title(main = paste("Estimated Periodicity:", round(mean(PERIOD, na.rm = TRUE), 1), "\u00B1", 
                     round(sd(PERIOD, na.rm = TRUE), 1)), line = -1.5, cex.main = 1.5)

  ## return statistics as list and for "textplot"ing
  lmStat <- c(signif(COEFS[1], 3), signif(COEFS[2], 3))
  names(lmStat) <- c("lin", "quad")
  randStat <- c(signif(RUNS$p.value, 3), signif(BOX$p.value, 3))
  names(randStat) <- c("Runs test", "Ljung-Box test")
  periodStat <- c(round(mean(PERIOD, na.rm = TRUE), 1), round(sd(PERIOD, na.rm = TRUE), 1))
  names(periodStat) <- c("Period (Mean)", "Period (SD)") 
  
  return(list(lmStat = lmStat, randStat = randStat, periodStat = periodStat, resid = RESID))
}