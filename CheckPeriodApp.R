library(qpcR)
library(lmtest)
library(lawstat)
library(gplots)

CheckPeriodApp <- function(x, cycle = 20, loess.span = 0.1, peak.window = 5, plot = TRUE, ...) {
  
  SEL <- x[, 2]
  POS <- x[, 1]
  
  ## make quadratic model
  nSample <- 1L:length(SEL)
  LM <- lm(SEL ~ I(nSample) + I(nSample^2))
  FITTED <- fitted(LM)
  RESID <- residuals(LM)
  COEFS <- as.numeric(summary(LM)$coefficients[2:3, 4])
  
  LOESS <- loess(RESID ~ nSample, span = loess.span)
  RUNS <- runs.test(RESID)
  BOX <- Box.test(RESID, lag = 1, type = "Ljung-Box")
  ACF <- acf(RESID, lag = length(RESID), plot = FALSE)
  
  FP <- findpeaks(as.numeric(ACF$acf), bw = peak.window)

  PERIOD <- diff(FP$xmax)
  
  ## return statistics as list and for "textplot"ing
  lmStat <- c(signif(COEFS[1], 3), signif(COEFS[2], 3))
  names(lmStat) <- c("lin", "quad")
  randStat <- c(signif(RUNS$p.value, 3), signif(BOX$p.value, 3))
  names(randStat) <- c("Runs test", "Ljung-Box test")
  periodStat <- c(round(mean(PERIOD, na.rm = TRUE), 1), round(sd(PERIOD, na.rm = TRUE), 1))
  names(periodStat) <- c("Period (Mean)", "Period (SD)") 
  
  list(SEL = SEL,
       POS = POS,
       FITTED = FITTED,
       COEFS = COEFS,
       RESID = RESID,
       LOESS = LOESS,
       RUNS = RUNS,
       BOX = BOX,
       PERIOD = PERIOD,
       FP = FP,
       ACF = ACF,
       #output for users
       user = list(lmStat = lmStat, randStat = randStat, periodStat = periodStat))
}



## plot data with fitted line from quadratic model 
## and show coefficients
plotFit <- function(PeriodApp) {
  COL <- rainbow(length(PeriodApp[["FITTED"]]))
  plot(PeriodApp[["SEL"]], col = COL, pch = 16)
  lines(PeriodApp[["FITTED"]], col = 1, lwd = 2)
  title(main = paste("lin:", signif(PeriodApp[["COEFS"]][1], 3), "      quad:", 
                     signif(PeriodApp[["COEFS"]][2], 3)), line = -1.5, cex.main = 1.5)
}

## plot residuals with Loess fit
## and do Runs test & Ljung-Box Q-test
plotRes <- function(PeriodApp) {
  COL <- rainbow(length(PeriodApp[["FITTED"]]))
  plot(PeriodApp[["RESID"]], col = COL, pch = 16)
  lines(fitted(PeriodApp[["LOESS"]]), lwd = 2)
  title(main = paste("Runs test:", signif(PeriodApp[["RUNS"]][["p.value"]], 3), 
                     "      Ljung-Box test:", signif(PeriodApp[["BOX"]][["p.value"]], 3)),
        line = -1.5, cex.main = 1.5)
}

## plot autocorrelation 
plotAc <- function(PeriodApp) {
  BP <- barplot(PeriodApp[["ACF"]][["acf"]][, , 1], col = "darkgrey", space = 2)
  xPos <- BP[, 1]
  points(xPos[PeriodApp[["FP"]][["xmax"]]], 1.2 * PeriodApp[["FP"]][["ymax"]], pch = 25,  col = "darkred")
  title(main = paste("Estimated Periodicity:", round(mean(PeriodApp[["PERIOD"]], na.rm = TRUE), 1), "\u00B1", 
                     round(sd(PeriodApp[["PERIOD"]], na.rm = TRUE), 1)), line = -1.5, cex.main = 1.5)
}

plotHm <- function(PeriodApp) {
  pos <- as.character(PeriodApp[["POS"]])
  pos_x <- factor(substr(pos, 0, 1))
  pos_y <- factor(as.numeric(vapply(pos, function(i) 
    substr(i, 2, nchar(i)), "a", USE.NAMES = FALSE)))
  levels(pos_y) <- levels(pos_y)[order(as.numeric(levels(pos_y)))]
  df <- data.frame(x = pos_x, y = pos_y, val = PeriodApp[["SEL"]], lab = PeriodApp[["POS"]])
  ggplot(df, aes(x = x, y = y, fill = val, label = lab)) +
    geom_tile() +
    geom_text() + 
    scale_x_discrete("") + 
    scale_y_discrete("") +
    scale_fill_continuous("Cq", low = "lightblue", high = "springgreen4") + 
    theme_bw() +
    theme(plot.background=element_blank(),
          panel.border = element_blank())
}


## Taken from https://rtricks.wordpress.com/2009/05/03/an-algorithm-to-find-local-extrema-in-a-vector/
findpeaks <- function(vec, bw = 1 , x.coo = c(1:length(vec)))
{
  pos.x.max <- NULL
  pos.y.max <- NULL
  pos.x.min <- NULL
  pos.y.min <- NULL   
  for(i in 1:(length(vec) - 1)) {
    if((i + 1 + bw) > length(vec)) sup.stop <- length(vec) else sup.stop <- i + 1 + bw
    if((i - bw) < 1) inf.stop <- 1 else inf.stop <- i - bw
    
    subset.sup <- vec[(i + 1):sup.stop]
    subset.inf <- vec[inf.stop:(i - 1)]
    
    is.max   <- sum(subset.inf > vec[i]) == 0
    is.nomin <- sum(subset.sup > vec[i]) == 0
    no.max   <- sum(subset.inf > vec[i]) == length(subset.inf)
    no.nomin <- sum(subset.sup > vec[i]) == length(subset.sup)
    
    if(is.max & is.nomin){
      pos.x.max <- c(pos.x.max,x.coo[i])
      pos.y.max <- c(pos.y.max,vec[i])
    }
    if(no.max & no.nomin){
      pos.x.min <- c(pos.x.min,x.coo[i])
      pos.y.min <- c(pos.y.min,vec[i])
    }
  }
  return(list(xmax = pos.x.max, ymax = pos.y.max, xmin = pos.x.min, ymin = pos.y.min))
}



