library(qpcR)
library(lmtest)
library(lawstat)
library(chipPCR)
library(lattice)

################# Figure 1: The discovery of periodicity ####################
## Raw data
ml1 <- modlist(reps384, model = spl3)
COL <- rainbow(length(ml1))
plot(ml1, par2D = list(xlab = "", ylab = "", pch = 16, cex = 0.5), col = COL)

VEC <- rep(1, length(ml1))
names(VEC) <- as.character(1:379)
barplot(VEC, width = 0.2, space = 0, , border = NA, col = COL, las = 1)

## Baselined by linear model
ml2 <- modlist(reps384, model = spl3, basecyc = 1:10, baseline = "lin") 
COL <- rainbow(length(ml2))
plot(ml2, par2D = list(xlab = "", ylab = "", pch = 16, cex = 0.5), col = COL)

## plot of fluorescence values each 5 cycles
DATA <- sapply(ml2, function(x) x$DATA[, 2])
plot(1:379, DATA[10, ], col = COL, pch = 16, xlab = "", ylab = "", las = 1)
plot(1:379, DATA[20, ], col = COL, pch = 16, xlab = "", ylab = "", las = 1)
plot(1:379, DATA[30, ], col = COL, pch = 16, xlab = "", ylab = "", las = 1)
plot(1:379, DATA[40, ], col = COL, pch = 16, xlab = "", ylab = "", las = 1)

######################### Figure 2: The analysis pipeline #################
## fitted curve of quadratic model
LM <- lm(DATA[20, ] ~ I(1:379) + I((1:379)^2))
plot(1:379, DATA[20, ], col = COL, pch = 16, xlab = "", ylab = "", las = 1)
FITTED <- fitted(LM)
lines(1:379, FITTED, col = 1, lwd = 2)

## plot residuals and add loess
RESID <- residuals(LM)
plot(1:379, RESID, col = COL, pch = 16, xlab = "", ylab = "", las = 1)
LOESS <- loess(RESID ~ I(1:379), span = 0.1)
lines(fitted(LOESS), col = "black", lwd = 2)

## autocorrelation
ACF <- acf(RESID, lag = 379)
BP <- barplot(as.numeric(ACF$acf), col = "darkgrey", space = 2)
xPos <- BP[, 1]
FP <- findpeaks(as.numeric(ACF$acf), bw = 5)
points(xPos[FP$xmax], 1.2 * FP$ymax, pch = 25,  col = "darkred")
PERIOD <- diff(FP$xmax)

#################### Figure 3: Periodicity analysis on all datasets #################
## reps384 from CFX384
ml1 <- modlist(reps384, model = spl3, baseline = "lin", basecyc = 1:10)
CheckPeriod(ml1, cycle = 20)

## S27 from Rotorgene
ml2 <- modlist(S27, model = spl3, baseline = "lin", basecyc = 1:10)
CheckPeriod(ml2, cycle = 18)

## VIM.CFX96 from CFX96
ml3 <- modlist(VIM.CFX96, model = spl3, baseline = "lin", basecyc = 1:10)
CheckPeriod(ml3, peak.window = 5, cycle = 18)

## VIM.iQ5 from iQ5
ml4 <- modlist(VIM.iQ5, model = spl3, baseline = "lin", basecyc = 1:10)
CheckPeriod(ml4, peak.window = 5, cycle = 18)

########## Supplemental Data 2: Periodicity analysis on all datasets ##########
## GAPDH.StepOne from StepOne
ml5 <- modlist(GAPDH.StepOne, model = spl3, baseline = "lin", basecyc = 1:10)
CheckPeriod(ml5, peak.window = 5, cycle = 18)

## GAPDH.LC96.Multi from LC96 with multi-channel pipettor
ml6 <- modlist(GAPDH.LC96.Multi, model = spl3, baseline = "lin", basecyc = 1:10)
CheckPeriod(ml6, peak.window = 5, cycle = 22)

## GAPDH.LC96.Single from LC96 with single-channel pipettor
ml7 <- modlist(GAPDH.LC96.Single, model = spl3, baseline = "lin", basecyc = 1:10)
CheckPeriod(ml7, peak.window = 5, cycle = 22)

################## Figure 4: Effect of scale on Cq value #################
## create 100 curves with plateau value between 1000 and 10000
DAT <- reps384[, 1:2]
ml8 <- modlist(DAT, 1, 2, model = spl3, baseline = "lin", basecyc = 1:10)
DAT2 <- ml8[[1]]$DATA
DAT2[1, 2] <- -20  ## set first cycle to low value
SEQ <- seq(5000, 10000, by = 100)

for (i in SEQ) {
  DAT3 <- DAT2[, 2]
  DAT3 <- qpcR:::rescale(DAT3, 0, i)  ## normalize within [0, max]
  DAT2 <- cbind(DAT2, DAT3)
}
DAT2 <- DAT2[, -2]

COL <- rainbow(ncol(DAT2) -1)
matplot(DAT2[, -1], pch = 16, cex = 0.5, col = COL, type = "o", lty = 1)

## Calculate and visualize Cq at thresh 500
ml9 <- modlist(DAT2, model = l5)
Cq.thresh <- sapply(ml9, function(x) predict(x, newdata = data.frame(Fluo = 500), which = "x"))
Cq.thresh <- unlist(Cq.thresh)
par(xpd = FALSE)
matplot(DAT2[, -1], pch = 16, cex = 0.5, col = COL, type = "p", lty = 1,
        xlim = c(16, 22), ylim = c(0, 4000))
FITTED <- sapply(ml9, fitted)
matlines(FITTED, lty = 1, col = COL)
par(xpd = TRUE)
USR <- par()$usr
segments(USR[1], 500, USR[2], 500, lwd = 2)
segments(Cq.thresh, USR[3], Cq.thresh, 500, col = COL)
plot(Cq.thresh, Fmax, col = COL, pch = 16)

## Calculate and visualize Cq at SDM
Cq.SDM <- sapply(ml9, function(x) efficiency(x, type = "cpD2", plot = FALSE)$cpD2)
Cq.SDM <- unlist(Cq.SDM)
par(xpd = FALSE)
matplot(DAT2[, -1], pch = 16, cex = 0.5, col = COL, type = "p", lty = 1,
        xlim = c(16, 22), ylim = c(0, 4000))
FITTED <- sapply(ml9, fitted)
matlines(FITTED, lty = 1, col = COL)
par(xpd = TRUE)
USR <- par()$usr
segments(Cq.SDM, USR[3], Cq.SDM, 500, col = COL)
par(xpd = FALSE)

################## Figure 5: Cq versus Fmax for all datasets ##################
## reps384
ml10 <- modlist(reps384, model = l5, baseline = "lin", basecyc = 1:10)
CQ <- sapply(ml10, function(x) predict(x, newdata = data.frame(Fluo = 500), which = "x"))
CQ <- unlist(CQ)
MAX <- sapply(ml10, function(x) coef(x)[3])
plot(CQ, MAX, pch = 16, col = 1)
LM <- lm(MAX ~ CQ)
abline(LM, col = "darkgrey", lwd = 2)
summary(LM)

## GAPDH.StepOne
ml11 <- modlist(GAPDH.StepOne, model = l5, baseline = "lin", basecyc = 1:10)
CQ <- sapply(ml11, function(x) predict(x, newdata = data.frame(Fluo = 50000), which = "x"))
CQ <- unlist(CQ)
MAX <- sapply(ml11, function(x) coef(x)[3])
plot(CQ, MAX, pch = 16, col = 1)
LM <- lm(MAX ~ CQ)
abline(LM, col = "darkgrey", lwd = 2)
summary(LM)

## VIM.CFX96
ml12 <- modlist(VIM.CFX96, model = l5, baseline = "lin", basecyc = 1:10)
CQ <- sapply(ml12, function(x) predict(x, newdata = data.frame(Fluo = 250), which = "x"))
CQ <- unlist(CQ)
MAX <- sapply(ml12, function(x) coef(x)[3])
plot(CQ, MAX, pch = 16, col = 1)
LM <- lm(MAX ~ CQ)
abline(LM, col = "darkgrey", lwd = 2)
summary(LM)

## VIM.iQ5
ml13 <- modlist(VIM.iQ5, model = l5, baseline = "lin", basecyc = 1:10)
CQ <- sapply(ml13, function(x) predict(x, newdata = data.frame(Fluo = 2500), which = "x"))
CQ <- unlist(CQ)
MAX <- sapply(ml13, function(x) coef(x)[3])
plot(CQ, MAX, pch = 16, col = 1)
LM <- lm(MAX ~ CQ)
abline(LM, col = "darkgrey", lwd = 2)
summary(LM)

################### Figure 6: Different methods of Ruijter et al. (2013) ##############
## Data is in RuijterCqs
## interpolate missing values
RuijterCqsCompl <- sapply(RuijterCqs, function(x) na.spline(x))

par(mfrow = c(6, 1))
par(mar = c(0, 3, 1, 0))
statMat <- NULL
for (i in (1:7)[-3]) {
  CP <- CheckPeriod(RuijterCqsCompl[, i], peak.window = 10, plot = FALSE)  
  ACF <- acf(CP$resid, lag = length(CP$resid), plot = FALSE)
  BP <- barplot(as.numeric(ACF$acf), col = "darkgrey", space = 2)    
  title(main = sub("X", "", colnames(RuijterCqsCompl)[i]), line = -1)
  statMat <- rbind(statMat, c(CP$lmStat, CP$randStat, CP$periodStat)) 
}
rownames(statMat) <- colnames(RuijterCqsCompl)[-3]

write.table(statMat, file = "clipboard", sep = "\t")

## Data is in RuijterEffs
## interpolate missing values
RuijterEffsCompl <- sapply(RuijterEffs, function(x) na.spline(x))

par(mfrow = c(6, 1))
par(mar = c(0, 3, 1, 0))
statMat <- NULL
for (i in 1:6) {
  CP <- CheckPeriod(RuijterEffsCompl[, i], peak.window = 10, plot = FALSE)  
  ACF <- acf(CP$resid, lag = length(CP$resid), plot = FALSE)
  BP <- barplot(as.numeric(ACF$acf), col = "darkgrey", space = 2)    
  title(main = sub("X", "", colnames(RuijterEffsCompl)[i]), line = -1)
  statMat <- rbind(statMat, c(CP$lmStat, CP$randStat, CP$periodStat))  
}
rownames(statMat) <- colnames(RuijterEffsCompl)

write.table(statMat, file = "clipboard", sep = "\t")

## Data is in RuijterN0s
## interpolate missing values
RuijterN0sCompl <- sapply(RuijterN0s, function(x) na.spline(x))

par(mfrow = c(6, 1))
par(mar = c(0, 3, 1, 0))
statMat <- NULL
for (i in 1:6) {
  CP <- CheckPeriod(RuijterN0sCompl[, i], peak.window = 10, plot = FALSE)  
  ACF <- acf(CP$resid, lag = length(CP$resid), plot = FALSE)
  BP <- barplot(as.numeric(ACF$acf), col = "darkgrey", space = 2)    
  title(main = sub("X", "", colnames(RuijterN0sCompl)[i]), line = -1)
  statMat <- rbind(statMat, c(CP$lmStat, CP$randStat, CP$periodStat))  
}
rownames(statMat) <- colnames(RuijterN0sCompl)[-7]

write.table(statMat, file = "clipboard", sep = "\t")

## For mak2 model
statMat <- NULL
for (i in 7:7) {
  CP <- CheckPeriod(RuijterN0sCompl[, i], peak.window = 10, plot = FALSE)  
  ACF <- acf(CP$resid, lag = length(CP$resid), plot = FALSE)
  BP <- barplot(as.numeric(ACF$acf), col = "darkgrey", space = 2)    
  title(main = sub("X", "", colnames(RuijterN0sCompl)[i]), line = -1)
  statMat <- rbind(statMat, c(CP$lmStat, CP$randStat, CP$periodStat))  
}
rownames(statMat) <- colnames(RuijterN0sCompl)[7]

write.table(statMat, file = "clipboard", sep = "\t")

## For Cy0 model
statMat <- NULL
for (i in 3:3) {
  CP <- CheckPeriod(RuijterCqsCompl[, i], peak.window = 10, plot = FALSE)  
  ACF <- acf(CP$resid, lag = length(CP$resid), plot = FALSE)
  BP <- barplot(as.numeric(ACF$acf), col = "darkgrey", space = 2)    
  title(main = sub("X", "", colnames(RuijterCqsCompl)[i]), line = -1)
  statMat <- rbind(statMat, c(CP$lmStat, CP$randStat, CP$periodStat)) 
}
rownames(statMat) <- colnames(RuijterCqsCompl)[3]

write.table(statMat, file = "clipboard", sep = "\t")

################# Figure 7: Effect of Normalization on Cq #####################
## reps384 from CFX384
ml10a <- modlist(reps384, model = spl3, baseline = "lin", basecyc = 1:10, norm = FALSE)
Cq.10a <- unlist(sapply(ml10a, function(x) predict(x, newdata = data.frame(Fluo = 500), which = "x")))
ml10b <- modlist(reps384, model = spl3, baseline = "lin", basecyc = 1:10, norm = TRUE)
Cq.10b <- unlist(sapply(ml10b, function(x) predict(x, newdata = data.frame(Fluo = 0.1), which = "x")))
par(xpd = FALSE)
plot(ml10b)
abline(h = 0.1, col = "red", lwd = 1)
USR <- par()$usr
segments(Cq.10b, USR[3], Cq.10b, 0.1, col = "red")
CheckPeriod(Cq.10a)
CheckPeriod(Cq.10b)

## VIM.CFX96 from CFX96
ml11a <- modlist(VIM.CFX96, model = spl3, baseline = "lin", basecyc = 1:10, norm = FALSE)
Cq.11a <- unlist(sapply(ml11a, function(x) predict(x, newdata = data.frame(Fluo = 250), which = "x")))
ml11b <- modlist(VIM.CFX96, model = spl3, baseline = "lin", basecyc = 1:10, norm = TRUE)
Cq.11b <- unlist(sapply(ml11b, function(x) predict(x, newdata = data.frame(Fluo = 0.1), which = "x")))
par(xpd = FALSE)
plot(ml11b)
abline(h = 0.1, col = "red", lwd = 1)
USR <- par()$usr
segments(Cq.11b, USR[3], Cq.11b, 0.1, col = "red")
CheckPeriod(Cq.11a)
CheckPeriod(Cq.11b)

## VIM.iQ5 from iQ5
ml12a <- modlist(VIM.iQ5, model = l5, baseline = "lin", basecyc = 1:10, norm = FALSE)
Cq.12a <- unlist(sapply(ml12a, function(x) predict(x, newdata = data.frame(Fluo = 2500), which = "x")))
ml12b <- modlist(VIM.iQ5, model = l5, baseline = "lin", basecyc = 1:10, norm = TRUE)
Cq.12b <- unlist(sapply(ml12b, function(x) predict(x, newdata = data.frame(Fluo = 0.1), which = "x")))
par(xpd = FALSE)
plot(ml12b)
abline(h = 0.1, col = "red", lwd = 1)
USR <- par()$usr
segments(Cq.12b, USR[3], Cq.12b, 0.1, col = "red")
CheckPeriod(Cq.12a)
CheckPeriod(Cq.12b)

## GAPDH.StepOne from StepOne
ml13a <- modlist(GAPDH.StepOne, model = l5, baseline = "lin", basecyc = 1:10, norm = FALSE)
Cq.13a <- unlist(sapply(ml13a, function(x) predict(x, newdata = data.frame(Fluo = 2500), which = "x")))
ml13b <- modlist(GAPDH.StepOne, model = l5, baseline = "lin", basecyc = 1:10, norm = TRUE)
Cq.13b <- unlist(sapply(ml13b, function(x) predict(x, newdata = data.frame(Fluo = 0.1), which = "x")))
par(xpd = FALSE)
plot(ml13b)
abline(h = 0.1, col = "red", lwd = 1)
USR <- par()$usr
segments(Cq.13b, USR[3], Cq.13b, 0.1, col = "red")
CheckPeriod(Cq.13a)
CheckPeriod(Cq.13b)

####################### Figure 8: Mapping of Cq's to MTP #############################
COL <- colorRampPalette(c("darkblue", "white", "darkred"))
## reps384
ml30 <- modlist(reps384, model = l5, baseline = "lin", basecyc = 1:10)
CQ30 <- sapply(ml30, function(x) predict(x, newdata = data.frame(Fluo = 500), which = "x"))
CQ30 <- unlist(CQ30)
RESID <- CheckPeriod(CQ30)$resid
MAT <- matrix(RESID, ncol = 24, byrow = TRUE)
colnames(MAT) <- as.character(1:24)
rownames(MAT) <- LETTERS[1:16]
levelplot(t(MAT[nrow(MAT):1, ]), col.regions = COL(100))

## CFX96
ml31 <- modlist(VIM.CFX96, model = spl3, baseline = "lin", basecyc = 1:10)
CQ31 <- sapply(ml31, function(x) predict(x, newdata = data.frame(Fluo = 250), which = "x"))
CQ31 <- unlist(CQ31)
RESID <- CheckPeriod(CQ31)$resid
MAT <- matrix(RESID, ncol = 12, byrow = TRUE)
colnames(MAT) <- as.character(1:12)
rownames(MAT) <- LETTERS[1:8]
levelplot(t(MAT[nrow(MAT):1, ]), col.regions = COL(100))

## StepOne
ml32 <- modlist(GAPDH.StepOne, model = l5, baseline = "lin", basecyc = 1:10)
CQ32 <- sapply(ml32, function(x) predict(x, newdata = data.frame(Fluo = 50000), which = "x"))
CQ32 <- unlist(CQ32)
RESID <- CheckPeriod(CQ32)$resid
MAT <- matrix(RESID, ncol = 12, byrow = TRUE)
colnames(MAT) <- as.character(1:12)
rownames(MAT) <- LETTERS[1:8]
levelplot(t(MAT[nrow(MAT):1, ]), col.regions = COL(100))

## LC96 multichannel
ml33 <- modlist(GAPDH.LC96.Single, model = spl3, baseline = "lin", basecyc = 1:10)
CQ33 <- sapply(ml33, function(x) predict(x, newdata = data.frame(Fluo = 0.1), which = "x"))
CQ33 <- unlist(CQ33)
RESID <- CheckPeriod(CQ33)$resid
MAT <- matrix(RESID, ncol = 12, byrow = TRUE)
colnames(MAT) <- as.character(1:12)
rownames(MAT) <- LETTERS[1:8]
levelplot(t(MAT[nrow(MAT):1, ]), col.regions = COL(100))
