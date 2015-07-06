# library(qpcR)
# library(lmtest)
# library(lawstat)
# library(chipPCR)
# library(lattice)
# library(ggplot2)

source("findpeaks.R")
source("CheckPeriod.R")

ml2 <- modlist(reps384, model = spl3, basecyc = 1:10, baseline = "lin")
DATA <- sapply(ml2, function(x) x$DATA[, 2])

#Fig 2

#plot B
LM <- lm(DATA[20, ] ~ I(1:379) + I((1:379)^2))
FITTED <- fitted(LM)

#plot C
RESID <- residuals(LM)
LOESS <- loess(RESID ~ I(1:379), span = 0.1)

#plot D
ACF <- acf(RESID, lag = 379)



COL <- rainbow(length(ml2))
theme_app <- theme(plot.background=element_rect(fill = "transparent",
                                                colour = "transparent"),
                   plot.margin = unit(c(1,1,1,1), "cm"),
                   axis.text.x = element_text(size=7), 
                   axis.text.y = element_text(size=7),
                   axis.title.x = element_text(size=11, vjust = -1), 
                   axis.title.y = element_text(size=11, vjust = 1),
                   strip.text = element_text(size=8),
                   legend.text = element_text(size=7), 
                   legend.title = element_text(size=11),
                   plot.title = element_text(size=17),
                   panel.grid.major = element_line(colour="grey"),
                   panel.grid.minor = element_line(colour="lightgrey", linetype = "dashed"),
                   panel.background = element_rect(fill = "transparent",colour = "black"),
                   legend.background = element_rect(fill="NA"))


CheckPeriod(DATA[20, ], cycle = 20, plot = FALSE)

#plot B
plot(1:379, DATA[20, ], col = COL, pch = 16, xlab = "", ylab = "", las = 1)
lines(1:379, FITTED, col = 1, lwd = 2)

#plot C
plot(1:379, RESID, col = COL, pch = 16, xlab = "", ylab = "", las = 1)
lines(fitted(LOESS), col = "black", lwd = 2)

#plot D
ACF.dat <- data.frame(pos = 1L:dim(ACF[["acf"]])[1], acf = ACF[["acf"]][, , 1], peak = rep(NA, dim(ACF[["acf"]])[1]))
peak.list <- findpeaks(ACF.dat[["acf"]], bw = 5)
ACF.dat[peak.list[["xmax"]], "peak"] <- 1.1*peak.list[["ymax"]]


ggplot(ACF.dat, aes(x = pos, y = acf)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous("") +
  scale_y_continuous("Autocorrelation of residuals") +
  theme_app +
  ggtitle("Autocorrelation plot") +
  geom_point(aes(y = peak), pch = 25, colour = "darkred", size = 4)

findpeaks(as.numeric(ACF$acf), bw = 5)

ACF.dat[["acf"]]
BP <- barplot(as.numeric(ACF$acf), col = "darkgrey", space = 2)
xPos <- BP[, 1]
FP <- findpeaks(as.numeric(ACF$acf), bw = 5)
points(xPos[FP$xmax], 1.2 * FP$ymax, pch = 25,  col = "darkred")
PERIOD <- diff(FP$xmax)