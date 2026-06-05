# 03-resultPlots.R
# Plot results of matching ranges data analysis - models from Heller's
#   (2020) framework, model comparison to test a range effect
#
# Input: analysisFunctions.R, PSEs.txt, analysisResults.RData
# Output: plots/results-adapRanges.pdf
# Packages: - 
#
# Last mod: 12/29/2025, AO

# Get necessary functions, data, and fitted models
source("analysisFunctions.R")
dat <- read.table("PSEs.txt", header=TRUE, as.is=FALSE)
dat$range <- factor(dat$range, levels=c("low", "medium", "high"))
dat$subj <- factor(dat$subj)
load("analysisResults.RData")   # available as "testResults"


pdf("plots/results-adapRanges.pdf", pointsize=10, height=8, width=3.5)
par(mfrow=c(length(levels(dat$subj)), 2), mar=c(0, 0, 0, 0), mgp=c(2, .7, 0), 
    oma=c(4, 7, 3, 1))
for(i in levels(dat$subj)) {
  dat_i <- dat[dat$subj == i, ]
  # M0
  plot_gpm(dat_i, unlist(testResults[[as.numeric(i)]]$selectedModels$m0$pars),
           ranges=FALSE) 
  if(i == "1") {
    legend("topleft", bty="n", pch=c(16, 21, 17, NA), 
           legend=c("low", "medium ", "high"))
  }
  axis(2, at=56.5, tick=FALSE, labels=paste0("Subj. 0", i), font=2, outer=TRUE,
       line=2, cex=1.5)
  if(i == "1") {
    axis(3, at=73, tick=FALSE, expression(M[0]), font=2, outer=TRUE, line=0.1, 
         cex=2)
  }
  # M1:
  plot_gpm(dat_i, unlist(testResults[[as.numeric(i)]]$selectedModels$m1$pars), 
           ranges=TRUE) 
  if(i == "1") {
      axis(3, at=73, tick=FALSE, expression(M[1]), font=2, outer=TRUE, 
          line=0.1, cex=2)
  }
  mtext("Standard (dB lambert)", side=1, line=2.3, outer=TRUE)
  mtext("PSE (dB SPL)", side=2, line=4.3, outer=TRUE)
}
dev.off()


# Determine effect sizes (in dB SPL):
pse <- aggregate(yhat ~ range + stan + subj, dat, mean)
effSizes <- lapply(split(pse, pse$subj), function(x) {
    junct <- x[duplicated(x$stan) != duplicated(x$stan, fromLast=TRUE), ]
    sapply(split(junct, junct$stan), function(x2) diff(x2$yhat))
})
summary(unlist(effSizes))
