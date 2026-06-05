# 05-bayesianEstimation.R
# Bayesian estimation of global psychophysics model
#
# Input: PSEs.txt; gpm_ranges.stan, gpm_bayesResults.RData
# Output: plots/results-bayes.pdf
# Packages: rstan 
#
# Last mod: 12/30/2025, AO

library(rstan)
rstan_options(auto_write=TRUE)
rstan_options(threads_per_chain=1)

dat <- read.table("PSEs.txt", header=TRUE)
dat$range <- factor(dat$range, levels=c("low", "medium", "high"))

m_i <- lapply(split(dat, dat$subj), function(x) {
  sds <- aggregate(yhat ~ range + stan, x, sd)
  sds$range <- factor(sds$range, levels=c("low", "medium", "high"))
  sds <- sds[order(sds$range), ]
  m_ind <- stan(file="gpm_ranges.stan", 
             data=list(ntotal=nrow(x), 
                       nx_bl=length(unique(paste0(x$stan, "_", x$range))),
                       nranges=length(unique(x$range)),
                       x_bl=x$stan,
                       y_bl=x$yhat,
                       x_bl_idx=rep(1:15, each=4),
                       range_idx=as.numeric(x$range),
                       sig_bl=rep(sds$yhat, each=4),
                       rho_b=60   # fix light referency intensity at arbitrary 
                                  #    (perceptible) value for estimation
             )
  )
  cat(paste0("\n\ncompleted subject 0", unique(x$subj), "\n\n"))
  m_ind
})

save(m_i, file="gpm_bayesResults.RData")
# load("gpm_bayesResults.RData")  # available as "m_i"

# View estimates and diagnostic summary
out <- lapply(m_i, function(x) {
   print(x, pars=c("alpha_l", "beta_b", "beta_l", "w_1", "const_bl", "mu_bl",
                   "rho_l_db"), prob=c(.025, .975), digits_summary=2)
})


# Plot data and predictions (compared to ML estimation)
source("analysisFunctions.R")
load("analysisResults.RData")   # available as "testResults"
dat$subj <- factor(dat$subj)

pdf("plots/results-bayes.pdf", pointsize=10, height=8, width=3.5)
par(mfrow=c(length(levels(dat$subj)), 2), mar=c(0, 0, 0, 0), mgp=c(2, .7, 0), 
    oma=c(4, 7, 3, 1))
for(i in levels(dat$subj)) {
  dat_i <- dat[dat$subj == i, ]
  plot_gpm(dat_i, unlist(testResults[[as.numeric(i)]]$selectedModels$m1$pars), 
           ranges=TRUE)
  if(i == "1") {
    legend("topleft", bty="n", pch=c(16, 21, 17, NA), 
           legend=c("low", "medium ", "high"))
  }
  axis(2, at=56.5, tick=FALSE, labels=paste0("Subj. 0", i), font=2, outer=TRUE,
       line=2, cex=1.5)
  if(i == "1") {
      axis(3, at=73, tick=FALSE, "MLE", font=2, outer=TRUE, line=0.1, cex=2)
  }
  mtext("Standard (dB lambert)", side=1, line=2.3, outer=TRUE)
  mtext("PSE (dB SPL)", side=2, line=4.3, outer=TRUE)
  # Bayesian:
  plot_gpm(dat_i, c(sapply(extract(m_i[[as.numeric(i)]]), mean)[1:4],
                    setNames(colMeans(extract(m_i[[as.numeric(i)]])$const_bl), 
                             paste0("const_bl", 1:3))),
           ranges=TRUE)
  if(i == "1") {
      axis(3, at=72, tick=FALSE, "Bayesian", font=2, outer=TRUE, 
          line=0.1, cex=2)
  }
  if(i == "5") {
      legend("bottomright", bty="n", legend=c(expression(c[bl]^low), 
          expression(c[bl]^medium), expression(c[bl]^high)),
          lty=1, col=c("slateblue3", "mediumaquamarine", "forestgreen"))
  }
}
dev.off()


# Plot dependencies of reference intensities
source("analysisFunctions.R")
rho_b <- function(alpha_l, beta_b, beta_l, w_1, c_bl, rhoL) {
    db_lambert(((c_bl - alpha_l*db_inv_spl(rhoL)^beta_l) / -w_1) ^(1/beta_b))
}
rho_l <- function(alpha_l, beta_b, beta_l, w_1, c_bl, rhoB) {
    db_spl(((c_bl + w_1*db_inv_lambert(rhoB)^beta_b) / alpha_l) ^ (1/beta_l))
}

m_par <- as.data.frame(t(rbind(sapply(1:5, function(i) {
    c(sapply(extract(m_i[[as.numeric(i)]]), mean)[1:4], 
      setNames(colMeans(extract(m_i[[as.numeric(i)]])$const_bl), 
               paste0("const_bl", 1:3)))
  })
)))


pdf("plots/constDepen-bayes.pdf", height=4, width=6, pointsize=10)
par(mfrow=c(2, 3), mgp=c(2, .7, 0), mar=c(0, 0, 0, 0), oma=c(4, 4.5, 1, 1))
for(i in 1:nrow(m_par)) {
  rho_b1 <- with(m_par[i, ], rho_b(alpha_l, beta_b, beta_l, w_1, const_bl1, 
                 rhoL=seq(0, 110, by=.2)))
  rho_b2 <- with(m_par[i, ], rho_b(alpha_l, beta_b, beta_l, w_1, const_bl2,
                 rhoL=seq(0, 110, by=.2)))
  rho_b3 <- with(m_par[i, ], rho_b(alpha_l, beta_b, beta_l, w_1, const_bl3,
                 rhoL=seq(0, 110, by=.2)))
  plot(seq(0, 110, by=0.2) ~ rho_b1, type="l", xlim=c(48, 90), ylim=c(0, 100),
      axes=FALSE, xlab="", ylab="", col="slateblue3", lwd=2)
  title(paste0("Subj. 0", i), line=-1.2)
  points(seq(0, 110, by=0.2) ~ rho_b2, type="l", col="mediumaquamarine", lwd=2)
  points(seq(0, 110, by=0.2) ~ rho_b3, type="l", col="forestgreen", lwd=2)
  box()
  axis(1, at=seq(50, 90, by=5), labels=FALSE, tcl=-.2)
  axis(2, at=seq(0, 110, by=10), labels=FALSE, tcl=-.2)
  axis(2, at=seq(0, 110, by=20), labels=FALSE)
  if(i == 1 | i == 4) axis(2)
  if(i == 3 | i == 4 | i == 5) axis(1)
}
mtext(expression(paste(rho["l<-b"], " (dB SPL)")), side=2, line=2.7, 
      outer=TRUE)
mtext(expression(paste(rho["b->l"], " (dB Lambert)")), side=1, line=2.7, 
      outer=TRUE)
plot(0, 0, xlim=c(50, 90), ylim=c(0, 100), axes=FALSE, xlab="", ylab="")
legend(60, 80, bty="n", legend=c(expression(c[bl]^low), 
       expression(c[bl]^medium), expression(c[bl]^high)),
       lty=1, col=c("slateblue3", "mediumaquamarine", "forestgreen"), cex=1.5, 
       lwd=2)
dev.off()
