# 04-refIntensity-dependencies.R
# Examine dependencies of reference intensities (rho's) within constant c(b,l)
#   for given parameter estimates
#
# Input: analysisFunctions.R; analysisResults.RData
# Output: plots/constDepen.pdf
# Packages: - 
#
# Last mod: 12/29/2025, AO

# Analysis functions
source("analysisFunctions.R")
# Get fitted models
load("analysisResults.RData")   # available as "testResults"

# Table for parameter estimates
# m0_par <- do.call(rbind, 
#                   lapply(testResults, function(x) x$selectedModels$m0$pars))
m1_par <- do.call(rbind, 
                  lapply(testResults, function(x) x$selectedModels$m1$pars))

rho_b <- function(alpha_l, beta_b, beta_l, w_1, c_bl, rhoL) {
    db_lambert(((c_bl - alpha_l*db_inv_spl(rhoL)^beta_l) / -w_1) ^(1/beta_b))
}
rho_l <- function(alpha_l, beta_b, beta_l, w_1, c_bl, rhoB) {
    db_spl(((c_bl + w_1*db_inv_lambert(rhoB)^beta_b) / alpha_l) ^ (1/beta_l))
}

pdf("plots/constDepen.pdf", height=4, width=6, pointsize=10)
    par(mfrow=c(2, 3), mgp=c(2, .7, 0), mar=c(0, 0, 0, 0), oma=c(4, 4.5, 1, 1))
    for(i in 1:nrow(m1_par)) {
        b1 <- with(m1_par[i, ], 
            rho_b(alpha_l, beta_b, beta_l, w_1, m1_par[i, "const_bl1"],
                seq(0, 90, by=0.2)))
        l1 <- with(m1_par[i, ], 
            rho_l(alpha_l, beta_b, beta_l, w_1, m1_par[i, "const_bl1"], 
                seq(40, 90, by=0.2)))
        b2 <- with(m1_par[i, ], 
            rho_b(alpha_l, beta_b, beta_l, w_1, m1_par[i, "const_bl2"],
                seq(0, 90, by=0.2)))
        l2 <- with(m1_par[i, ], 
            rho_l(alpha_l, beta_b, beta_l, w_1, m1_par[i, "const_bl2"], 
                seq(40, 90, by=0.2)))
        b3 <- with(m1_par[i, ], 
            rho_b(alpha_l, beta_b, beta_l, w_1, m1_par[i, "const_bl3"],
                seq(0, 90, by=0.2)))
        l3 <- with(m1_par[i, ], 
            rho_l(alpha_l, beta_b, beta_l, w_1, m1_par[i, "const_bl3"], 
                seq(40, 90, by=0.2)))
        plot(seq(0, 90, by=0.2) ~ b1, type="l", xlim=c(48, 90), ylim=c(0, 100),
            axes=FALSE, xlab="", ylab="", col="slateblue3", lwd=2)
        title(paste0("Subj. 0", i), line=-1.2)
        points(seq(0, 90, by=0.2) ~ b2, type="l", col="mediumaquamarine", lwd=2)
        points(seq(0, 90, by=0.2) ~ b3, type="l", col="forestgreen", lwd=2)
        points(l1 ~ seq(40, 90, by=0.2), type="l", col="slateblue3", lwd=2)
        points(l2 ~ seq(40, 90, by=0.2), type="l", col="mediumaquamarine", lwd=2)
        points(l3 ~ seq(40, 90, by=0.2), type="l", col="forestgreen", lwd=2)
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
