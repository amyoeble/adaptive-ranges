# power_ranges.R
# Power global psychophysics model (Heller, 2020)
# Range effect based on Naumann's (2023) MLEs of matching data, subject 04
#
# Input: ./../analysis/analysisFunctions.R
# Output: simulation_example.pdf
# Packages: -
#
# Last mod: 01/02/2026, AO

source("./../analysis/analysisFunctions.R")

# Parameter estimates for Subj04 (from matching data of Naumann, 2021)
par04 <- list(w_1=0.68, alpha=17.06, beta_l=0.15, beta_b=0.35, const_bl=5.75)

# Set sample size (PSEs per standard), SD is based on pilot (-> something in 
#   between 2 and 0.8), and minimally relevant effect (subjective JND ~ 3 dB)
n <- 4; se <- 2; mineff <- 2

pred <- data.frame(
    # Light standards
    stan=c(seq(71, 79, by=2), seq(77, 85, by=2), seq(83, 91, by=2)),
    range=factor(rep(1:3, each=5), labels=c("low", "medium", "high"))
)
# Predicted points (for Subj04 in matching)
pred$pred0 <- bright_to_loud(pred$stan, par04$alpha, par04$beta_l, 
        par04$beta_b, par04$w_1, par04$const_bl)
# Add "normalization effect" of 2 dB
pred$effpred <- pred$pred0 + c(rep(mineff, 5), rep(0, 5), rep(-mineff, 5))

# Simulate power (caution: loooong runtime, ~30 replicates/h)
# Width of CI for simulated power is: 
# qnorm(.95) * sqrt((0.8*0.2)/nRepl) * 2    # -> 0.0701 for nRepl=500

system.time(pval <- replicate(1000, {
    # Use "effpred" (predictions with effect) as mu to simulate new PSEs from 
    # normal distribution with SD from pilot
    dat <- data.frame(stan=rep(pred$stan, each=n),
                      range=rep(pred$range, each=n),
                      yhat=unlist(lapply(pred$effpred, rnorm, n=n, sd=se)))
    # fit restricted and unrestricted model, with multiple start values
    m0 <- llFit(dat, ranges=FALSE, nStarts=50)
    m1 <- llFit(dat, ranges=TRUE, nStarts=50)
    # select best model fit with minimal NLL, calculate difference in likelihood
    selectedModels <- list(m0=llSelect(m0), m1=llSelect(m1))
    ldiff <- abs(diff(c(selectedModels$m0$value, selectedModels$m1$value)))
    # Likelihood ratio test, determine p-value (2 DF)
    1 - pchisq(2*ldiff, 2)
}))
mean(pval < .05)    
# Minimally relevant effect of 2 dB -> n=4: 0.876 (n=5: 0.934)

# To visualize simulated PSES and models for one loop:
dat <- data.frame(stan=rep(pred$stan, each=n),
                  range=rep(pred$range, each=n),
                  yhat=unlist(lapply(pred$effpred, rnorm, n=n, sd=se)))
# fit restricted and unrestricted model, with multiple start values
m0 <- llFit(dat, ranges=FALSE, nStarts=50)
m1 <- llFit(dat, ranges=TRUE, nStarts=50)
# select best model fit with minimal NLL, calculate difference in likelihood
selectedModels <- list(m0=llSelect(m0), m1=llSelect(m1))


pdf("simulation-example.pdf", pointsize=10, height=3.5, width=7)
par(mfrow=c(1, 2), mar=c(3, 3, 1, 1), mgp=c(2, .7, 0))
pse_hat <- aggregate(yhat ~ stan + range, dat, mean)
# M0:
# Add data points
plot(yhat ~ stan, dat[dat$range=="low", ], pch=16, col="gray", xlim=c(70, 92), 
    ylim=c(40, 75), xlab="Standard (dB lambert)", ylab="Simulated PSE (dB SPL)")
points(yhat ~ stan, dat[dat$range=="medium", ], pch=21, col="gray")
points(yhat ~ stan, dat[dat$range=="high", ], pch=17, col="gray")
# Add estimated PSEs
points(yhat ~ stan, pse_hat[pse_hat$range=="low", ], pch=16, cex=1.5)
points(yhat ~ stan, pse_hat[pse_hat$range=="medium", ], pch=21, cex=1.5)
points(yhat ~ stan, pse_hat[pse_hat$range=="high", ], pch=17, cex=1.5)
p0 <- with(selectedModels$m0$pars,
    bright_to_loud(seq(65, 95, by=.2), alpha_l, beta_l, beta_b, w_1, const_bl))
points(p0 ~ seq(65, 95, by=.2), type="l", col="maroon", lwd=2)
legend("topleft", bty="n", pch=c(16, 21, 17, NA), 
       legend=c("low", "medium ", "high"))
# M1:
# Add data points
plot(yhat ~ stan, dat[dat$range=="low", ], pch=16, col="gray", xlim=c(70, 92), 
    ylim=c(40, 75), xlab="Standard (dB lambert)", ylab="Simulated PSE (dB SPL)")
points(yhat ~ stan, dat[dat$range=="medium", ], pch=21, col="gray")
points(yhat ~ stan, dat[dat$range=="high", ], pch=17, col="gray")
# Add PSE estimates
points(yhat ~ stan, pse_hat[pse_hat$range=="low", ], pch=16, cex=1.5)
points(yhat ~ stan, pse_hat[pse_hat$range=="medium", ], pch=21, cex=1.5)
points(yhat ~ stan, pse_hat[pse_hat$range=="high", ], pch=17, cex=1.5)
p1_1 <- with(selectedModels$m1$pars,
    bright_to_loud(seq(65, 95, by=.2), alpha_l, beta_l, beta_b, w_1, const_bl1))
points(p1_1 ~ seq(65, 95, by=.2), type="l", col="slateblue3", lwd=2)
p1_2 <- with(selectedModels$m1$pars,
    bright_to_loud(seq(65, 95, by=.2), alpha_l, beta_l, beta_b, w_1, const_bl2))
points(p1_2 ~ seq(65, 95, by=.2), type="l", col="mediumaquamarine", lwd=2)
p1_3 <- with(selectedModels$m1$pars,
    bright_to_loud(seq(65, 95, by=.2), alpha_l, beta_l, beta_b, w_1, const_bl3))
points(p1_3 ~ seq(65, 95, by=.2), type="l", col="forestgreen", lwd=2)
legend("bottomright", bty="n", legend=c(expression(c[bl]^low), 
    expression(c[bl]^medium), expression(c[bl]^high)),
    lty=1, col=c("slateblue3", "mediumaquamarine", "forestgreen"))
dev.off()

