# analysisFunctions.R
# Functions for model fitting and testing for adaptive range data
#
# Input: -
# Output: -
#
# Last mod: 12/29/2025, AO

# Convert cd/m2 to dB Lambert
db_lambert <- function(x_b) 10 * log10(x_b * pi / (10^(-6)))
# Convert sound pressure in Pa to dB SPL re 20 muPa
db_spl <- function(x_l) 20 * log10(x_l / (2 * 10^(-5)))
# Convert dB Lambert to cd/m2
db_inv_lambert <- function(db_b) 10^(db_b / 10) * 10^(-6) / pi
# Convert dB SPL to Pa
db_inv_spl <- function(db_l) 10^(db_l / 20) * 2 * 10^-5

# The model (for light standards, sound comparison)
# x_bright: standard light in dB lambert
# alpha: alpha_l; alpha_b = 1
# const_bl: - w_1 * alpha_b * db_inv_lambert(rho.sbl)^beta_b 
#            + alpha_l * db_inv_spl(rho.cbl)^beta_l)
bright_to_loud <- function(x_bright, alpha_l, beta_l, beta_b, w_1, const_bl) {
  db_spl(
     (w_1 * 1 / alpha_l * db_inv_lambert(x_bright)^beta_b +
      1 / alpha_l * const_bl)^(1 / beta_l)
  )
}

# Negative log-likelihood (also just brightness to loudness):
# NOTE: nll-functions require data argument to have a column with factor "range"
# (low vs medium vs high), "stan" (standard) and "yhat" (PSE estimates)
nll <- function(param, dat, ranges=TRUE) {
  # Correct number of constants according to range specification
  if (ranges == TRUE) {
    range_idx <- as.numeric(dat$range)
  } else range_idx <- NULL
  # Estimate sigma from data
  sig_bl <- round(rep(aggregate(yhat ~ stan + range, dat, sd)$yhat, 
                      each=nrow(dat)/15), 2)
  # Calculate likelihood
    -sum(dnorm(dat$yhat,
         mean=bright_to_loud(dat$stan, param["alpha_l"], param["beta_l"], 
           param["beta_b"], param["w_1"], param[paste0("const_bl", range_idx)]),
         sd=sig_bl, log=TRUE)
    )
}

# Draw set of start values for optimization:
# NOTE: With ranges=TRUE, equal starting values are set for all constants, 
# three ranges are hardcoded
draw_start_par <- function(ranges=TRUE) {
  w_1 <- runif(1, 0.5, 1.5)
  alpha_l <- runif(1, 0, 20)    # alpha_b = 1
  beta_l <- rbeta(1, 6, 3)
  beta_b <- rbeta(1, 3, 6)
  rho_sbl <- runif(ifelse(ranges == TRUE, 3, 1), 50, 80)
  rho_cbl <- runif(ifelse(ranges == TRUE, 3, 1), 20, 80)
  const_bl <- - w_1 * db_inv_lambert(rho_sbl)^beta_b +
                        alpha_l * db_inv_spl(rho_cbl)^beta_l
  unlist(list(w_1=w_1, alpha_l=alpha_l, beta_l=beta_l, beta_b=beta_b, 
              const_bl=const_bl)
  )
}

# Fit both models with different start values
llFit <- function(data, ranges=TRUE, nStarts=100) {
  m <- lapply(1:nStarts, function(i) {
      start_par <- draw_start_par(ranges=ranges)
      out <- tryCatch(
          optim(par=start_par, fn=nll, gr=NULL, data, ranges=ranges, 
                method="Nelder-Mead", control=list(maxit=2000)),
          error=function(e) list(NA))
      out$startvals <- start_par
      out
  })
  # use successful runs only
  m <- m[unlist(lapply(m, function(x) !is.null(x$par)))]
  list(values=unlist(lapply(m, "[[", "value")),
       pars=data.frame(do.call(rbind, lapply(m, "[[", "par"))),
       startvals=data.frame(do.call(rbind, lapply(m, "[[", "startvals")))
  )
}

# Select and return best model
llSelect <- function(Mfit) {
  idx <- which.min(Mfit$values)
  list(value=Mfit$values[idx], pars=Mfit$pars[idx, ])
}

# Visualize data and model
plot_gpm <- function(data, pars, ranges=TRUE) {
  pse_hat <- pse_hat <- aggregate(yhat ~ range + stan, data, mean)
  # Add data points
  plot(yhat ~ stan, data[data$range=="low", ], pch=16, col="dimgray", 
       xlim=c(58, 86), ylim=c(18, 95), axes=FALSE)
  axis(2, at=seq(25, 85, by=15)); axis(1); box()
  axis(2, at=seq(20, 95, by=5), labels=FALSE, tcl=-0.2)
  points(yhat ~ stan, data[data$range=="medium", ], pch=21, col="dimgray")
  points(yhat ~ stan, data[data$range=="high", ], pch=17, col="dimgray")
  # Add estimated PSEs
  points(yhat ~ stan, pse_hat[pse_hat$range=="low", ], pch=16, cex=1.5)
  points(yhat ~ stan, pse_hat[pse_hat$range=="medium", ], pch=21, cex=1.5)
  points(yhat ~ stan, pse_hat[pse_hat$range=="high", ], pch=17, cex=1.5)
  # Add model predictions
  if(ranges == FALSE) {
    p0 <- bright_to_loud(seq(55, 90, by=.2), pars["alpha_l"], pars["beta_l"], 
                         pars["beta_b"], pars["w_1"], pars["const_bl"])
    points(seq(55, 90, by=.2), p0, type="l", col="maroon", lwd=2)
  } else {
    p1_1 <- bright_to_loud(seq(55, 90, by=.2), pars["alpha_l"], pars["beta_l"], 
                         pars["beta_b"], pars["w_1"], pars["const_bl1"])
    points(p1_1 ~ seq(55, 90, by=.2), type="l", col="slateblue3", lwd=2)
    p1_2 <- bright_to_loud(seq(55, 90, by=.2), pars["alpha_l"], pars["beta_l"], 
                         pars["beta_b"], pars["w_1"], pars["const_bl2"])
    points(p1_2 ~ seq(55, 90, by=.2), type="l", col="mediumaquamarine", lwd=2)
    p1_3 <- bright_to_loud(seq(55, 90, by=.2), pars["alpha_l"], pars["beta_l"], 
                         pars["beta_b"], pars["w_1"], pars["const_bl3"])
    points(p1_3 ~ seq(55, 90, by=.2), type="l", col="forestgreen", lwd=2)
  }
}
