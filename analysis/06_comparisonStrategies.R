# 06-comparisonStrategies.R
# Check whether participants compare consecutive sounds across trials instead
# of cross-modal intensity within trials
#
# Input: ../data/subj*.txt
# Output: 
# Packages: - 
#
# Last mod: 05/19/2026, AO

## SHORT PREPROCESSING
# Read data
dat0 <- do.call(rbind, 
    lapply(dir("../data/", "Subj.*\\.txt", full.names=TRUE),
        read.table, header=TRUE, as.is=FALSE))
dat0$range <-  factor(dat0$range, levels=c("low", "medium", "high"))
dat0$subj <- factor(dat0$subj)

# Exclude training trials
dat <- dat0[dat0$stair != "training", ]

# Flat matching curves appear especially for Subject 02's and Subject 03's 
#   medium and high session
# Calculate consecutive sound differences
d <- lapply(split(dat, dat$subj), function(d_i) {
        do.call(rbind, lapply(split(d_i, d_i$range), function(d_ir) {
            d_ir$sound_diff=c(NA, diff(d_ir$comparison))
            # Exclude first trials after pauses (every 20th trial)
            d_ir$trialIndex <- 1:nrow(d_ir)
            # Exclude 'same' responses for GLM
            d_ir[!is.na(d_ir$sound_diff) & d_ir$response != "space" &
                 d_ir$trialIndex %% 20 != 1, ]
        }))
})        

glm_i <- lapply(d, function(d_i) {
       glm(response ~ sound_diff*range, d_i, family=binomial(logit))
})
lapply(glm_i, summary)
# Generally: the larger the difference in consecutive sounds, the larger the 
#   likelihood of responding 'j' ('loudness more intense')
lapply(d, \(x) aggregate(sound_diff ~ range, x, quantile, 
                         c(0, .05, .25, .75, .95, 1), digits=3))
# Even the seemingly small coefficients cause pretty large effects when 
#   considering that consecutive sounds can differ A LOT (Q1 and Q3 look less
#   alarming)

# Calculate actual effects (in logits) for each subject and session
c_i <- lapply(glm_i, function(g_i) {
    coef_idx <- grep("sound_diff", names(coef(g_i)))
    c_i <- cbind(estimate=coef(g_i), confint(g_i))[coef_idx, ]
    apply(c_i, 2, function(est) {
      c(low=est[1], medium=sum(est[1:2]), high=sum(est[c(1, 3)]))
    }) |> as.data.frame()
})
for(i in 1:5) {
    c_i[[i]]$Subject <- paste0("0", i)
    c_i[[i]]$Range <- factor(levels(dat$range), levels=levels(dat$range))
}
glm_coefs <- do.call(rbind, c_i)[, c(4:5, 1:3)]
rownames(glm_coefs) <- NULL
glm_coefs[, 3:5] <- round(glm_coefs[, 3:5], 3)
glm_coefs$Estimate <- paste0(format(glm_coefs[, 3], digits=3), " [", 
                             format(glm_coefs[, 4], digits=3), ", ", 
                             format(glm_coefs[, 5], digits=3), "]")
glmTab <- reshape(glm_coefs[, c("Subject", "Range", "Estimate")], 
                  direction="wide", timevar="Range", idvar="Subject", 
                  varying=paste(c("Low", "Medium", "High"), "Range"))

# Calculate bias for actually observed sound differences
d <- lapply(d, function(d_i) {
    d_i$responseBias <- exp(d_i$sound_diff * 
        c_i[[as.numeric(d_i$subj)[1]]]$estimate[as.numeric(
                                        d[[as.numeric(d_i$subj)[1]]]$range)])
    d_i
})
# Observed SPL differences and estimated effects
lapply(d, \(x) rbind(quantile(x$responseBias, c(0, .05, .95, 1)),
                     quantile(x$sound_diff, c(0, .05, .95, 1))) |> round(3))
