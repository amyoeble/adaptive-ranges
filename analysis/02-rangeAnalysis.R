# 02-rangeAnalysis.R
# Analysis of matching ranges data - fit and comparison of models from Heller's
#   (2020) framework for testing of range effect
#
# Input: analysisFunctions.R, PSEs.txt
# Output: analysisResults.RData
# Packages: - 
#
# Last mod: 12/29/2025, AO

# Get necessary functions
source("analysisFunctions.R")

# Get data
dat <- read.table("PSEs.txt", header=TRUE, as.is=FALSE)
dat$range <- factor(dat$range, levels=c("low", "medium", "high"))
dat$subj <- factor(dat$subj)

# Analysis
system.time(testResults <- lapply(split(dat, dat$subj), function(data) {
    # fit restricted and unrestricted model, with multiple start values
    m0 <- llFit(data, ranges=FALSE, nStarts=500)
    m1 <- llFit(data, ranges=TRUE, nStarts=500)
    # select best model fit with minimal NLL, calculate difference in likelihood
    selectedModels <- list(m0=llSelect(m0), m1=llSelect(m1))
    ldiff <- abs(diff(c(selectedModels$m0$value, selectedModels$m1$value)))
    # Likelihood ratio test, determine p-value (2 DF)
    list(selectedModels=selectedModels, ldiff=ldiff, 
         pval=1 - pchisq(2*ldiff, 2), allModels=list(m0=m0, m1=m1))
}))

selectedModels <- lapply(testResults, "[[", "selectedModels")

# Save results to R object
save(testResults, file="analysisResults.RData")
