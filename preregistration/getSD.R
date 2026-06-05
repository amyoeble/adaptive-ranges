# getSD.R
# Get SD of PSE estimates from adaptive matching pilot data (Naumann, 2023)
#
# Input: Subjtest_cond15_20230609_0908.txt
# Output: -
# Packages: -
#
# Last mod: 09/06/2023, AO

dat <- read.table("Subjtest_cond15_20230609_0908.txt", header=TRUE)

f <- dat$stair  # save grouping factor before splitting dataframe
# Mark reversal points
dat <- unsplit(lapply(split(dat, dat$stair), function(data) {
    # if comparison intensity difference has a different sign in two consecutive
    #   trials there has been a reversal
    isRev <- c(diff(data$comparison) < 0, 2) == c(2, diff(data$comparison) > 0)
    # endpoint of staircase is always reversal
    isRev[length(isRev)] <- TRUE
    # ignore first three reversals (step sizes still varying)
    reversal <- isRev
    reversal[which(reversal)[1:3]] <- FALSE
    cbind(data, isRev=isRev, reversal=reversal)
}), f)
rm(f)   # clear grouping factor

# ceiling effect for snd_comparison_70 staircase, see:
dat[dat$stair == "snd_comparison_70", ]

# Get SE for PSE estimates with last 8 reversals
aggregate(comparison ~ stair, dat[dat$reversal == TRUE, ], 
    function(x) sd(x)/sqrt(length(x)))  # ~ 2
# Get SE for PSE estimates with "same" responses 
xtabs( ~ stair + response, dat)     # many "same" responses
aggregate(comparison ~ stair, dat[dat$response == "space", ], 
    function(x) sd(x)/sqrt(length(x)))   # ~ 0.8
# Conclusion: SD for "same" responses estimate is way lower, and therefore 
#   better (but much less controllable...)
