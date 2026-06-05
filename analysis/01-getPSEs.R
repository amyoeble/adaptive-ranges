# 01-getPSEs.R
# Preprocessing: estimate PSEs from adaptive ranges data
#
# Input: ../data/subj*.txt
#        ../data/demographics-matchRanges.txt
# Output: PSEs.txt
#         ./plots/stairs.pdf
# Packages: -
#
# Last mod: 12/29/2025, AO

# Preamble: demographics
dem <- read.table("../data/demographics-matchRanges.txt", 
        col.names=c("gender", "age"))
table(dem$gender)
mean(dem$age); sd(dem$age)

# Actual data
dat0 <- do.call(rbind, 
    lapply(dir("../data/", "Subj.*\\.txt", full.names=TRUE),
        read.table, header=TRUE, as.is=FALSE))
dat0$range <-  factor(dat0$range, levels=c("low", "medium", "high"))
dat0$subj <- factor(dat0$subj)

# check amount of training trials
aggregate(standard ~ subj + range, dat0[dat0$stair == "training", ], length)

# exclude training trials
dat <- dat0[dat0$stair != "training", ]

# Add reversal factor
dat <- do.call(rbind, lapply(split(dat, list(dat$subj, dat$stair), drop=TRUE), 
    function(data) {
        # if comparison intensity difference has a different sign in two 
        # consecutive trials there has been a reversal
        isRev <- c(diff(data$comparison) <= 0, 2) == 
                                c(2, diff(data$comparison) > 0)
        # endpoint of staircase is always reversal
        isRev[length(isRev)] <- TRUE
        # ignore first three reversals (step sizes still varying)
        reversal <- isRev
        reversal[which(reversal)[1:3]] <- FALSE
        cbind(data, isRev=isRev, reversal=reversal)
}))
rownames(dat) <- NULL

# Look at staircase trials (if interested... a LOT of plots)
pdf("plots/stairs.pdf", pointsize=10, height=5, width=6)
for(i in levels(dat$subj)) {
    par(mfrow=c(5, 4), mar=c(3, 3, 1, 1), mgp=c(2, .7, 0), oma=c(0, 0, 4, 0))
    lapply(split(dat[dat$subj == i, ], dat[dat$subj == i, ]$stair, drop=TRUE), 
        function(x) {
          plot(x$comparison ~ x$trial, type="o", pch=".", main=x$stair[1], 
              xlim=c(0,  max(dat[dat$subj == i, ]$trial)), ylim=c(18, 90), 
              xlab="Trial", ylab="Comparison (dB)")
          abline(h=mean(x[x$response == "space", ]$comparison), col="darkblue")
          abline(h=mean(x[x$reversal == TRUE, ]$comparison), col="red")
          mtext(paste0("Subj0", i), side=3, outer=TRUE, line=2)
        }  
    )
    legend("bottomright", col="red", lty=1, cex=0.8, bty="n",
        legend="mean (last 8 rev)")
}
dev.off()

# PSE estimates per stair: mean of reversals on smallest step size
stairPSEs <- aggregate(comparison ~ subj + range + stair + standard, 
    dat[dat$reversal == TRUE, ], mean)
stairPSEs <- stairPSEs[order(stairPSEs$subj, stairPSEs$range), ]
# SE of staircase PSE estimates
stairPSEs$se <- aggregate(comparison ~ subj + range + stair + standard, 
     dat[dat$reversal == TRUE, ], sd)$comparison / sqrt(8)
summary(stairPSEs$se)
# Rename PSE column (for estimation function)
colnames(stairPSEs)[which(colnames(stairPSEs) == "standard")] <- "stan"
colnames(stairPSEs)[which(colnames(stairPSEs) == "comparison")] <- "yhat"
# PSE estimates: mean of "same" responses (just for comparison)
#   NOTE: some staircases do not have any "same" responses!
aggregate(comparison ~ subj + range + stair + standard,
    dat[dat$response == "space", ], mean)

# PSE estimates
aggregate(yhat ~ subj + range + stan, stairPSEs, mean)
# SD of PSE estimates
summary(aggregate(yhat ~ subj + range + stan, stairPSEs, sd)$yhat)

# Export PSE file (one PSE per staircase)
write.table(stairPSEs, "PSEs.txt", quote=FALSE, row.names=FALSE)

# Staircase descriptives:
trialN <- aggregate(trial ~ stair + subj, dat, max)
# Mean length of staircases
mean(trialN$trial)  # 30.86 trials
# Percentage of 50 trial max reached
mean(trialN$trial == 50)    # 5 % of staircases

# Proportions of 'same' responses
proportions(aggregate(response ~ subj, dat, table)[, 2], 1)[,3] |> round(4)
