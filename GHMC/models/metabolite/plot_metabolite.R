library(ggplot2)
library(tidyr)
dat <- read.table(file.path("saved_samples", "metabolite_errors_per_missing.txt"))
dat2 <- read.table(file.path("saved_samples", "metabolite_errors_per_missing_7_8.txt"))
colnames(dat) <- c("miss_percent", "Stiefel", "VB", "PPCA")
colnames(dat2) <- c("miss_percent", "Stiefel", "VB", "PPCA")
dat2$PPCA <- NULL
dat$PPCA <- NULL

dat <- gather(dat, key=method, value=error, 2:3)
dat2 <- gather(dat2, key=method, value=error, 2:3)
dat <- rbind(dat, dat2)

dat <- dat[dat$miss_percent < 0.9,]
dat$miss_percent <- factor(dat$miss_percent)

pdf(file=file.path("plots", "metabolite.pdf"), width=6, height=5)
ggplot(dat, aes(x=miss_percent)) +
  geom_boxplot(aes(y=error, fill=factor(method, levels=c("VB", "Stiefel", "PPCA")))) +
  scale_fill_manual(name="", values=c("grey100","grey70","grey45")) +
  theme_classic() + ylim(0.05, .17) + 
  ylab("mean absolute error") + xlab("proportion missing")
dev.off()