load(file.path("saved_samples", "tobamo_FA_errs_grass.Rdata"))

datBAY <- as.data.frame(t(aaply(BAY_error, .margins=c(2,3), .fun=mean)))
datMLE <- as.data.frame(t(aaply(MLE_error, .margins=c(2,3), .fun=mean)))
datBAY <- gather(datBAY, key=predictands, value=error)
datMLE <- gather(datMLE, key=predictands, value=error)
datBAY$model <- "Bayes"
datMLE$model <- "MLE"

dat <- join(datMLE, datBAY, type="full")
dat$predictands <- strtoi(dat$predictands)
levels(dat$model)

pdf("tobamovirus_errs_FA.pdf", width=6, height=4)
ggplot(data=dat, aes(factor(predictands), error)) +
  geom_boxplot(aes(fill=factor(model, levels=c("MLE", "Bayes")))) + 
  scale_fill_manual(name="", values=c("grey45", "grey90")) +
  theme_classic() + xlab("number of predictands") + ylab("mean absolute error") +
  coord_cartesian(ylim=c(.7, 1.8))
dev.off()