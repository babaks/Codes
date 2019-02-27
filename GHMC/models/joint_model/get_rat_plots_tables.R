load(file.path("saved_samples", "thinned_joint_samples_seq.Rdata"))

betas <- laply(joint_samples, .fun=function(x) x$beta)

#png
pdf(file.path("plots", "trace_plots_rat.pdf"), width=4, height=3)
par(mar=c(3,3,2,2),oma= c(0,0,0,0), mgp=c(1.5,.5,0))
plot(betas[,5],type='l',ylim=c(-10,10 ),
     xlab="Iteration",
     ylab="Trace values")
lines(betas[,4],lty=2,col="grey50")
lines(betas[,3],lty=2,col="grey50")
lines(betas[,2],lty=2,col="grey50")
lines(betas[,1],lty=2,col="grey50")
legend(x='topright',legend=c('Factor 1','Others'),
       col=c('black','grey50'),
       lty=c(1,2),
       bg="white")
dev.off()
