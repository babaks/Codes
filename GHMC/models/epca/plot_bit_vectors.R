# Create plots for bit vector results

load(file.path("saved_samples", "bit_vec_Stiefel_results_3.Rdata"))
png(file.path("plots", "bit_images.png"), width=1400, height=470)
par(mfrow=c(1,3), mar=rep(2,4), cex=1.5)
image(t(true_bits[seq(1, 600, by=2),]),axes = FALSE, col = grey(seq(0, 1, length = 256)), main="Original")
Y = (2 * Y) %% 3
image(t(Y[seq(1, 600, by=2),]),axes = FALSE, col = grey(seq(0, 1, length = 256)), main="Corrupted + missing")
image(t(res$reconstruction[seq(1, 600, by=2),]),axes = FALSE, col = grey(seq(0, 1, length = 256)), main="Reconstruction")
dev.off()
print(res$mean_err)

pdf(file.path("plots", "bit_results.pdf"), width=6.5, height=5)
par(mar=2 * c(2, 2, 1.5, 1))
plot(x=1:length(vs_true), vs_true, type='l', lwd=1.2, col='grey1', log='x', 
     xlim=c(1,length(vs_true)), ylim=c(.15, .58), xlab="HMC sample #", ylab="reconstruction error")
lines(vs_corrupt, lty=1, lwd=1.2, col='grey60')

materr <- readMat(file.path("data", "bit_vector_mat_errs.mat"))
m_true <- as.vector(materr$err.true)
m_corr <- as.vector(materr$err.corr)
lines(x=1:length(m_corr), m_true, lty=4, col='grey1')
lines(x=1:length(m_corr), m_corr, lty=4, col='grey60')
legend(x='topright', legend=c("error vs original", "vs corrupted"), 
       lty=c(4,4), lwd=c(1.2,1.5), col=c('grey1', 'grey60'), title="Unconstrained par.",
       title.adj = .2)
legend(x='bottomleft', legend=c("error vs original", "vs corrupted"), 
       lty=c(1,1), lwd=c(1.2,1.2), col=c('grey1', 'grey60'), title="Stiefel parameterization",
       title.adj = .3)
dev.off()

pdf(file.path("plots", "bit_dists.pdf"), width=6.5, height=5)
plot(x=1:10000, y=dists[1:10000], 'l', xlab="HMC sample #", ylab="chordal distance to mean")
dev.off()

