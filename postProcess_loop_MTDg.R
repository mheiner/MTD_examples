# rm(list=ls())
library("coda")
ids = commandArgs(trailingOnly = TRUE)

### user input

simregexp = paste0(ids[1], ".*_MTDg_.*", ids[2], ".*\\.rda")
dctory = "./postsim"

### end user input

library("lattice")
library("grDevices")
colscale = colorRampPalette(c("white", "navyblue"), bias=1)
ncollevels = 32
cols = colscale(ncollevels)
(levels=seq(0.0, 1.0, length=ncollevels))

for ( i in 1:nfiles ) {
  load( paste0("./postsim/", selfiles[i]) )
  (nsim = nrow(sims_lam))
  whichiter = floor(seq(1, nsim, length=2000))

  pdf(paste0("plots/trace/", mesg, ".pdf"), height=4, width=9)

  traceplot(as.mcmc(sims_llik[whichiter]), main="log-likelihood")

  for ( j in 0:L ) {
    traceplot(as.mcmc(sims_lam[whichiter,j+1]), ylim=c(0,1))
    title(main=paste0("lambda ", j))
  }

  for ( k in 1:(K) ) {
    traceplot(as.mcmc(sims_Q0[whichiter,k]), ylim=c(0,1))
    title(main=paste0("Q0 ", k))
  }

  traceplot(as.mcmc(sims_Q[whichiter,1,2,1]))
  traceplot(as.mcmc(sims_Q[whichiter,2,1,1]))
  traceplot(as.mcmc(sims_Q[whichiter,1,2,2]))
  traceplot(as.mcmc(sims_Q[whichiter,2,1,2]))

  dev.off()

  (pmlam = colMeans(sims_lam))
  (pmQ0 = colMeans(sims_Q0))
  (pmQQ = apply(sims_Q, 2:4, mean))
  for (j in 1:L) {
    print( pmQQ[j,,] )
  }

  pdf(paste0("plots/postmeanQ_", mesg, ".pdf"), height=5, width=6)
  for (j in 1:L) {
    print(levelplot(t(pmQQ[j,,]), at=levels, col.regions=cols, main=paste0("posterior mean: Q ",j)))
  }
  dev.off()

  pdf(paste0("plots/postdenslambda_", mesg, ".pdf"), height=4, width=10)
  par(mfrow=c(2,4))
  for ( ell in 0:L ) {
    plot(density(sims_lam[,ell+1]), xlab=bquote(lambda[.(ell)]),
         main=bquote("Lag"~.(ell)), axes=FALSE, xlim=c(0,1))
    axis(side=1)
    axis(side=2)
  }
  dev.off()

  cat("\r file", i, "of", nfiles)
}
