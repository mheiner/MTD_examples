# rm(list=ls())
library("coda")
ids = commandArgs(trailingOnly = TRUE)

### user input

simregexp = paste0(ids[1], ".*_MTD_.*", ids[2], ".*\\.rda")
dctory = "./postsim"

### end user input

allfiles = list.files(dctory)
(selfiles = allfiles[grep(simregexp, allfiles)])
(nfiles = length(selfiles))

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

  traceplot(as.mcmc(sims_lam[whichiter,1]))
  traceplot(as.mcmc(sims_lam[whichiter,2]))
  traceplot(as.mcmc(sims_lam[whichiter,3]))

  traceplot(as.mcmc(sims_Q[whichiter,1]))
  traceplot(as.mcmc(sims_Q[whichiter,2]))
  traceplot(as.mcmc(sims_Q[whichiter,3]))

  dev.off()

  (pmQQ = matrix(colMeans(sims_Q), nrow=K))

  pdf(paste0("plots/postmeanQ_", mesg, ".pdf"), height=5, width=6)
  print(levelplot(t(pmQQ), at=levels, col.regions=cols, main="posterior mean: Q"))
  dev.off()

  pdf(paste0("plots/postdenslambda_", mesg, ".pdf"), height=4, width=10)
  par(mfrow=c(2,4))
  for ( ell in 1:ncol(sims_lam) ) {
    plot(density(sims_lam[,ell]), xlab=bquote(lambda[.(ell)]),
         main=bquote("Lag"~.(ell)), axes=FALSE, xlim=c(0,1))
    axis(side=1)
    axis(side=2)
  }
  dev.off()

  cat("\r file", i, "of", nfiles)
}
