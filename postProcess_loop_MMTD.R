# rm(list=ls())
library("coda")
ids = commandArgs(trailingOnly = TRUE)

### user input

simregexp = paste0(ids[1], ".*_MMTD_.*", ids[2], ".*\\.rda")
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
# cols[floor(ncollevels/2)] = "#FFFFFF"
(levels=seq(0.0, 1.0, length=ncollevels))

for ( i in 1:nfiles ) {
  load( paste0(dctory, "/", selfiles[i]) )
  nsim = length(sims_llik)
  nuse = 2000
  whichiter = floor(seq(1, nsim, length=nuse))

  pdf(paste0("plots/trace/", mesg, ".pdf"), height=4, width=9)

  traceplot(as.mcmc(sims_llik[whichiter]), main="log-likelihood")

  for ( j in 1:(R+1) ) {
    traceplot(as.mcmc(sims_Lam[whichiter,j]), main=paste0("Lambda ", j-1))
  }
  (pmLam = colMeans(sims_Lam))
  ordR = order(pmLam, decreasing=TRUE)
  (Rselect = ifelse(which.max(pmLam) == 1, ordR[2] - 1, ordR[1] - 1))

  (pmlam = colMeans(sims_lam[[Rselect]]))
  (lamselect = which.max(pmlam))

  traceplot(as.mcmc(sims_lam[[Rselect]][whichiter,lamselect]),
            main=paste0("lambda for lags: ",
                        paste(lamindx[[Rselect]][[lamselect]], collapse=",")))

  pmQ0 = colMeans(sims_Q0)

  for(k in 1:ncol(sims_Q0)) {
  	traceplot(as.mcmc(sims_Q0[whichiter,k]), main=paste0("Q0: ", k))
  }

  (ncolQ = ncol(sims_Q[[Rselect]]))
  (whichQ = floor(seq(1, ncolQ, length=3)))
  traceplot(as.mcmc(sims_Q[[Rselect]][whichiter, whichQ[1]]), main="from Q")
  traceplot(as.mcmc(sims_Q[[Rselect]][whichiter, whichQ[2]]), main="from Q")
  traceplot(as.mcmc(sims_Q[[Rselect]][whichiter, whichQ[3]]), main="from Q")
  # traceplot(as.mcmc(sims_Q1[whichiter,15]))
  # plot(as.mcmc(sims_p1))

  # traceplot(as.mcmc(sims_p1Q[[Rselect]][whichiter,1]), main="one of the p1Q")

  dev.off()

  (pmQQ = matrix(colMeans(sims_Q[[Rselect]]), nrow=K, ncol=K^Rselect))

  pdf(paste0("plots/postmeanQ_", mesg, ".pdf"), height=1.0*K, width=1.0*K^Rselect)

  print(levelplot(matrix(pmQ0, ncol=1), at=levels, col.regions=cols,
                  main=paste0("posterior mean: Q 0")))

  froms_ls = list()
  for (r in 1:Rselect) {
    froms_ls[[r]] = 1:K
  }
  fromsQ = apply(expand.grid(froms_ls), 1, paste, collapse="")

  print(levelplot(t(pmQQ), at=levels, col.regions=cols,
                  main=paste0("posterior mean: Q ", Rselect),
                  xlab="lagged states", ylab=expression(s[t]),
                  scales=list(x=list(at=1:K^Rselect, labels=fromsQ, rot=90), cex=1.0)
                  ))

  dev.off()

  pdf(paste0("plots/postdenslambda_", mesg, ".pdf"), height=4, width=10)
  par(mfrow=c(1,R))
  for ( r in 1:(R+1) ) {
    plot(density(sims_Lam[,r]), xlab=bquote(Lambda[.(r-1)]),
         main=bquote("Order"~.(r-1)), axes=FALSE, xlim=c(0,1))
    axis(side=1)
    axis(side=2)
  }
  par(mfrow=c(2,4))
  for ( ell in 1:ncol(sims_lam[[Rselect]]) ) {
    plot(density(sims_lam[[Rselect]][,ell]), xlab=bquote(lambda[.(ell)]),
         paste0("Lags ", main=paste(lamindx[[Rselect]][[ell]], collapse=",")),
         axes=FALSE, xlim=c(0,1))
    axis(side=1)
    axis(side=2)
  }
  dev.off()

  ### lag importance across all orders

  Lamlam_sims = list()
  for(jj in 1:R) {
    Lamlam_sims[[jj]] = matrix(0.0, nrow=nsim, ncol=lamlens[jj])
    for(ii in 1:nuse) {
      Lamlam_sims[[jj]][ii,] = sims_Lam[whichiter[ii],jj+1] * sims_lam[[jj]][whichiter[ii],]
    }
  }
  (pmLamlam = lapply(Lamlam_sims, colMeans))
  sum(sapply(pmLamlam, sum)) + mean(sims_Lam[,1])

  lagweight = matrix(0.0, nrow=nuse, ncol=L+1)
  lagweight[,1] = sims_Lam[whichiter,1]
  for (ell in 1:L) {
    for (r in 1:R) {
      for (ii in 1:lamlens[r]) {
        if (ell %in% lamindx[[r]][[ii]]) {
          lagweight[,ell+1] = lagweight[,ell+1] + Lamlam_sims[[r]][,ii]
        }
      }
    }
  }
  head(lagweight)
  rowSums(lagweight[1:5,])
  pmLagweight = colMeans(lagweight)
  (hpds = HPDinterval(as.mcmc(lagweight)))
  quantile(lagweight[,2])

  setEPS()
  postscript(file=paste0("plots/PosteriorMeanLagInclusion_", mesg, ".eps"), width=L/1.5, height=3.5)
    bp = barplot(pmLagweight, ylim=c(0,1), names=0:L, ylab="Inclusion index", xlab="Lag")
    for(ell in 1:(L+1)) arrows(bp[ell,1], hpds[ell,1], bp[ell,1], hpds[ell,2], code=3, angle=90, length=0.05)
  dev.off()

  pdf(file=paste0("plots/PosteriorMeanLagInclusion_", mesg, ".pdf"), width=L*0.5, height=3.0)
  par(mar=c(4,4,1,1)+0.1)
  bp = barplot(pmLagweight, ylim=c(0,1), names=0:L, ylab="Inclusion index", xlab="Lag")
  for(ell in 1:(L+1)) arrows(bp[ell,1], hpds[ell,1], bp[ell,1], hpds[ell,2], code=3, angle=90, length=0.05)
  dev.off()

  cat("\r file", i, "of", nfiles)
}
