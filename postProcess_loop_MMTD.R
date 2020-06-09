library("coda")
ids = commandArgs(trailingOnly = TRUE)

### user input

# simregexp = "pink.*22619.rda"
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

  for ( j in 1:(R+1) ) {
    traceplot(as.mcmc(sims_Dlevelwgt[whichiter,j]), main=paste0("Decomposition level weight ", j-1))
  }
  for ( j in 1:L ) {
    traceplot(as.mcmc(sims_Dlaginvolv[whichiter,j]), main=paste0("Decomposition lag involvement ", j))
  }

  # pmQ0 = colMeans(sims_Q0)
  pmQ0 = apply(sims_Q0, 2, median)

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


  pdf(paste0("plots/postmedianQ0_", mesg, ".pdf"), height=0.75*K, width=0.8*K)
  print(
    levelplot(matrix(pmQ0, ncol=1), at=levels, col.regions=cols,
              main=paste0("posterior median: Q 0"))
  )
  dev.off()

  for (rr in 1:R) {
    # (pmQQ = matrix(colMeans(sims_Q[[rr]]), nrow=K, ncol=K^rr))
    (pmQQ = matrix(apply(sims_Q[[rr]], 2, median), nrow=K, ncol=K^rr))

    froms_ls = list()
    for (i in 1:rr) {
      froms_ls[[i]] = 1:K
    }
    fromsQ = apply(expand.grid(froms_ls), 1, paste, collapse="")

    # pdf(paste0("plots/postmeanQ", rr, "_", mesg, ".pdf"), height=1.5*K, width=1.25*K^rr)
    pdf(paste0("plots/postmedianQ", rr, "_", mesg, ".pdf"), height=.75*K, width=.75*K^rr) # for pink salmon plots

    print(
      levelplot(t(pmQQ), at=levels, col.regions=cols,
                main=paste0("posterior median: Q ", rr),
                xlab="lagged states", ylab=expression(s[t]),
                scales=list(x=list(at=1:K^rr, labels=fromsQ, rot=90), cex=1.0))
    )
    dev.off()
  }

  pdf(paste0("plots/postdensLLambda_", mesg, ".pdf"), height=3, width=3.25)
  par(mar=c(4,4,3,1)+0.1)
  d1 = list()
  d2 = list()
  for ( j in 1:(R+1) ) {
    d1[[j]] = density(sims_Dlevelwgt[,j])
    d2[[j]] = density(sims_Lam[,j])
  }
  ylim = c(0, max(sapply(d1, function(z) max(z$y)), sapply(d2, function(z) max(z$y))))
  for (j in 1:(R+1)) {
    plot(d1[[j]], xlab=bquote(Lambda[.(j-1)]), ylim=ylim,
         main=bquote("Level"~.(j-1)), axes=FALSE, xlim=c(0,1), lty=2)
    lines(d2[[j]])#, xlab=bquote(Lambda[.(j-1)]),
    #main=bquote("Level"~.(j-1)), axes=FALSE, xlim=c(0,1))
    axis(side=1)
    axis(side=2)
  }
  for (j in min(R+2,L+1):(L+1)) {
    d1[[j]] = density(sims_Dlevelwgt[,j])
    plot(d1[[j]], xlab=paste("level weight: ", j-1),
         main=bquote("Level"~.(j-1)), axes=FALSE, xlim=c(0,1), lty=2)
    axis(side=1)
    axis(side=2)
  }
  dev.off()

  pdf(paste0("plots/postdenslambda_", mesg, ".pdf"), height=4, width=10)
  par(mfrow=c(1,R+1))
  for ( j in 1:(R+1) ) {
    d1 = density(sims_Dlevelwgt[,j])
    d2 = density(sims_Lam[,j])
    ylim = c(0, max(d1$y, d2$y))
    plot(d1, xlab=bquote(Lambda[.(j-1)]), ylim=ylim,
         main=bquote("Level"~.(j-1)), axes=FALSE, xlim=c(0,1), lty=2)
    lines(d2)#, xlab=bquote(Lambda[.(j-1)]),
         #main=bquote("Level"~.(j-1)), axes=FALSE, xlim=c(0,1))
    axis(side=1)
    axis(side=2)
  }

  for (r in 1:R) {
    whichsim_cutoff = median(sims_Lam[,r+1])
    whichsim = which(sims_Lam[,r+1] > whichsim_cutoff)
    par(mfrow=c(2,4))
    for ( ell in 1:ncol(sims_lam[[r]]) ) {
      plot(density(sims_lam[[r]][whichsim,ell]), xlab=bquote(lambda[.(ell)]),
           paste0("Lags ", main=paste(lamindx[[r]][[ell]], collapse=",")),
           axes=FALSE, xlim=c(0,1))
      axis(side=1)
      axis(side=2)
    }
  }
  dev.off()

  ### lag importance across all levels

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
  for (r in 1:L) {
    for (jj in 1:R) {
      for (ii in 1:lamlens[jj]) {
        if (r %in% lamindx[[jj]][[ii]]) {
          lagweight[,r+1] = lagweight[,r+1] + Lamlam_sims[[jj]][whichiter,ii] # + pmLamlam[[j]][i]
        }
      }
    }
    #cat(r, "\n")
  }
  head(lagweight)
  rowSums(lagweight[1:5,])
  colnames(lagweight) = paste(0:L)

  # pairs(lagweight, xlim=c(0,1), ylim=c(0,1), pch=".")

  Dlaginvolv = cbind(sims_Dlevelwgt[,1], sims_Dlaginvolv)
  colnames(Dlaginvolv) = paste(0:L)

  # setEPS()
  # postscript(file=paste0("plots/PosteriorMeanLagInclusion_", mesg, ".eps"), width=L/1.5, height=3.5)
  #   bp = barplot(pmLagweight, ylim=c(0,1), names=0:L, ylab="Inclusion index", xlab="Lag")
  #   for(r in 1:(L+1)) arrows(bp[r,1], hpds[r,1], bp[r,1], hpds[r,2], code=3, angle=90, length=0.05)
  # dev.off()

  library("tidyverse")
  library("ggplot2")

  plot_laginvolv = function(X) {
    ncols = ncol(X)
    colnames(X) = paste(0:(ncols-1))
    df = as.data.frame(X)
    dfg = gather(df)
    colnames(dfg) = c("lag", "weight")
    dfg$lag = factor(dfg$lag, levels=colnames(X))
    ggplot(dfg, aes(x=lag, y=weight)) + xlab("Lag") + ylab("Total weight") + ylim(c(0,1)) +
      geom_violin(color="lightblue3", fill="lightblue3", scale="width") + #, draw_quantiles=c(0.5)) + # draw_quantiles based on kde
      stat_summary(fun="median", geom="point", size=3) +
      # stat_summary(mapping = aes(x=lag, y=weight), fun.min=function(z){quantile(z, 0.025)}, fun.max=function(z){quantile(z, 0.975)}, fun=median) +
      stat_summary(mapping = aes(x=lag, y=weight), fun.min=function(z){HPDinterval(as.mcmc(z))[1,"lower"]}, fun.max=function(z){HPDinterval(as.mcmc(z))[1,"upper"]}, fun=median) +
      theme_bw()
  }

  # pdf(file=paste0("plots/postLagInclusion_", mesg, ".pdf"), width=max(.5*L, sqrt(L)+1), height=2.5)
  pdf(file=paste0("plots/postLagInclusion_", mesg, ".pdf"), width=max(.25*L, sqrt(.65*L)+1), height=1.8)

  par(mar=c(4,4,1,1)+0.1)
  p = plot_laginvolv(lagweight)
  print(p)
  # bp = barplot(pmLagweight, ylim=c(0,1), names=0:L, ylab="Inclusion index", xlab="Lag")
  # for(r in 1:(L+1)) arrows(bp[r,1], hpds[r,1], bp[r,1], hpds[r,2], code=3, angle=90, length=0.05)

  par(mar=c(4,4,1,1)+0.1)
  p = plot_laginvolv(Dlaginvolv)
  print(p)
  # bp = barplot(pmDlaginvolv, ylim=c(0,1), names=0:L, ylab="Inclusion index", xlab="Lag")
  # for(r in 1:(L+1)) arrows(bp[r,1], hpdD[r,1], bp[r,1], hpdD[r,2], code=3, angle=90, length=0.05)
  dev.off()

  colnames(sims_Lam) = paste(0:R)
  colnames(sims_Dlevelwgt) = paste(0:L)

  pdf(file=paste0("plots/postLagPairs_", mesg, ".pdf"), width=L, height=L)
    pairs(sims_Lam, xlim=c(0,1), ylim=c(0,1), pch=".", lower.panel=NULL, main="Model level pairs")
    pairs(sims_Dlevelwgt[,1:(R+1)], xlim=c(0,1), ylim=c(0,1), pch=".", lower.panel=NULL, main="Decomp. level pairs")

    pairs(lagweight, xlim=c(0,1), ylim=c(0,1), pch=".", lower.panel=NULL, main="Model lag pairs")
    pairs(Dlaginvolv, xlim=c(0,1), ylim=c(0,1), pch=".", lower.panel=NULL, main="Decomp. lag pairs")
  dev.off()

  cat("\r file", i, " of", nfiles)
}
