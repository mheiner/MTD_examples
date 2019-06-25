# MTD Example Code

This folder contains code to reproduce results from the associated article by Heiner and Kottas. All MCMC is run in [Julia](https://julialang.org/) v1.1+. One can work directly with the .jl scripts, or call them with the .sh scripts.

## Packages

Before running the .jl scripts, it is necessary to install the packages that appear after all instances of `using` in the .jl scripts. To do this, open the Julia REPL and enter package mode by pressing `]`. In package mode, enter `add` followed by each package name, separated by spaces. The SparseProbVec and MTD packages are unregistered and can be installed with

```julia
pkg> add "https://github.com/mheiner/SparseProbVec.jl"

pkg> add "https://github.com/mheiner/MTD.jl"
```

Post processing is run with the R scripts, which require the coda package.

## Data

The data folder contains the two simulations described in the paper (collected in sets of 1,000 observations and associated true transition probabilities), as well as the pink salmon data set (in the public domain, see <https://inport.nmfs.noaa.gov/inport/item/17256>).

## Models

The .jl scripts contain various options for prior and MCMC settings/initialization. Models can be run interactively by un-commenting appropriate lines in the first block and ignoring the `ARGS` assignments.

## Folders

- BayesFactors: Files reporting Bayes factors (from the .jl scripts) are saved in this folder. Factors are calculated for full, unrestricted Markov chains using the lag configurations listed, and are reported as the ratio of the highest marginal likelihood to the current.
- data: Contains data.
- forecast_validation: Files reporting out-of-sample validation metrics (calculated in the GoF sections of .jl scripts) are saved in this folder. The metrics are as follows.
  * L1 loss: Mean (across time points and MCMC iterations) L1 loss of one-step forecast probabilities against observations (a vector of zeros with a one at the index of the observed state).
  * Misclass rate: Misclassification rate (across time points and MCMC iterations) with respect to highest one-step forecast probabilities.
  * NLL: Mean (across time points and MCMC iterations) negative log-likelihood of one-step forecast probabilities evaluated at validation observations.
  * SqErr: Mean (across time points and MCMC iterations) squared error between observed state (as an integer) and its expected value calculated from one-step forecast probabilities. Assumes meaningful ordering of states.
  * postmean L1 loss, Ptrue: Mean (across time points and MCMC iterations) L1 loss of one-step forecast probabilities against the true probabilities.
  * L1 loss, Ptrue postmean forec: One-step forecast L1 loss calculated from posterior mean forecast probabilities against the true probabilities. (This is the metric used in the article.)
- plots: R post processing saves plots in this folder.
- postsim: Posterior simulations are saved to this folder.
- postsim_progress: Files reporting MCMC progress and statistics are saved in this folder.
