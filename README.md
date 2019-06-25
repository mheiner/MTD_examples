# MTD Example Code

This folder contains code to reproduce results from the associated article by Heiner and Kottas. All MCMC is run in [Julia](https://julialang.org/) v1.1+. One can work directly with the .jl scripts, or call them with the .sh scripts.

## Packages

Before running the .jl scripts, it is necessary to install the packages that appear after all instances of `using` in the .jl scripts. To do this, open the Julia REPL and enter package mode by pressing `]`. In package mode, enter `add` followed by each package name, separated by spaces. The SparseProbVec and MTD packages are unregistered and can be installed with

```julia
pkg> add "https://github.com/mheiner/SparseProbVec.jl.git"

pkg> add "https://github.com/mheiner/MTD.jl.git"
```

Post processing is run with the R scripts, which require the coda package.

## Data

The data folder contains the two simulations described in the paper (collected in sets of 1,000 observations and associated true transition probabilities), as well as the pink salmon data set (in the public domain, see <https://inport.nmfs.noaa.gov/inport/item/17256>).

## Models

The .jl scripts contain various options for prior and MCMC settings/initialization. Models can be run interactively by un-commenting appropriate lines in the first block and ignoring the `ARGS` assignments.

## Folders

- BayesFactors: Files reporting Bayes factors (from the .jl scripts) are saved in this folder.
- data: Contains data.
- forecast_validation: Files reporting out-of-sample validation metrics (calculated in the GoF sections of .jl scripts) are saved in this folder.
- plots: R post processing saves plots in this folder.
- postsim: Posterior simulations are saved to this folder.
- postsim_progress: Files reporting MCMC progress and statistics are saved in this folder.
