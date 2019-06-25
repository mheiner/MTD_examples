# const priλ = "Dir"
# const priλ = "SDM"
# const priλ = "SBM"
#
# const priQ = "Dir"
# const priQ = "SDM"
# const priQ = "SBM"
#
# const srnd = 13119
# const TT = 50
# const L = 5
#
# whichdat = "sim1_sets"
# whichdat = "sim2_sets"
# whichdat = "pinkSalmon_first30_K4"

### if ARGS used

const priλ = ARGS[1]
const priQ = ARGS[2]
const srnd = parse(Int64, ARGS[3])
const TT = parse(Int64, ARGS[4])
const L = parse(Int64, ARGS[5])
const whichdat = ARGS[6]

### after user selections:

using BSON
using StatsBase

# pkg> add "https://github.com/mheiner/SparseProbVec.jl"
using SparseProbVec

# pkg> add "https://github.com/mheiner/MTD.jl"
using MTD

using Random
using Dates

# srnd = 10091
Random.seed!(srnd)

### Read in data
simdata = occursin(r"sim\d+_sets", whichdat) # run validation at end?

if simdata
    BSON.@load "data/$(whichdat).bson" y X P
    const S = deepcopy(y[1])
else
    BSON.@load "data/$(whichdat).bson" S
end


### Setup
const SS = S[1:TT]
println(sort(unique(SS[(L+1):TT])))
println(counts(SS))
const K = maximum(S)

const n_keep = Int(2e3)
const n_burn = Int(200e3)
const thin = 200

const mesg = "$(whichdat)_MTDg_T$(TT)_L$(L)_K$(K)_$(priλ)lam_$(priQ)Qmarg_$(srnd)"


### Prior

fieldnames(PriorMTDg)
symm_prior = symmetricDirPrior_mtdg(1.0, 1.0, 1.0, L, K)

# β_λ = log(float(TT))
β_λ = sqrt(float(TT))
# β_λ = float(TT) / 4.0
# β_λ = 1.001

η_λ = 1.0e3
π_small_λ = vcat(0.0, fill(0.50, L-1))
π_large_λ = vcat(0.0, fill(0.20, L-1))
γ_λ, δ_λ = shape_Dir2GenDir(symm_prior[1])
δc_λ = del_correctionSBM(π_small_λ[2], π_large_λ[2], L)
γ_λ0 = 1.0
δ_λ0 = 1.0


# β_Q = float(TT) / 4.0 # for SDM on Q
β_Q = float(K)

η_Q = 1.0e3
π_small_Q = 1.0 / float(K)
π_large_Q = π_small_Q / 5.0
γ_Q, δ_Q = shape_Dir2GenDir(symm_prior[2][:,1])
δc_Q = del_correctionSBM(π_small_Q, π_large_Q, K)



if priλ == "Dir"
        prior_λ = deepcopy(symm_prior[1])
    elseif priλ == "SDM"
        prior_λ = SparseDirMix(symm_prior[1], β_λ)
    elseif priλ == "SBM"
        prior_λ = SBMprior(L+1, η_λ, π_small_λ, π_large_λ, vcat(γ_λ0, γ_λ), vcat(δ_λ0, δc_λ.*δ_λ))
end

prior_Q0 = deepcopy(symm_prior[2])

if priQ == "Dir"
        prior_Q = deepcopy(symm_prior[3])
    elseif priQ == "SDM"
        prior_Q = [ [ SparseDirMix(symm_prior[3][j][:,kk], β_Q) for kk = 1:K ] for j = 1:L ]
    elseif priQ == "SBM"
        prior_Q = [ [ SBMprior(K, η_Q, π_small_Q, π_large_Q, γ_Q, δc_Q*δ_Q) for kk = 1:K ] for j = 1:L ]
end

# bfact = bfact_MC(SS, L, 1, K, prior_Q0, prior_Q)
#
# using Printf
# bfact_file = open("BayesFactors/bfact_$(mesg).txt", "a+")
# write(bfact_file, "order 0, lags 0: $(@sprintf("%.2e", bfact[2][1,3])), llik: $(@sprintf("%.2e", bfact[2][1,4])) \n")
# for i in 1:bfact[1].nZζ
#     Znow, ζnow = bfact[1].Zζindx[i,:]
#     write(bfact_file, "order $(Znow), lags $(bfact[1].indxs[Znow][ζnow]): $(@sprintf("%.2e", bfact[2][i+1,3])), llik: $(@sprintf("%.2e", bfact[2][i+1,4])) \n")
# end
# close(bfact_file)


### Inits
fieldnames(ParamsMTDg)
inits = ParamsMTDg(
    SparseProbVec.rand(prior_λ, logout=true), # λ
    [ StatsBase.sample(0:L) for i = 1:(TT-L) ], # ζ
    SparseProbVec.rDirichlet(prior_Q0, logout=true),
    [ reshape(vcat([ SparseProbVec.rDirichlet(ones(K), logout=true) for k in 1:K ]...), fill(K, 2)... ) for j = 1:L ]
  # log( transTens_MLE( count_trans_L(SS, K, 1) + 1.0 ) ) # lQ (mle)
  )

### Create model
fieldnames(ModMTDg)
model = ModMTDg(L, K, TT, SS,
  PriorMTDg(prior_λ, prior_Q0, prior_Q),
  inits)

logfilename = "postsim_progress/out_prog_$(mesg).txt"
timestart = Dates.now()

timemod!(5, model, 100, logfilename)

## Burn in
mcmc!(model, Int(0.15*n_burn)-500, save=false, report_filename=logfilename,
      thin=1, jmpstart_iter=5, report_freq=1000)
mcmc!(model, Int(0.85*n_burn), save=false, report_filename=logfilename,
      thin=1, jmpstart_iter=10, report_freq=10000)

## Estimate time remaining
etr(timestart, n_keep, thin, logfilename)

### MCMC
sims = mcmc!(model, n_keep, save=true, report_filename=logfilename,
             thin=thin, jmpstart_iter=10)

### Output
bson("./postsim/mcmc_$(mesg).bson", sims=deepcopy(sims),
    model=deepcopy(model), whichdat=deepcopy(whichdat),
    mesg=deepcopy(mesg))

### Send simulations to R
using RCall

nsim = length(sims)

R"rm(list=ls())"
priorinfo = "$(model.prior)"
sims_llik = [ sims[i][:llik] for i=1:nsim ]

sims_lam = permutedims(hcat([ exp.(sims[i][:lλ]) for i=1:nsim ]...))
sims_Q0 = permutedims(hcat([ exp.(sims[i][:lQ0]) for i=1:nsim ]...))
sims_Q = [ exp(sims[i][:lQ][r][j,jj]) for i=1:nsim, r=1:model.L, j=1:model.K, jj=1:model.K ]

@rput sims_lam sims_Q0 sims_Q sims_llik priorinfo K L TT mesg

R"ls()"
R"save.image(file=paste0('postsim/mcmc_', $(mesg), '.rda'))"


### GoF
if simdata

    using Random
    using Dates
    using Statistics

    modeltype = "MTDg"
    report_filename = "forecast_validation/gof_$(mesg).txt"
    test_indx = 2

    report_file = open(report_filename, "a+")

    nprime = length(y[test_indx])
    niter = length(sims)
    nsim = 2000

    Random.seed!(srnd)
    simind = sort(StatsBase.sample(1:niter, nsim, replace=false))

    (FL1, MC, NLL, SQE, PL1, FOREC) = meanForecLoss(y[test_indx], X[test_indx][:,1:model.L],
        P[test_indx], nprime, lossL1, sims,
        model.TT, model.L, K,
        simind, modeltype=modeltype)

    PMforec = reshape(mean(FOREC, dims=1), (nprime, K) )
    PML1P = mean(abs.(PMforec .- P[test_indx])) # lower because of Jensen's inequality

    write(report_file, "L1 loss: $(round(100.0*mean(FL1), digits=2))\n" )
    write(report_file, "Misclass rate: $(round(100.0*mean(MC), digits=2))\n" )
    write(report_file, "NLL: $(round(mean(sum(NLL, dims=2)), digits=2))\n" )
    write(report_file, "SqErr: $(round(100.0*mean(SQE), digits=2))\n" )
    write(report_file, "postmean L1 loss, Ptrue: $(round(100.0*mean(PL1), digits=2))\n" )
    write(report_file, "L1 loss, Ptrue postmean forec: $(round(100.0*PML1P, digits=2))\n" )

    close(report_file)
end
