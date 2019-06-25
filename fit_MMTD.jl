# const priΛ = "Dir"
# const priΛ = "SDM"
# const priΛ = "SBM"

# const priλ = "Dir"
# const priλ = "SDM"
# const priλ = "SBM"
# const priλ = "SBMSDM" # SBM on single lag, SDM on higher orders

# const priQ = "Dir"
# const priQ = "SDM"
# const priQ = "SBM"

# const srnd = 22219
# const TT = 50
# const L = 3
# const R = 2
#
# whichdat = "sim1_sets"
# whichdat = "sim2_sets"
# whichdat = "pinkSalmon_first30_K4"

### if ARGS used

const priΛ = ARGS[1]
const priλ = ARGS[2]
const priQ = ARGS[3]
const srnd = parse(Int64, ARGS[4])
const TT = parse(Int64, ARGS[5])
const L = parse(Int64, ARGS[6])
const R = parse(Int64, ARGS[7])
const whichdat = ARGS[8]

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
const λ_indx = build_λ_indx(L, R)

const n_keep = Int(2e3)
const n_burn = Int(200e3)
const thin = 200

const mesg = "$(whichdat)_MMTD_T$(TT)_L$(L)_R$(R)_K$(K)_$(priΛ)Lam_$(priλ)lam_$(priQ)Qmarg_$(srnd)"


### Prior

fieldnames(PriorMMTD)
symm_prior = symmetricDirPrior_mmtd(1.0, 1.0, 1.0, 1.0, L, R, K, λ_indx)

η_Λ = 1.0e3
π_small_Λ = vcat(0.0, fill(0.25, R-1))
# π_small_Λ = vcat(0.0, fill(1.0/float(R+1), R-1))
π_large_Λ = π_small_Λ .* 1.0
α_Λ = [ 5.0 / 2.0^float(r-1) for r = 1:(R+1) ]
# α_Λ = vcat(float(R), fill(1.0/float(R), R))
# γ_Λ, δ_Λ = shape_Dir2GenDir(α_Λ)
γ_Λ = 1.0
δ_Λ = 1.0
# δc_Λ = del_correctionSBM(π_small_Λ[2], π_large_Λ[2], L)
δc_Λ = 1.0

# β_λ = log(float(TT))
β_λ = sqrt(float(TT))
# β_λ = float(TT) / 10.0
# β_λ = 1.001

η_λ = 1.0e3
π_small_λ = 0.75
π_large_λ = 0.10
γ_λ, δ_λ = shape_Dir2GenDir(symm_prior[2][1]) # this works for SBMSDM only (will need one for each set otherwise)
δc_λ = del_correctionSBM(π_small_λ, π_large_λ, L)
# γ_λ = 1.5
# δ_λ = 1.5
# δc_λ = 1.0

# β_Q = float(TT) / 8.0 # for SDM on Q
β_Q = float(K) / 2.0

η_Q = 1.0e3
π_small_Q = 1.0 / float(K)
π_large_Q = π_small_Q / 5.0
γ_Q, δ_Q = shape_Dir2GenDir(symm_prior[4][1][:,1])
δc_Q = del_correctionSBM(π_small_Q, π_large_Q, K)


if priΛ == "Dir"
        prior_Λ = deepcopy(symm_prior[1]) # Dirichlet
    elseif priΛ == "DirLowM"
        prior_Λ = deepcopy(α_Λ) # Dirichlet with decreasing α
    elseif priΛ == "SDM"
        prior_Λ = SparseDirMix(symm_prior[1], β_Λ)
    elseif priΛ == "SBM"
        prior_Λ = SBMprior(R+1, η_Λ, π_small_Λ, π_large_Λ, γ_Λ, δc_Λ.*δ_Λ)
end

if priλ == "Dir"
        prior_λ = deepcopy(symm_prior[2])
    elseif priλ == "SDM"
        prior_λ = [ SparseDirMix(symm_prior[2][r], β_λ) for r in 1:R ]
    elseif priλ == "SBM"
        prior_λ = [ SBMprior(λ_indx.lens[r], η_λ, π_small_λ, π_large_λ, γ_λ, δc_λ.*δ_λ) for r in 1:R ]
    elseif priλ == "SBMSDM"
        prior_λ = Vector{Union{SparseDirMix, SBMprior}}(undef, R)
        prior_λ[1] = SBMprior(λ_indx.lens[1], η_λ, π_small_λ, π_large_λ, γ_λ, δc_λ.*δ_λ)
        for r in 2:R
            prior_λ[r] = SparseDirMix(symm_prior[2][r], β_λ)
        end
end

prior_Q0 = deepcopy(symm_prior[3])

if priQ == "Dir"
        prior_Q = deepcopy(symm_prior[4]) # for Dirichlet
    elseif priQ == "SDM"
        prior_Q = [ [ SparseDirMix(symm_prior[4][1][:,1], β_Q) for kk in 1:(K^r) ] for r in 1:R ] # SDM
        prior_Q = [ reshape(prior_Q[r], fill(K, r)... ) for r in 1:R ]
    elseif priQ == "SBM"
        prior_Q = [ [ SBMprior(K, η_Q, π_small_Q, π_large_Q, γ_Q, δc_Q.*δ_Q) for kk in 1:(K^r) ] for r in 1:R ] # SBM
        prior_Q = [ reshape(prior_Q[r], fill(K, r)... ) for r in 1:R ]
end

bfact = bfact_MC(SS, L, R, K, prior_Q0, prior_Q)

using Printf
bfact_file = open("BayesFactors/bfact_$(mesg).txt", "a+")
write(bfact_file, "order 0, lags 0: $(@sprintf("%.2e", bfact[2][1,3])), llik: $(@sprintf("%.2e", bfact[2][1,4])) \n")
for i in 1:bfact[1].nZζ
    Znow, ζnow = bfact[1].Zζindx[i,:]
    write(bfact_file, "order $(Znow), lags $(bfact[1].indxs[Znow][ζnow]): $(@sprintf("%.2e", bfact[2][i+1,3])), llik: $(@sprintf("%.2e", bfact[2][i+1,4])) \n")
end
close(bfact_file)


### Inits
fieldnames(ParamsMMTD)
inits = ParamsMMTD(
  SparseProbVec.rDirichlet(α_Λ, logout=true), # Λ
  # log(Λ), # Λ
  [ SparseProbVec.rDirichlet(symm_prior[2][r], logout=true) for r in 1:R], # λ
  # [ log(λ[r]) for r in 1:R ], # λ
  [ StatsBase.sample(1:(λ_indx.nZζ+1)) - 1 for i in 1:(TT-L) ], # Zζ
  # ζ[1:(TT-L),:], # ζ
  SparseProbVec.rDirichlet(prior_Q0, logout=true),
  # [ log(fill(1.0/K, fill(K, r+1)...)) for r in 1:R ] # lQ (const)
  [ reshape(vcat([ SparseProbVec.rDirichlet(ones(K), logout=true) for k in 1:K^(r) ]...), fill(K, r+1)... ) for r in 1:R ] # lQ (rand)
  # [ log( transTens_MLE( count_trans_L(SS, K, r) + 1.0 ) ) for r in 1:R ] # lQ (mle)
  # [ log(Q[r]) for r in 1:R ] # lQ
  )

### Create model
fieldnames(ModMMTD)
model = ModMMTD(L, R, K, TT, SS,
  PriorMMTD(prior_Λ, prior_λ, prior_Q0, prior_Q),
  inits, λ_indx)

logfilename = "postsim_progress/out_prog_$(mesg).txt"
timestart = Dates.now()

# timing for benchmarks
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

lamindx = model.λ_indx.indxs
lamlens = model.λ_indx.lens
lamZzetaindx = model.λ_indx.Zζindx
lamNZzeta = model.λ_indx.nZζ

sims_Lam = permutedims(hcat([ exp.(sims[i][:lΛ]) for i=1:nsim ]...))
sims_lam = [ permutedims(hcat([ exp.(sims[i][:lλ][r]) for i=1:nsim]...)) for r in 1:R ]
sims_Q0 = permutedims(hcat([ exp.(sims[i][:lQ0]) for i=1:nsim ]...))
sims_Q = [ permutedims(hcat([ exp.(vec(sims[i][:lQ][r])) for i=1:nsim]...)) for r in 1:R ]

@rput sims_Lam sims_lam sims_Q0 sims_Q sims_llik priorinfo K L R TT lamindx lamlens lamZzetaindx lamNZzeta mesg

R"ls()"
R"save.image(file=paste0('postsim/mcmc_', $(mesg), '.rda'))"


### GoF
if simdata

    using Random
    using Dates
    using Statistics

    modeltype = "MMTD"
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
        simind, λ_indx=model.λ_indx, modeltype=modeltype)

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
