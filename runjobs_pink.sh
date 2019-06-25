nice julia fit_MTD.jl SDM Dir 22619 30 5 pinkSalmon_first30_K4 &

nice julia fit_MTDg.jl SBM Dir 22619 30 5 pinkSalmon_first30_K4 &

nice julia fit_MMTD.jl SBM Dir Dir 22619 30 5 2 pinkSalmon_first30_K4 &
nice julia fit_MMTD.jl SBM SDM Dir 22619 30 5 2 pinkSalmon_first30_K4 &

wait

Rscript --vanilla postProcess_loop_MTD.R pink 22619 &
Rscript --vanilla postProcess_loop_MTDg.R pink 22619 &
Rscript --vanilla postProcess_loop_MMTD.R pink 22619 &

wait
