nice julia fit_MTDg.jl SBM Dir 22219 100 7 sim2_sets &
nice julia fit_MTDg.jl SBM Dir 22219 200 7 sim2_sets &
nice julia fit_MTDg.jl SBM Dir 22219 200 6 sim1_sets &

# wait

# nice julia fit_MTDg.jl SBM Dir 22219 500 6 sim1_sets &
# nice julia fit_MTDg.jl SBM Dir 22219 500 7 sim2_sets &

wait

Rscript --vanilla postProcess_loop_MTDg.R sim 22219 &

wait
