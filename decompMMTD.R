rm(list=ls())

library("JuliaCall")
julia = julia_setup()
julia_install_package_if_needed("https://github.com/mheiner/MTD.jl")
julia_installed_package("MTD")
julia_library("MTD")
source("decomp_functions.R")


Omeg = array( c(0.5, 0.5, 0.4, 0.6, 0.9, 0.1, 0.2, 0.8), dim=c(2,2,2) )

x = decompMMTD(Omeg)
# x = decompMMTD(Omeg, random=TRUE)

# decomp_to_julia(x)

print(x)
reconstruct(x)




