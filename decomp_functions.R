# library("JuliaCall")
# julia = julia_setup()
# julia_library("MTD")


decompMMTD = function(Omega, print=TRUE, totalweight_thresh=0.001, maxiter=100, leftoversum_tol=1.0e-9, random=FALSE) {
  julia_assign("Omeg", Omega)
  julia_command(paste0("r = reduce_transTens(TransTens(Omeg),",
                       "totalweight_thresh=", totalweight_thresh, ",",
                       "maxiter=", maxiter, ",",
                       "leftoversum_tol=", leftoversum_tol, ",",
                       "random=", ifelse(random, "true", "false"),
                       ")"), 
                show_value=FALSE)
  
  out = list()
  out$Omega = Omega
  out$order = julia_eval("r.order")
  out$totalweight_tol = julia_eval("r.totalweight_tol")
  out$levels = julia_eval("r.levels")
  out$configs = as.list(julia_eval("r.configs"))
  out$Q = as.list(julia_eval("r.Q"))
  out$involv = julia_eval("r.involv")
  out$level_lam_total = julia_eval("r.level_lam_total")
  
  attr(out, "class") = "decompMMTD"
  
  if (print) julia_eval("print(r)")  
  return(out)
}

decomp_to_julia = function(x) {
  UseMethod("decomp_to_julia")
}

decomp_to_julia.decompMMTD = function(decmp) {
  stopifnot( class(x) == "decompMMTD" )
  
  julia_assign("dcp", decmp)
  julia_command(" dcpJ = DecompMMTDg( dcp[:Omega], 
             dcp[:totalweight_tol], dcp[:levels], 
             [ [Int.(dcp[:configs][i])...] for i in 1:length(dcp[:configs]) ], 
             [ dcp[:Q][i] for i in 1:length(dcp[:Q]) ], 
             dcp[:involv], dcp[:level_lam_total]) ",
                show_value=FALSE)
  
  return(NULL)
}

reconstruct = function(x) {
  UseMethod("reconstruct")
}

reconstruct.decompMMTD = function(decmp) {
  decomp_to_julia.decompMMTD(decmp)
  julia_eval("dcpJ.Omega")
}

print.decompMMTD = function(decmp) {
  decomp_to_julia.decompMMTD(decmp)
  julia_command("print(dcpJ)")
}

### test
# Omeg = array( c(0.5, 0.5, 0.9, 0.1, 0.8, 0.2, 0.4, 0.6), dim=c(2,2,2) )
# 
# x = decompMMTD(Omeg)
# x = decompMMTD(Omeg, random=TRUE)
# 
# # decomp_to_julia(x)
# 
# print(x)
# reconstruct(x)




