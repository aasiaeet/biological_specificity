# Because we have a mixture model, to compute the posterior of the Beta model hyperparameters we need to marginalize out (integrate out) the other parameter $\pi$ where $(\pi, 1-\pi)$ is the mixture proportion vector. 
# To do this numerically, we consider discrete points with dPi $=\Delta\pi=1/100$ distance from each other and use the below function to compute logarithm of $\int_0^1 P(\pi) \prod_t \Big( \pi \frac{\Gamma(\alpha_1+\beta_1)}{\Gamma(\alpha_1)\Gamma(\beta_1)}  \frac{\Gamma(\alpha_1+y[t])\Gamma(\beta_1+n[t]-y[t])}{\Gamma(\alpha_1+\beta_1+n[t])} +  (1-\pi) \frac{\Gamma(\alpha_2+\beta_2)}{\Gamma(\alpha_2)\Gamma(\beta_2)}  \frac{\Gamma(\alpha_2+y[t])\Gamma(\beta_2+n[t]-y[t])}{\Gamma(\alpha_2+\beta_2+n[t])} \Big)$.
# Note that taking log and subtracting maximum and then exponentiate is the standard approach to avoid numerical instability. 

logIntegralOfPi <- function(val1, val2, dPi = 1/10)
{
  # Sum over the pi dimension of a 3D array.
  return(log(rowSums(posteriorOfPi(val1, val2, dPi), dims = 2)))
}

posteriorOfPi <- function(V1, V2, dPi = 1/10, isgrid = T){
  if (!isgrid) {stopifnot(nrow(V2) == nrow(V1))}
  
  piV <- seq(dPi, 1 - dPi, by = dPi)
  dims <- if(isgrid) c(nrow(V1), nrow(V2), length(piV)) else c(nrow(V1), length(piV));
  sumLogs <- array(rep(NA, prod(dims)), dim = dims)
  cat("loop over pi: ")
  for (k in 1:length(piV)) {
    cat(k, " ") #7s for all inner loops
    for (i in 1:nrow(V1)) {
      range <- if(isgrid) 1:nrow(V2) else i
      for (j in range) {
        accSum <- 0
        for (t in 1:ncol(V1)) {
          accSum <- accSum + logTrick(V1[i, t], V2[j, t], piV[k])
        }
        if (isgrid) {sumLogs[i, j, k] <- accSum}
        else{sumLogs[i, k] <- accSum}
      }
    }
  }
  cat("\n")
  if (isgrid) {outcome <- exp(sumLogs - max(sumLogs))} #log trick for the whole matrix
  else{outcome <- t(apply(sumLogs, 1, function(x) exp(x - max(x))))} #log trick for each row
  return(outcome)
}

# Note that in the `posteriorOfPi` function I take care of two similar computations in one place. 
# If the input is a grid and we are integrating the $\pi$ out we have one case and we should use the exp(log-max) trick for the whole matrix. 
# But there is a similar computation needed for computing the posterior of pi where $(\alpha_1, \beta_1)$ and $(\alpha_2, \beta_2)$ are already sampled (no grid) and the exp(log-max) should be done for each row.

#log trick for numerical stability
logTrick <- function(v1, v2, piV){
  t1 <- v1 + log(piV)
  t2 <- v2 + log(1 - piV)
  mnT <- min(t1, t2); mxT <- max(t1,t2);
  return(mxT + log(1 + exp(mnT - mxT)))
}