# The two parameter sets $(\alpha_1, \beta_1)$ and $(\alpha_2, \beta_2)$ are assumed to be independent and with $\propto (\alpha+\beta)^{-5/2}$ prior. 
# Translating this to the grid parameter space we get an extra term from the Jacobian $P(\log(\frac{\alpha}{\beta}), \log(\alpha + \beta)) \propto \alpha \beta (\alpha+\beta)^{-5/2}$ as the prior. 
# For each point on the grid of parameters $(\alpha_1, \beta_1, \alpha_2, \beta_2)$ we compute the log hyperprior below. 
# Note that we set the prior to zero (or log of it it -Inf [or the smallest possible number: `.Machine$longdouble.min.exp`]) when the mean1 > mean2 and when either of the Beta distributions become bimodal.

logHyperPrior <- function(P1, P2, sigma= 1){
  logPriors <- array(NA, dim=c(nrow(P1),nrow(P2)))
  lp1 <- nrow(P1); lp2 <- nrow(P2);
  cnt <- 0;
  minLog <- .Machine$longdouble.min.exp
  # minLog <- .Machine$double.min.exp
  for(i in 1:lp1){
    for(j in 1:lp2){
      a1 <- P1$alpha[i]; b1 <- P1$beta[i];
      a2 <- P2$alpha[j]; b2 <- P2$beta[j];
      if((a1 < 1 &&  b1 < 1) || (a2 < 1 && b2 < 1)){
        logPriors[i, j] <- minLog; next;
      }
      size1 <- a1 + b1;   size2 <- a2 + b2;
      mean1 <- a1 / size1;mean2 <- a2 / size2;
      if(mean1 < mean2){ #to help with identifiablity in later Fisher's criterion computation
        logPriors[i, j] <-
          (log(a1) + log(b1) + log(a1 + b1) * (-5/2)) +
          (log(a2) + log(b2) + log(a2 + b2) * (-5/2)) #+
        # (-1 / (2 * sigma)) * ((mean1 - mean2) ^ 2 + (size1^(-1/2) - size2^(-1/2))^2)
      } else {logPriors[i, j] <- minLog}
    }
    # cat(cnt, '\n')
  }
  return(logPriors)
}