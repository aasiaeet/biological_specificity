# To save processing time, we compute the gamma ratios and store them. 
# Given a $(x, y) = (\log(\frac{\alpha}{\beta}), \log(\alpha+\beta))$ point and two vectors of number or successes (K) and number of trials (N) per tissue, the below code computes the ratio $\frac{\Gamma(\alpha+\beta)}{\Gamma(\alpha)\Gamma(\beta)} \frac{\Gamma(\alpha+y[t])\Gamma(\beta+n[t]-y[t])}{\Gamma(\alpha+\beta+n[t])}$ for all $T$ tissues and returns a row of those values. 
# We keep track of corresponding $\alpha$ and $\beta$ parameters in the `params` matrix and their gamma ratios in the `vals` matrix. 
# We also save the grid coordinates along with the index of those parameters in the given parameter set. 
# Note that we support both grid and pairs of parameters as inputs. 
# Finally, to avoid numerical issues we compute the log of the ratios. 

logRatios <- function(x, y, K, N){
  # transform back x and y to alpha and beta scale 
  alpha <- exp(x + y)/(1 + exp(x))
  beta <- exp(y)/(1 + exp(x))
  # compute the base ratio which is common among all tissues 
  baseRatio <- lgamma(alpha + beta) - (lgamma(alpha) + lgamma(beta))
  
  outputRow <- c()
  for(t in 1:length(K)){
    outputRow <- c(outputRow, baseRatio + 
                     (lgamma(alpha + K[t]) + lgamma(beta + N[t] - K[t]) - lgamma(alpha + beta + N[t])))
  }
  return(c(x, y, alpha, beta, outputRow))
}

allLogRatios <- function(K, N, x, y, isgrid = T){
  gammaTable <- c()
  if(isgrid) { #grid of params
    for (i in 1:length(x)) {
      for (j in 1:length(y)) {
        gammaTable <-
          rbind(gammaTable, c(i, j, logRatios(x[i], y[j], K, N)))
      }
    }
  } else { #arbit. params
    stopifnot(length(x) == length(y))
    for (i in 1:length(x)){
      gammaTable <-
        rbind(gammaTable, c(i, i, logRatios(x[i], y[i], K, N)))
    }
  }
  
  gammaTable <- as.data.frame(gammaTable)
  colnames(gammaTable)[1:6] <- c("xIndex", "yIndex", 
                                 "logitMu", "logSize", 
                                 "alpha", "beta")
  gRatios <- list();
  gRatios$index  <- gammaTable[,1:2]
  gRatios$coords <- gammaTable[,3:4]
  gRatios$params <- gammaTable[,5:6]
  gRatios$vals <- gammaTable[,-(1:6)]
  return(gRatios)
}