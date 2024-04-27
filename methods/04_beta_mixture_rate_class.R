# Since we do not have any prior information about the balance between the two components and their parameters for that matter, the `guessCenter` is as before and generates one set of parameters which we use for both Betas and define the box around that value. 
myGuessCenter <- function (v) 
{
  m <- mean(v)
  s2 <- var(v)
  temp <- m * (1 - m)/s2 - 1
  if (temp < 0) 
    temp <- 0.001
  a <- m * temp
  b <- (1 - m) * temp
  X <- log(a/b)
  Y <- log(a + b)
  list(x1 = X, x2 = X, y1 = Y, y2 = Y, a1 = a, b1 = b, a2 = a, b2 = b)
}

# Now we are ready for sampling from the posterior of the hierarchical model. 
# To do so, we first need to sample from the posterior of the hyperparameters; the rest is easier. 
# So in the constructor of the BetaMixtureRates object below, we make the approximation of of the transformed Beta posteriors as step functions over a sufficiently large grid and save it in the object. 
# The rest of the computation will be done whenever the `sampleFromPosterior` is called.
BetaMixtureRates <- function(K, N, dPi=1/10, massTrsh=0.87, resolution=10)
{
  
  V <- K/N
  guess <- myGuessCenter(V)
  W <- 3
  deltas <- list(lx1=W, rx1=W, ly1=W, ry1=W, lx2=W, rx2=W, ly2=W, ry2=W)
  logitMu1 <- logitMu2 <- logSize1 <- logSize2 <- list(mn=0, mx=0)
  cnt <- 0
  L <- resolution
  x1 <- y1 <- x2 <- y2 <- c()
  repeat{
    cnt <- cnt + 1;
    if(cnt == 3) break;
    print(paste("Finding a big enough grid, trial", cnt))
    logitMu1[1:2] <- c(guess$x1 - deltas$lx1, guess$x1 + deltas$rx1)
    logSize1[1:2] <- c(guess$y1 - deltas$ly1, guess$y1 + deltas$ry1)
    logitMu2[1:2] <- c(guess$x2 - deltas$lx2, guess$x2 + deltas$rx2)
    logSize2[1:2] <- c(guess$y2 - deltas$ly2, guess$y2 + deltas$ry2)
    
    x1 <- seq(logitMu1$mn, logitMu1$mx, length=L)
    y1 <- seq(logSize1$mn, logSize1$mx, length=L)
    x2 <- seq(logitMu2$mn, logitMu2$mx, length=L)
    y2 <- seq(logSize2$mn, logSize2$mx, length=L)
    
    
    gRatios1 <- allLogRatios(K, N, x1, y1)
    gRatios2 <- allLogRatios(K, N, x2, y2)
    H <- logIntegralOfPi(gRatios1$vals, gRatios2$vals, dPi)
    P <- logHyperPrior(gRatios1$params, gRatios2$params, sigma = 1/10)
    results <- exp(H + P - max(H + P))
    
    maxPost <- max(results) #this is always 1 
    print(paste("How many violate?", sum(results > 0.0087 * maxPost), "/", L^4))
    
    faces <- 100*c(lx1 =max(results[gRatios1$coords$logitMu == logitMu1$mn, ]),
                   rx1 =max(results[gRatios1$coords$logitMu == logitMu1$mx, ]),
                   ly1 =max(results[gRatios1$coords$logSize == logSize1$mn, ]),
                   ry1 =max(results[gRatios1$coords$logSize == logSize1$mx, ]),
                   lx2 =max(results[, gRatios2$coords$logitMu == logitMu2$mn]),
                   rx2 =max(results[, gRatios2$coords$logitMu == logitMu2$mx]),
                   ly2 =max(results[, gRatios2$coords$logSize == logSize2$mn]),
                   ry2 =max(results[, gRatios2$coords$logSize == logSize2$mx]))
    
    if (all(faces < massTrsh * maxPost)) break #if nothing on edges violates
    if (sum(results > 0.0087 * maxPost)/ L^4 < 0.01) break #if less than 1% of points violate 
    temp <- names(which(faces >= massTrsh * maxPost))
    for (n in temp) deltas[[n]] <- deltas[[n]] * 2 # make a wider box
  }
  results <- results/sum(results)
  
  
  #Moving everything into a 4D array
  betaPost <- array(0, rep(L,4))
  for(u in 1:L^2){
    for(v in 1:L^2){
      betaPost[gRatios1$index[u, 1],
               gRatios1$index[u, 2],
               gRatios2$index[v, 1],
               gRatios2$index[v, 2]] <- results[u, v]
    }
  }
  
  bestD <- which(betaPost == max(betaPost), arr.ind = TRUE)
  cat("Best (alpha1, beta1)", unlist(xy2ab(x1[bestD[1]], y1[bestD[2]])), '\n')
  cat("Best (alpha2, beta2)", unlist(xy2ab(x2[bestD[3]], y2[bestD[4]])), '\n')
  
  
  
  new("BetaMixtureRates", 
      k = K, n = N, l = L,
      x1 = x1, y1 = y1,
      x2 = x2, y2 = y2, 
      betaPost = betaPost)
}

setClass("BetaMixtureRates", slots = list(k="numeric", n="numeric", l="numeric", 
                                          x1="numeric", y1="numeric", 
                                          x2="numeric", y2="numeric",
                                          betaPost="array"))

ab2xy <- function(a,b){
  out <- list()
  out$x <- log(a/b)
  out$y <- log(a+b)
  return(out)
}

xy2ab <- function(x,y){
  out <- list();
  out$a <- exp(x + y)/(1 + exp(x))
  out$b <- exp(y)/(1 + exp(x))
  return(out)
}