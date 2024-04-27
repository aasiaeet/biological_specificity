# Now that we have a posterior over a grid of $(\alpha_1, \beta_1, \alpha_2, \beta_2)$, we can sample from the joint posterior. We start by sampling from $P(\alpha_1, \beta_1, \alpha_2, \beta_2|y)$:

library(dplyr)

sampPdf <- function(n, pdf){
  cdf <- cumsum(pdf)
  unlist(sapply(runif(n), function(i) {1 + sum(cdf < i)}))
}

samplePosteriorMixture <- function(br, nsamp){
  P <- br@betaPost
  #Compute distribution of x1=logit(alpha1/beta1) (betaPost already normalized)
  marPdf <- apply(P, 1, sum)
  sam <- sampPdf(nsamp, marPdf)
  #Sampling the rest of dimensions
  #Each row of sam is a sample and its columns grow in the loop
  for(i in 2:4){
    tabSam <- as.data.frame(as.data.frame(sam) %>% group_by_all() %>% count)
    sam <- c()
    for(r in 1:nrow(tabSam)){
      R <- as.numeric(tabSam[r,]);
      cnt <- R[i]; #last col is the count
      inds <- list()
      for(j in 1:4){
        if(j < i) {inds[[j]] <- R[j];}
        else {inds[[j]] <- 1:dim(P)[j];}
      }
      subX <- P[unlist(inds[[1]]),unlist(inds[[2]]) ,unlist(inds[[3]]) ,unlist(inds[[4]])]
      
      if(is.null(dim(subX))){marPdf <- subX / sum(subX)} #if it is vector
      else{marPdf <- apply(subX, 1 , sum) / sum(subX)}
      
      tSamp <- sampPdf(cnt, marPdf)
      sam <- rbind(sam, cbind(t(matrix(rep(R[1:(i-1)], cnt), nrow = i-1)), tSamp))
    }
  }
  #Note that sample is actually sample index and we need to read out the params:
  return(cbind(br@x1[sam[,1]], br@y1[sam[,2]], br@x2[sam[,3]], br@y2[sam[,4]]))
}


# Next, we sample $P(\pi|\alpha_1, \beta_1, \alpha_2, \beta_2, y)$. 
samplePosteriorPi <- function(K, N, dPi=1/10, betaPar){
  ratios1 <- allLogRatios(K, N, betaPar[,1], betaPar[,2], isgrid = F)
  ratios2 <- allLogRatios(K, N, betaPar[,3], betaPar[,4], isgrid = F)
  unNormPi <- posteriorOfPi(ratios1$vals, ratios2$vals, dPi, isgrid = F)
  piPdf <- apply(unNormPi, 1, function(x) x / sum(x)); #each col of piPdf is a Pdf
  piSeq <- seq(dPi, 1 - dPi, by = dPi)
  piSamp <- apply(piPdf, 2, function(x) piSeq[sampPdf(1, x)])
  #TODO: add jitter to the samples
  # piSamp <- apply(piPdf, 2, function(x) piSeq[sampPdf(1, x)]+jitter(x, amount = dPi))
  return(list(piSamp = piSamp, ratios1 = ratios1, ratios2 = ratios2))
}


# Then, we need to sample the posterior of the latent variable $P(Z|\alpha_1, \beta_1, \alpha_2, \beta_2, \pi, y)$. 
#returns a matrix of nsamp by T(i.e., br@n)
samplePosteriorZ <- function(output){
  #vals are nsamp by T dimensional 
  
  unNorm1 <- as.data.frame(t(output$ratios1$vals + log(output$piSamp)))
  unNorm2 <- as.data.frame(t(output$ratios2$vals + log(1 - output$piSamp)))
  nT <- nrow(unNorm1)
  sampZ <- mapply(function(x,y) {
    pMx <- pmax(x,y)
    tx <- exp(x - pMx)
    ty <- exp(y - pMx)
    normZ <- cbind(tx, ty) / (tx + ty)
    return(c(2,1)[(runif(nT) < normZ[,1])+1])
  },  
  unNorm1, unNorm2)
  
  
  return(t(sampZ))
}

# Finally, we can sample the rates $P(\theta_j|\alpha_1, \beta_1, \alpha_2, \beta_2, \pi, z, y)$ and compute the Fisher's Score The below function calls the previous ones and completes the sampling: 

samplePosteriorRates <- function (br, dPi=1/10, nsamp = 2000) {
  #samples are in transformed scale (logit and log)
  sampleBetaParams <- samplePosteriorMixture(br, nsamp)
  #generates many useful intermediate results that returns
  output <- samplePosteriorPi(br@k, br@n, dPi, sampleBetaParams)
  samplePi <- output$piSamp
  #returns a nsamp * length(br@n) matrix
  sampleZ <- samplePosteriorZ(output)
  
  alpha1 <- output$ratios1$params$alpha
  beta1 <- output$ratios1$params$beta
  alpha2 <- output$ratios2$params$alpha
  beta2 <- output$ratios2$params$beta
  theta <- matrix(NA, nrow = nsamp, ncol = length(br@n))
  for (i in 1:length(br@n)) {
    comp1 <- sampleZ[,i] == 1
    comp2 <- sampleZ[,i] == 2
    theta[comp1, i] <- rbeta(sum(comp1), alpha1 + br@k[i], beta1 + br@n[i] - br@k[i])
    theta[comp2, i] <- rbeta(sum(comp2), alpha2 + br@k[i], beta2 + br@n[i] - br@k[i])
  }
  c1 <- theta[which(sampleZ == 1, arr.ind = T)]
  c2 <- theta[which(sampleZ == 2, arr.ind = T)]
  theta <- data.frame(theta)
  
  brand <- names(br@n)
  if (!is.null(brand)) {
    if (all(!is.na(brand))) {
      colnames(theta) <- names(br@n)
    }
  }
  
  list(xy = data.frame(alpha1 = alpha1, beta1 = beta1,
                       alpha2 = alpha2, beta2 = beta2,
                       piV = samplePi),
       fisherC = (mean(c2) - mean(c1))^2 / (var(c1)+var(c2)),
       z = sampleZ, theta = theta,
       clusters = list(c1 = c1, c2 = c2))
}

