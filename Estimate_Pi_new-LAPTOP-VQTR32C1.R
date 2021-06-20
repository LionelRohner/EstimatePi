# Library -----------------------------------------------------------------

library(rbenchmark)

library(fitdistrplus)


# Empirical Methods -------------------------------------------------------

# generate uniformly distributed points
generate_points <- function(n){
  
  # put them in a matrix
  XY <- matrix(runif(2*n), nrow = n, ncol = 2)
  
  return(XY)

}

# calculate vector norm and check whether point is within radius (i.e. < 1)
get_distance <- function(XY){
  
  # calculate vector norm and test if its bigger than the radius (i.e. 1)
  dist <- sqrt(XY[,1]^2 + XY[,2]^2) < 1
  return(dist)
}

# approximate pi using 4*points_within_circle/total_points
approx_pi <- function(dist){
  # pi*r^2 / 4 == numb_circle / numb_total (divided by for because we have just one quadrant)
  return(4*sum(dist)/length(dist))
}



# Resampling Methods ------------------------------------------------------

# calculate ratio of points within
calc_ratio <- function(n, samplingSize, plot){
  distMatrix <- matrix(get_distance(generate_points(n)), nrow = samplingSize)
  ratioVector <- rowSums(distMatrix)/ncol(distMatrix)
  
  if (plot){
    hist(ratioVector, freq = F)
  }
  
  return(ratioVector)
}

# fit a gamma distribution to ratio data set and sample from gamma
generate_gamma <- function(ratioVector, outputLength, plot){
  thetaGamma <- fitdistr(ratioVector, "gamma")$estimate 
  
  set.seed(sample(1:1e6,1))
  
  # histogram
  if (plot){
    hist(rgamma(outputLength,shape = thetaGamma[1],rate = thetaGamma[2]),
         add = T, col = rgb(0.9,0.1,0.1,0.2), freq = F)
  }
  return(rgamma(outputLength,shape = thetaGamma[1],rate = thetaGamma[2]))  
}

# from: https://stats.stackexchange.com/questions/376634/how-to-pick-starting-parameters-for-massfitdist-with-the-beta-distribution
beta_mom <- function(x) {
  
  m_x <- mean(x, na.rm = TRUE)
  s_x <- sd(x, na.rm = TRUE)
  
  alpha <- m_x*((m_x*(1 - m_x)/s_x^2) - 1)
  beta <- (1 - m_x)*((m_x*(1 - m_x)/s_x^2) - 1)
  
  return(list(shape1 = alpha, shape2 = beta))
  
}

# fit a beta distribution to ratio data set and sample from this distribution
generate_beta <- function(ratioVector, outputLength, plot){
  
  # remove 1s and 0s (the latter is very unlikely)
  ratioVector <- ratioVector[which(ratioVector < 1)]

  
  # use methods of moments to get an estimate for alpha (shape1) and beta (shape2)
  # fitdistr needs starting values for beta-distribution
  thetaBetaStartValues <- beta_mom(ratioVector)
  
  # MLE of beta distribution given ratioVector data 
  thetaBeta <- fitdistr(ratioVector, "beta", start=thetaBetaStartValues)$estimate 
  
  
  # set.seed(sample(1:1e6,1))
  
  # histogram
  if (plot){
    hist(rbeta(outputLength,shape1 = thetaBeta[1],shape2 = thetaBeta[2]),
         add = T, col = rgb(0.9,0.1,0.1,0.2), freq = F)
  }

  # generate data from the beta-distribution with MLE-estimates
  return(rbeta(outputLength, shape1 = thetaBeta[1], shape2 = thetaBeta[2]))  
}


# fit a beta distribution to ratio data set and sample from this distribution
get_beta_dist <- function(ratioVector){
  
  # remove 1s and 0s (the latter is very unlikely)
  ratioVector <- ratioVector[which(ratioVector < 1)]
  
  
  # use methods of moments to get an estimate for alpha (shape1) and beta (shape2)
  # fitdistr needs starting values for beta-distribution
  thetaBetaStartValues <- beta_mom(ratioVector)
  
  # MLE of beta distribution given ratioVector data 
  return(fitdistr(ratioVector, "beta", start=thetaBetaStartValues)$estimate) 
}


# approximate pi (same as above).  
approx_pi_resample <- function(randVec){
  withinCircle <- mean(randVec)
  outsideCircle <- 1 - withinCircle 
  
  return(4*(withinCircle/(withinCircle+outsideCircle)))
}
  

# Main Function -----------------------------------------------------------

# empirical version just generates data between 0 and 1 and calculates the
# to the radius of the unit circle.
# More details here : https://en.wikipedia.org/wiki/Approximations_of_%CF%80#Summing_a_circle's_area
estimate_pi_empirical <- function(n){
  dist <- get_distance(generate_points(n = n))
  return(approx_pi(dist = dist))
}


# this function uses a probabilistic approach. The ratio of points that are
# within the radius is sampled from a fitted gamma distribution. Currently, the 
# mean of randomly generated data from this vector is used to 
estimate_pi_resampled <- function(n,
                                  outputLength = 1e6,
                                  samplingSize = 1e5,
                                  distr = "beta",
                                  plot = F){
  
  if(!require(fitdistrplus)){ 
    message("Install 'fitdistrplus' first!")
    return(NULL)
  }
  

  # generate ratios from n
  # create gamma distribution
  if(distr == "beta"){
    rand <- generate_beta(calc_ratio(n, samplingSize = samplingSize, plot = plot),
                        plot = plot, outputLength = outputLength)
    
  } else {
    rand <- generate_gamma(calc_ratio(n, samplingSize = samplingSize, plot = plot),
                        plot = plot, outputLength = outputLength)
    
  }
  
  # use resampled data to approx pi
  return(approx_pi_resample(rand))
}




MCMC_Pi <- function(nInit = 1e6,nSD = 1000, nIter){
  
  # 0.) initialize result vector and get a value for d
  x <- rep(0,nIter)
  
  d <- sd(calc_ratio(nInit,nSD,F))
  
  # 1.) create a prior distribution
  
  # beta parameter from 1e6 ratios >> prior distribution params
  thetaBetaMLE <- get_beta_dist(calc_ratio(1e6,1e4,F))
  
  # mean of beta distribution
  meanBeta <- unname(1/(1+(thetaBetaMLE[2]/thetaBetaMLE[1])))
  
  pbeta(meanBeta,shape1 = thetaBetaMLE[1], shape2 = thetaBetaMLE[2])
  
  # 2.) Propose first move (start with mean of prior)
  x[1] <- meanBeta
  
  # 3.) Initiate loop
  
  for (iter in 2:nIter){
    
    # 4.) propose a move with uniform proposal kernel
    x[iter] <- runif(n = 1, min = x[iter-1]-d/2, max=x[iter-1]+d/2) 
    
    # 3.) Compute accuracy, accept move if accuracy is better else stay
    current <- abs(approx_pi_resample(x[iter]) - pi)
    previous <- abs(approx_pi_resample(x[iter-1]) - pi)
    
    if (current < previous){
      next
    } else {
      x[iter] <- x[iter-1]
    }
    
  }
  return(x)
}


# Auxiliary Functions -----------------------------------------------------

# the smaller the number the better
accuracy_pi_estimate <- function(pi_estimate){
  return(abs(pi_estimate-pi))
}

scoring <- function(res){
  if(res > 0.1){
    return(0)
  } else if (res > 0.01){
    return(1)
  } else if(res > 0.001){
    return(2)
  } else if(res > 0.0001){
    return(3)
  } else if(res > 0.00001){
    return(4)
  } else if(res > 0.000001){
    return(5)
  } else{
    return(6)
  } 
}

# compare time and score of accuracy, not for MCMC as we only consider the last values
test_accuracy <- function(n,
                          nIter,
                          type,
                          samplingSize = 1e4,
                          outputLength = 1e6,
                          distr = "beta"){
  
  # Scoring Output
  message("Scoring Calculation: Difference = estimated pi minus real pi \n")
  message("Scores:")
  message("Difference > 0.1 : Score = 0")
  message("Difference > 0.01 : Score = 1")
  message("Difference > 0.001 : Score = 2")
  message("Difference > 0.0001 : Score = 3")
  message("Difference > 0.00001 : Score = 4")
  message("Difference > 0.000001 : Score = 5")
  message("Difference < 0.000001 : Score = 6 \n")

  
  # create vector for accuracy measure
  vecAccuracy = c() 
  
  # loop for nIter
  if (type == "empirical"){
    message("Method : Empirical \n")
    for (i in 1:nIter){
      estimate <- accuracy_pi_estimate(estimate_pi_empirical(n))
      score <- scoring(estimate)
      vecAccuracy = c(vecAccuracy, score)
    }

  } else {
    message("Method : Resampled \n")
    for (i in 1:nIter){
      estimate <- accuracy_pi_estimate(estimate_pi_resampled(n,
                                                             samplingSize = samplingSize,
                                                             outputLength = outputLength,
                                                             distr = distr))
      score <- scoring(estimate)
      vecAccuracy = c(vecAccuracy, score)
    }

  }
  # output formatting
  DF <- data.frame(nIter, mean(vecAccuracy), sd(vecAccuracy), min(vecAccuracy), max(vecAccuracy))
  colnames(DF) <- c("Iterations", "Mean Score", "SD Score", "Min Score", "Max Score")
  
  message("Accuracy of the Pi Estimate: ")
  return(DF)
   
}

# descriptive statistics of the pi estimate collected from several iterations
mean_estimate <- function(n,
                          nIter,
                          type,
                          samplingSize = 1e4,
                          outputLength = 1e6,
                          distr = "beta"){
  
  vecAccuracy = c() 
  
  if (type == "empirical"){
    for (i in 1:nIter){
      estimate <- estimate_pi_empirical(n)
      vecAccuracy = c(vecAccuracy, estimate)
    }

  } else {
    for (i in 1:nIter){
      estimate <- accuracy_pi_estimate(estimate_pi_resampled(n,
                                                             samplingSize = samplingSize,
                                                             outputLength = outputLength,
                                                             distr = distr))
      vecAccuracy = c(vecAccuracy, estimate)
    }

  }
  # output formatting
  DF <- data.frame(nIter, mean(vecAccuracy), sd(vecAccuracy), min(vecAccuracy), max(vecAccuracy))
  colnames(DF) <- c("Iterations", "Mean Score", "SD Score", "Min Score", "Max Score")
  return(DF)
}


# Test functions ----------------------------------------------------------

# empirical
estimate_pi_empirical(1e6)

# resampled
estimate_pi_resampled(n = 1e6,
                      outputLength = 1e6,
                      samplingSize = 1e4,
                      plot = F,
                      distr = "beta")

# MCMC (super accurate), but fake
posterior <- MCMC_Pi(nIter = 1000)
plot(4*posterior, type = "l")
abline(h=pi, col = "firebrick")

approx_pi_resample(tail(posterior,n = 1))



# BenchMark (Speed) -------------------------------------------------------


n = 1e6

test <- benchmark("empirical" = {estimate_pi_empirical(n)},
                  "resampled" = {estimate_pi_resampled(n)},
                  "MCMC"      = {MCMC_Pi(nIter = 1e5)}
                  replications = 10)

test$meanTime <- test$elapsed/test$replications
test
 

# BenchMark (Accuracy) ----------------------------------------------------


# Accurcay Scores
gamma <- test_accuracy(n = 1e6, samplingSize = 1e4, outputLength = 1e6, distr = "gamma", nIter = 1000, type = "resampled")
beta <- test_accuracy(n = 1e6, samplingSize = 1e4, outputLength = 1e6, distr = "beta", nIter = 1000, type = "resampled")
norm <- test_accuracy(n = 1e6, samplingSize = 1e4, outputLength = 1e6, distr = "beta", nIter = 1000, type = "empirical")


gamma
beta
norm


# Accuracy of Estimate (raw)
mean_estimate(n = 1e6, nIter = 100, type = "empirical")
mean_estimate(n = 1e6, nIter = 100, type = "resampled")


# As expected the accuracies are more or less the same


# Garbage Code ------------------------------------------------------------

# find best params (single)
RES <- c()
seqN <- seq(1e3,1e5,1e3)
for (i in seqN){
  piVec <- c()
  for (j in 1:100){
    pi <- estimate_pi_resampled(n = 1e6,
                                outputLength = 1e6,
                                samplingSize = i,
                                plot = F)
    piVec <- c(piVec,pi) 
  }
  cat(sprintf("\r%.3f%%", (i / (1e5) * 100)))
  RES <- rbind(RES,abs(mean(piVec)-pi))
}; rownames(RES) <- seqN
