
# Library -----------------------------------------------------------------

library(rbenchmark)

library(fitdistrplus)


# Empirical Methods -------------------------------------------------------

generate_points <- function(n){
  
  # put them in a matrix
  XY <- matrix(runif(2*n), nrow = n, ncol = 2)
  
  return(XY)

}

get_distance <- function(XY){
  
  # calculate vector norm and test if its bigger than the radius (i.e. 1)
  dist <- sqrt(XY[,1]^2 + XY[,2]^2) < 1
  return(dist)
}


approx_pi <- function(dist){
  # pi*r^2 / 4 == numb_circle / numb_total (divided by for because we have just one quadrant)
  return(4*sum(dist)/length(dist))
}



# Resampling Methods ------------------------------------------------------



calc_ratio <- function(n, samplingSize, plot){
  distMatrix <- matrix(get_distance(generate_points(n)), nrow = samplingSize)
  ratioVector <- rowSums(distMatrix)/ncol(distMatrix)
  
  if (plot){
    hist(ratioVector)
  }
  
  return(ratioVector)
}

generate_gamma <- function(ratioVector, outputLength, plot){
  thetaGamma <- fitdistr(ratioVector, "gamma")$estimate 
  
  # histogram
  if (plot){
    hist(rgamma(outputLength,shape = thetaGamma[1],rate = thetaGamma[2]),
         add = T, col = rgb(0.9,0.1,0.1,0.2))
  }
  return(rgamma(outputLength,shape = thetaGamma[1],rate = thetaGamma[2]))  
}

approx_pi_resample <- function(rgammaVec, accuracyHeuristic){
  withinCircle <- mean(rgammaVec) * accuracyHeuristic
  outsideCircle <- accuracyHeuristic - withinCircle 
  
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

# more or less exactly the same speed. The idea was to average the dist matrix
# values by generating many small matrices instead of one big, since matrix/vecs 
# bigger than 1e7 get really slow to generate.
estimate_pi_v2 <- function(n){
  pis <- c()
  for (i in 1:1000){
    points <- generate_points_distance(n = n%/%1000)
    dist <- get_distance(points)
    approxed_pi <- approx_pi(dist = dist)
    pis <- c(pis,approxed_pi)
  }
  return(mean(pis))
}


# this function uses a probabilistic approach. The ratio of points that are
# within the radius is sampled from a fitted gamma distribution. Currently, the 
# mean of randomly generated data from this vector is used to 
estimate_pi_resampled <- function(n,
                                  outputLength = 1000,
                                  samplingSize = 100,
                                  accuracyHeuristic = 1e6,
                                  plot = F){
  
  if(!require(fitdistrplus)){
    message("Install 'fitdistrplus' first!")
    return(NULL)
  }
  
  # ratio of at least 100 data should be used to get ratios
  if(n/samplingSize < 100){
    message("Choose bigger n or smaller sampling size")
    return(NULL)
  }
  
  # generate ratios from n
  # create gamma distribution
  rG <- generate_gamma(calc_ratio(n, samplingSize = samplingSize, plot = plot),
                       plot = plot, outputLength = outputLength)
  
  # use resampled data to approx pi
  return(approx_pi_resample(rG, accuracyHeuristic = accuracyHeuristic))
}



# Test functions ----------------------------------------------------------



estimate_pi_resampled(n = 100000,
                      outputLength = 1e6,
                      samplingSize = 1000,
                      accuracyHeuristic = 1e9)





# BenchMark (Speed) -------------------------------------------------------



n = 1e7

test <- benchmark("v1" = {estimate_pi_empirical(n)},
                  "v2" = {estimate_pi_v2(n)},
                  "v3" = {estimate_pi_resampled(n)},
                  replications = 2)

test$meanTime <- test$elapsed/test$replications
test
 

# BenchMark (Accuracy) ----------------------------------------------------


# WRITE A FUNCTION THAT COMPARES ACCURACY OF PI

accuracy_pi_estimate <- function