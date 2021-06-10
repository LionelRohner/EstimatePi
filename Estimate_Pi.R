
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
  
  # histogram
  if (plot){
    hist(rgamma(outputLength,shape = thetaGamma[1],rate = thetaGamma[2]),
         add = T, col = rgb(0.9,0.1,0.1,0.2), freq = F)
  }
  return(rgamma(outputLength,shape = thetaGamma[1],rate = thetaGamma[2]))  
}

# approximate pi (same as above).  
approx_pi_resample <- function(rgammaVec){
  withinCircle <- mean(rgammaVec)
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

# more or less exactly the same speed. The idea was to average the dist matrix
# values by generating many small matrices instead of one big, since matrix/vecs 
# bigger than 1e7 get really slow to generate.
estimate_pi_v2 <- function(n){
  pis <- c()
  for (i in 1:1000){
    points <- generate_points(n = n%/%1000)
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
                                  outputLength = 1e6,
                                  samplingSize = 1e5,
                                  plot = F){
  
  if(!require(fitdistrplus)){ 
    message("Install 'fitdistrplus' first!")
    return(NULL)
  }
  
  # # ratio of at least 100 data should be used to get ratios
  # if(n/samplingSize < 100){
  #   message("Choose bigger n or smaller sampling size")
  #   return(NULL)
  # }
  
  
  # generate ratios from n
  # create gamma distribution
  rG <- generate_gamma(calc_ratio(n, samplingSize = samplingSize, plot = plot),
                       plot = plot, outputLength = outputLength)
  
  # use resampled data to approx pi
  return(approx_pi_resample(rG))
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

test_accuracy <- function(n, nIter, type){
  
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

  
  
  vecAccuracy = c() 
  
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
      estimate <- accuracy_pi_estimate(estimate_pi_resampled(n))
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


mean_estimate <- function(n, nIter, type){
  
  vecAccuracy = c() 
  
  if (type == "empirical"){
    for (i in 1:nIter){
      estimate <- estimate_pi_empirical(n)
      vecAccuracy = c(vecAccuracy, estimate)
    }

  } else {
    for (i in 1:nIter){
      estimate <- estimate_pi_resampled(n)
      vecAccuracy = c(vecAccuracy, estimate)
    }

  }
  # output formatting
  DF <- data.frame(nIter, mean(vecAccuracy), sd(vecAccuracy), min(vecAccuracy), max(vecAccuracy))
  colnames(DF) <- c("Iterations", "Mean Score", "SD Score", "Min Score", "Max Score")
  return(DF)
}


# Test functions ----------------------------------------------------------

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




estimate_pi_resampled(n = 1e6,
                      outputLength = 1e6,
                      samplingSize = 1000,
                      plot =T)



calc_ratio(100,10, plot = F)

# BenchMark (Speed) -------------------------------------------------------


n = 1e6

test <- benchmark("empirical" = {estimate_pi_empirical(n)},
                  "resampled" = {estimate_pi_resampled(n)},
                  replications = 10)

test$meanTime <- test$elapsed/test$replications
test
 

# BenchMark (Accuracy) ----------------------------------------------------


# Accurcay Scores
test_accuracy(n = 1e5, nIter = 100, type = "empirical")
test_accuracy(n = 1e5, nIter = 100, type = "resampled")


# Accuracy of Estimate (raw)
mean_estimate(n = 1e6, nIter = 100, type = "empirical")
mean_estimate(n = 1e6, nIter = 100, type = "resampled")


# As expected the accuracies are more or less the same
