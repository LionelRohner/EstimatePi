
# Library -----------------------------------------------------------------

library(rbenchmark)

# Methods -----------------------------------------------------------------


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

dist <- generate_points(100000)

plot(dist[,1],dist[,2])

plot_pi_ratio <- function(dist)

# Main Function -----------------------------------------------------------

estimate_pi <- function(n){
  dist <- get_distance(generate_points(n = n))
  return(approx_pi(dist = dist))
}


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


# BenchMark ---------------------------------------------------------------

estimate_pi(1000000)
estimate_pi_v2(1000000)

n = 1e8

test <- benchmark("v1" = {estimate_pi(n)},
                  "v2" = {estimate_pi_v2(n)},
                  replications = 2)

test$meanTime <- test$elapsed/test$replications
test
 
