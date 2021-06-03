
# Methods -----------------------------------------------------------------


generate_points_distance <- function(n){
  
  # # generate x,y coordinates
  # x <- runif(n)
  # y <- runif(n)
  
  # put them in a matrix
  XY <- matrix(runif(2*n), nrow = n, ncol = 2)
  
  # calculate vector norm and test if its bigger than the radius (i.e. 1)
  dist <- sqrt(XY[,1]^2 + XY[,2]^2) < 1
  return(dist)
}



generate_points_distance <- function(n){
  
  # # generate x,y coordinates
  # x <- runif(n)
  # y <- runif(n)
  
  # put them in a matrix
  XY <- matrix(runif(2*n), nrow = n, ncol = 2)
  
  # calculate vector norm and test if its bigger than the radius (i.e. 1)
  dist <- sqrt(XY[,1]^2 + XY[,2]^2) < 1
  return(dist)
}


approx_pi <- function(dist){
  # pi*r^2 / 4 == numb_circle / numb_total (divided by for because we have just one quadrant)
  return(4*sum(dist)/length(dist))
}


# Main Function -----------------------------------------------------------

estimate_pi <- function(n){
  dist <- generate_points_distance(n = n)
  return(approx_pi(dist = dist))
}

estimate_pi(10000000)
