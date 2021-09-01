# EstimatePi
Estimate π from uniformly drawn points within a square of length 1. π is estimated as the ratio of points within a unit circle within the square divided to all points in the square.

## Empirical Method

*The Idea:* 
Create uniformly distributed points between 0 and 1 and consider them as vectors pointing to this random coordinate. Then count the vectors that have a vector norm of less than 1, i.e. the points inside the unit circle. The ratio of points inside to points outside the unit circle is an approximation to π/4. The accuracy of the estimate depends on the number of uniform points generated. Below 1e6, the accuracy is quite poor.

## Resampling Method

*The Idea:* 
The idea is to avoid generating a large number of points to get an accurate approximation of pi. Instead, we compute rather imprecise estimates of the ratio of the points inside and outside the unit circle and generate a beta or gamma distribution of ratios. The parameter estimate of the beta distribution is based on the method-of-moments to find initial values for alpha and beta in ~ Beta(alpha, beta). Once the distribution is fixed, we sample from that distribution and take the mean of the distribution, which should be a good estimate of pi.

This method has many parameters, such as the number of points to calculate a ratio, the number of initial ratios to generate, and the number of ratios to draw from the newly created distribution, and so some optimization is required (grid search or something more intelligent?). But at the moment, this method is neither more accurate nor faster than the empirical method...

## MCMC-like Method

MCMC-like algorithm, but instead of the Hastings ratio I used the accuracy measure, which kind of defeats the purpose of the MCMC algorithm, since it is used when the true value is unknown. But heck, it's just for fun, right? A newer implementation has the Hastings ratio with a uniform proposal kernel (cancels out).

*The Idea:*
1.) Create a prior distribution (beta distribution of unit circle ratio (same as in the resampling method))
2.) Propose first move (start with mean of prior)
3.) Initiate chain
4.) Propose a move with uniform proposal kernel
5.) Compute hastings-ratio H, accept if H >= than Uniform([0,1]). Alternatively, Compute accuracy, accept move if accuracy is better else stay.
6.) End chain at maxIter.

