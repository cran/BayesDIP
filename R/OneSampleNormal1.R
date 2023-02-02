# ------------------------------------------------------------------------------
#     prior --> list containing the information for prior
#          [[1]] - the prior distribution type:
#                1 - DIP
#                2 - Normal(mu0,var/n0)
#          [[2]] - n0: the level of information contained in the prior and
#                      the contribution of null mean
# ------------------------------------------------------------------------------

#' One sample Normal model with one-parameter unknown, given variance
#'
#' For a given planned sample size, the efficacy and futility boundaries,
#' return the power, the type I error, the expected sample size and its
#' standard deviation, the probability of reaching the efficacy and futility boundaries.
#'
#'
#' @param prior A list of length 2 containing the distributional information of the prior.
#' The first element is a number specifying the type of prior. Options are
#' \enumerate{
#' \item DIP ;
#' \item Normal(mu0,var/n0), where mu0 = prior mean, var = the known variance}
#' The second elements of the list is the parameter n0.
#' @param N The planned sample size.
#' @param mu0 The null mean value, which could be taken as the standard or current mean.
#' @param mu1 The mean value of the new treatment.
#' @param var The variance
#' @param d The target improvement (minimal clinically meaningful difference).
#' @param ps The efficacy boundary (upper boundary).
#' @param pf The futility boundary (lower boundary).
#' @param alternative less (lower values imply greater efficacy) or greater (larger
#' values imply greater efficacy).
#' @param seed The seed for simulations.
#' @param sim The number of simulations.
#' @return A list of the arguments with method and computed elements.
#' @examples
#' # with traditional Bayesian prior Beta(1,1)
#' OneSampleNormal1(list(2,6), N = 100, mu0 = 100, mu1 = 95, var=15, d = 0.05,
#'                   ps = 0.95, pf = 0.05, alternative = "less",
#'                   seed = 202210, sim = 10)
#' OneSampleNormal1(list(1,0), N = 100, mu0 = 100, mu1 = 95, var=15, d = 0.05,
#'                   ps = 0.95, pf = 0.05, alternative = "less",
#'                   seed = 202210, sim = 10)
#' @importFrom stats rbeta rbinom rgamma rnorm rpois
#' @export OneSampleNormal1

OneSampleNormal1 <- function(prior, N = 100, mu0, mu1, var, d = 0,
                               ps = 0.95, pf = 0.05,
                               alternative = c("less", "greater"), seed = 202209, sim = 5000) {

  alternative <- match.arg(alternative)
  # Define the inputs
  if(prior[[1]] == 1){
     prior[[2]] <- NA
  }
  ## N limit
  if(!is.null(N) && (!is.numeric(N) || N <= 0 ))
    stop("N must be positive number and greater than 10")
  ## mu0 limit
  if(!is.null(mu0) && (!is.numeric(mu0)))
    stop("mu0 must be numeric")
  ## mu1 limit
  if(!is.null(mu1) && (!is.numeric(mu1)))
    stop("mu1 must be numeric")
  ## var limit
  if(!is.null(var) && (!is.numeric(var) || var <  0))
    stop("var must be positive numeric")
  ## d limit
  if(!is.null(d) && (!is.numeric(d) || (d < 0 | d > abs(mu1-mu0))))
    stop("d must be numeric in [0, |mu1-mu0|]")
  ## efficacy boundary limit
  if(!is.null(ps) && (!is.numeric(ps) || (ps < 0.8 | ps > 1)))
    stop("ps (efficacy boundary) must be numeric in [0.8,1]")
  ## futility boundary limit
  if(!is.null(pf) && (!is.numeric(pf) || (pf < 0 | pf > 0.2)))
    stop("pf (futility boundary) must be numeric in [0,0.2]")
  ## set.seed
  if(!is.numeric(seed))
    stop("seed must be numeric")
  ## number of simulation
  if(!is.numeric(sim))
    stop("simulation number  must be numeric")

  set.seed(seed)

  # Functions to calculate the posterior
  Normal <- function(n0, mu0, var, y){
    n <- length(y)
    y_mean <- mean(y)
    posterior<-rnorm(1000,(n0*mu0+n*y_mean)/(n0+n), sqrt(var)/sqrt(n0+n))
  }
  Normal.DIP <- function(N, mu0, var, y){
    n <- length(y)
    y_mean <- mean(y)
    posterior<-rnorm(1000, ((N-n)*mu0+n*y_mean)/N, sqrt(var/N))
  }

  # Simulated Data
  # calculate power
  n.enrolled <- NULL
  cat1s <- 0
  cat1f <- 0
  for (k in 1:sim) {
    y.data<-NULL
    j<-0
    cat<-0
    cats<-0
    catf<-0
    pp_stop<-0.5
    while(cat == 0)
    {
      j<-j+1
      y.data<-append(y.data,rnorm(1,mu1,sqrt(var)))
      if(j>=10)
      {
        if (prior[[1]] == 2){
          mu1_s<-Normal(n0 = prior[[2]], mu0, var, y = y.data)
        }else if (prior[[1]] == 1){
          mu1_s<-Normal.DIP(N, mu0, var, y = y.data)
        }

        if (alternative == "greater"){
          pp_stop<-sum(mu1_s>mu0+d)/length(mu1_s)
        }else if (alternative == "less"){
          pp_stop<-sum(mu1_s<mu0-d)/length(mu1_s)
        }
      }
      if(pp_stop>=ps){cats<-1}
      if(pp_stop<pf){catf<-1}
      cat<-cats+catf
      if(j==N){cat<-1}
    }
    if(cats==1){cat1s<-cat1s+1}
    if(cats==0){cat1s<-cat1s}
    if(catf==1){cat1f<-cat1f+1}
    if(catf==0){cat1f<-cat1f}

    # Recruited Sample Size
    n.enrolled <- append(n.enrolled, j)
  }
  ss <- round(mean(n.enrolled), digits = 1)
  sd <- round(sd(n.enrolled), digits = 2)
  fut.rate <- cat1f/sim
  power <- cat1s/sim


  # calculate type I error
  cat1s <- 0
  cat1f <- 0
  for (k in 1:sim) {
    y.data<-NULL
    j<-0
    cat<-0
    cats<-0
    catf<-0
    pp_stop<-0.5
    while(cat == 0)
    {
      j<-j+1
      y.data<-append(y.data,rnorm(1,mu0, sqrt(var))) # under the null hypothesis mu1 = mu0
      if(j>=10)
      {
        if (prior[[1]] == 2){
          mu1_s<-Normal(n0 = prior[[2]], mu0, var, y = y.data)
        }else if (prior[[1]] == 1){
          mu1_s<-Normal.DIP(N, mu0, var, y = y.data)
        }

        if (alternative == "greater"){
          pp_stop<-sum(mu1_s>mu0+d)/length(mu1_s)
        }else if (alternative == "less"){
          pp_stop<-sum(mu1_s<mu0-d)/length(mu1_s)
        }
      }
      if(pp_stop>=ps){cats<-1}
      if(pp_stop<pf){catf<-1}
      cat<-cats+catf
      if(j==N){cat<-1}
    }
    if(cats==1){cat1s<-cat1s+1}
    if(cats==0){cat1s<-cat1s}
    if(catf==1){cat1f<-cat1f+1}
    if(catf==0){cat1f<-cat1f}
  }
  t1error <- cat1s/sim

  # Outputs
  if (prior[[1]] == 1) {method = "DIP"
  } else if (prior[[1]] == 2) {method = paste("Normal(",mu0, ",", var/prior[[2]], ")", sep="")}

  z <- list(method = method, power = power, type_I_error = t1error,
            expected_sample_size = ss, expected_sample_size_std = sd,
            the_prob_efficacy = power, the_prob_futility = fut.rate)
  z
}
