# ------------------------------------------------------------------------------
#     prior --> list containing the information for prior
#          [[1]] - the prior distribution type:
#                1 - DIP
#                2 - Gamma(a,b)
#          [[2]] - a: shape parameter of the Gamma distribution
#          [[3]] - b: rate parameter of the Gamma distribution
# ------------------------------------------------------------------------------

#' One sample Poisson model
#'
#' For a given planned sample size, the efficacy and futility boundaries,
#' return the power, the type I error, the expected sample size and its
#' standard deviation, the probability of reaching the efficacy and futility boundaries.
#'
#'
#' @param prior A list of length 3 containing the distributional information of the prior.
#' The first element is a number specifying the type of prior. Options are
#' \enumerate{
#' \item DIP ;
#' \item Gamma(a,b), where a = shape, b = rate}
#' The second and third elements of the list are the parameters a and b, respectively.
#' @param N The planned sample size.
#' @param m0 The null event rate, which could be taken as the standard or current event rate.
#' @param m1 The event rate of the new treatment.
#' @param d The target improvement (minimal clinically meaningful difference).
#' @param ps The efficacy boundary (upper boundary).
#' @param pf The futility boundary (lower boundary).
#' @param alternative less (lower values imply greater efficacy) or greater (larger
#' values imply greater efficacy).
#' @param seed The seed for simulations.
#' @param sim The number of simulations.
#' @return A list of the arguments with method and computed elements
#' @examples
#' # with traditional Bayesian prior Gamma(0.5,0.001)
#' OneSamplePoisson(list(2,0.5,0.001), N = 100, m0 = 0.5, m1 = 0.4, d = 0.05,
#'                  ps = 0.95, pf = 0.05, alternative = "less",
#'                  seed = 202210, sim = 10)
#' # with DIP
#' OneSamplePoisson(list(1,0,0), N = 100, m0 = 0.5, m1 = 0.4, d = 0.05,
#'                  ps = 0.95, pf = 0.05, alternative = "less",
#'                  seed = 202210, sim = 10)
#' @importFrom stats rbeta rbinom rgamma rnorm rpois
#' @export OneSamplePoisson

OneSamplePoisson <- function(prior, N = 100, m0, m1, d = 0,
                              ps = 0.95, pf = 0.05,
                              alternative = c("less", "greater"), seed = 202209, sim = 5000) {

  alternative <- match.arg(alternative)
  # Define the inputs
  if(prior[[1]] == 1){
     prior[[2]] <- NA
     prior[[3]] <- NA
  }
  ## N limit
  if(!is.null(N) && (!is.numeric(N) || N <= 0 ))
    stop("N must be positive number and greater than 10")
  ## m0 limit
  if(!is.null(m0) && (!is.numeric(m0) || (m0 < 0)))
    stop("m0 must be numeric in [0,inf]")
  ## m1 limit
  if(!is.null(m1) && (!is.numeric(m1) || (m1 < 0)))
    stop("m1 must be numeric in [0,inf]")
  ## d limit
  if(!is.null(d) && (!is.numeric(d) || (d < 0 | d > abs(m1-m0))))
    stop("d must be numeric in [0, |m1-m0|]")
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
    stop("simulation number must be numeric")

  set.seed(seed)

  # Functions to calculate the posterior
  Poisson <- function(a,b,y){posterior<-rgamma(1000, a+sum(y), b+length(y))}
  Poisson.DIP <- function(m0, y, N){
    posterior<-rgamma(1000, 0.05+sum(y)+m0*(N-length(y)), 0.001+N)
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
      y.data<-append(y.data,rpois(1,m1))
      if(j>=10)
      {
        if (prior[[1]] == 2){
          m1_s<-Poisson(a = prior[[2]], b = prior[[3]], y = y.data)
        }else if (prior[[1]] == 1){
          m1_s <- Poisson.DIP(m0, y = y.data, N = N)
        }

        if (alternative == "greater"){
          pp_stop<-sum(m1_s>m0+d)/length(m1_s)
        }else if (alternative == "less"){
          pp_stop<-sum(m1_s<m0-d)/length(m1_s)
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
      y.data<-append(y.data,rpois(1,m0)) # under the null hypothesis m1 = m0
      if(j>=10)
      {
        if (prior[[1]] == 2){
          m1_s<-Poisson(a = prior[[2]], b = prior[[3]], y = y.data)
        }else if (prior[[1]] == 1){
          m1_s <- Poisson.DIP(m0, y = y.data, N = N)
        }

        if (alternative == "greater"){
          pp_stop<-sum(m1_s>m0+d)/length(m1_s)
        }else if (alternative == "less"){
          pp_stop<-sum(m1_s<m0-d)/length(m1_s)
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
  } else if (prior[[1]] == 2) {method = paste("Gamma(",prior[[2]], ",", prior[[3]], ")", sep="")}

  z <- list(method = method, power = power, type_I_error = t1error,
            expected_sample_size = ss, expected_sample_size_std = sd,
            the_prob_efficacy = power, the_prob_futility = fut.rate)
  z
}
