# ------------------------------------------------------------------------------
#     prior --> list containing the information for prior
#          [[1]] - the prior distribution type:
#                1 - DIP
#                2 - Gamma(a,b)
#          [[2]] - a: shape parameter of the Gamma distribution
#          [[3]] - b: rate parameter of the Gamma distribution
# ------------------------------------------------------------------------------

#' One sample Poisson model - Trial Design
#'
#' Calculate the minimum planned sample size under an admissible design.
#' The users decide the power and type-I-error, and pick the efficacy and futility boundaries.
#' If there are no admissible design based on controlled type-I-error, then default to output
#' the designs with the lowest type-I-error and at least the user-defined (e.g. 80\%) power.
#'
#'
#' @param prior A list of length 3 containing the distributional information of the prior.
#' The first element is a number specifying the type of prior. Options are
#' \enumerate{
#' \item DIP ;
#' \item Gamma(a,b), where a = shape, b = rate}
#' The second and third elements of the list are the parameters a and b, respectively.
#' @param nmin The start searching sample size
#' @param nmax The stop searching sample size
#' @param m0 The null event rate, which could be taken as the standard or current event rate.
#' @param m1 The event rate of the new treatment.
#' @param d The target improvement (minimal clinically meaningful difference).
#' @param ps The efficacy boundary (upper boundary).
#' @param pf The futility boundary (lower boundary).
#' @param power The expected power to achieve.
#' @param t1error The controlled type-I-error.
#' @param alternative less (lower values imply greater efficacy) or greater (larger
#' values imply greater efficacy).
#' @param seed The seed for simulations.
#' @param sim The number of simulations.
#' @return A list of the arguments with method and computed elements
#' @examples
#' \donttest{
#' # with traditional Bayesian prior Gamma(0.5,0.001)
#' OneSamplePoisson.Design(list(2,0.5,0.001), nmin = 10, nmax=100, m0 = 5, m1 = 4, d = 0,
#'                    ps = 0.95, pf = 0.05, power = 0.80, t1error=0.05, alternative = "less",
#'                    seed = 202210, sim = 10)
#' # with DIP
#' OneSamplePoisson.Design(list(1,0,0), nmin = 10, nmax=100, m0 = 5, m1 = 4, d = 0,
#'                    ps = 0.95, pf = 0.05, power = 0.80, t1error=0.05, alternative = "less",
#'                    seed = 202210, sim = 10)
#' }
#' @importFrom stats rbeta rbinom rgamma rnorm rpois
#' @export OneSamplePoisson.Design


OneSamplePoisson.Design <- function(prior, nmin = 10, nmax = 100, m0, m1, d = 0,
                                      ps, pf, power = 0.8, t1error = 0.05,
                                      alternative = c("less", "greater"), seed = 202209, sim = 1000){

  alternative <- match.arg(alternative)
  # Define the inputs
  if(prior[[1]] == 1){
    prior[[2]] <- NA
    prior[[3]] <- NA
  }
  ## nmin limit
  if(!is.null(nmin) && (!is.numeric(nmin) || nmin < 10 || nmin >= nmax))
    stop("nmin must be positive number and at least 10")
  ## nmax limit
  if(!is.null(nmax) && (!is.numeric(nmax) || nmax <= nmin || nmax >= 200))
    stop("nmax must greater than 'nmin' and less than 200")
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
  ## power limit
  if(!is.null(power) && (!is.numeric(power) || (power < 0 | power > 1)))
    stop("power must be numeric in [0,1]")
  ## t1error limit
  if(!is.null(t1error) && (!is.numeric(t1error) || (t1error < 0 |t1error > 1)))
    stop("type-I-error must be numeric in [0,1]")
  ## set.seed
  if(!is.numeric(seed))
    stop("seed must be numeric")
  if(!is.numeric(sim))
    stop("simulation number must be numeric")

  set.seed(seed)

  # Functions to calculate the posterior
  Poisson <- function(a,b,y){posterior<-rgamma(1000, a+sum(y), b+length(y))}
  Poisson.DIP <- function(m0, y, N){
    posterior<-rgamma(1000, 0.05+sum(y)+m0*(N-length(y)), 0.001+N)
  }


  # Simulated Data
  # calculate N that can achieve the power
  N_v <- NULL
  power_v <- NULL
  n_v <- NULL
  sd_v <- NULL
  for (N in seq(from=nmin, to=nmax, by=1)){
    cat1s <- 0
    cat1f <- 0
    n.enrolled <- NULL
    for (k in 1:sim) {
      y.data <- NULL
      j <- 0
      cat <- 0
      cats <- 0
      catf <- 0
      pp_stop <- 0.5
      while(cat == 0){
        j <- j+1
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
    power.cal <- cat1s/sim

    jitter <- 0.01
    if (power.cal >= power-jitter){
      N_v <- append(N_v, N)
      power_v <- append(power_v, power.cal)
      n_v <- append(n_v, round(mean(n.enrolled), 0))
      sd_v <- append(sd_v, round(sd(n.enrolled), 1))
    }
    result1 <- cbind(N_v, power_v, n_v, sd_v)
  } # End of power calculation


  if (is.null(result1)){
    message("Suggest: please adjust your input values!")
    stop(paste("No sample size in the range [",nmin,",",nmax,"] can achieve ", power*100, "% power", sep=""))
  }


  # calculate type I error
  nmin1 <- N_v[which.min(N_v)] # start minimum sample size in calculation of exact type I error
  N_v <- NULL
  t1error_v <- NULL
  for (N in seq(from=nmin1, to=nmax, by=1)){
    cat1s <- 0
    cat1f <- 0
    for (k in 1:sim) {
      y.data <- NULL
      j <- 0
      cat <- 0
      cats <- 0
      catf <- 0
      pp_stop <- 0.5
      while(cat == 0){
        j <- j+1
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
    t1error.cal <- cat1s/sim
    N_v <- append(N_v, N)
    t1error_v <- append(t1error_v, t1error.cal)
    result2 <- cbind(N_v, t1error_v)

  } # End of Type-I-error calculation


  # Outputs
  if (!is.null(result1) & !is.null(result2)){
    result <- merge(result1, result2, by=c("N_v"))
    final <- as.data.frame(result)

    # select the lowest/best-controlled type I error
    final$diff <- abs(final$t1error_v - t1error)
    final <- final[order(final$diff, final$t1error_v, final$power_v, final$N_v), ]
    ff <- final[1,]
    planN <- ff$N_v
    exact.power <- ff$power_v
    exact.t1 <- ff$t1error_v
    ss <- ff$n_v
    sd <- ff$sd_v
    if (prior[[1]] == 1) {method = "DIP"
    } else if (prior[[1]] == 2) {method = paste("Gamma(",prior[[2]], ",", prior[[3]], ")", sep="")}

    z <- list(method = method, planned_sample_size = planN,
              efficacy_boundary = ps, futility_boundary = pf,
              exact_power = exact.power, exact_type_I_error = exact.t1,
              expected_sample_size = ss, expected_sample_size_std = sd)
    z
  } # End of Outputs

}


