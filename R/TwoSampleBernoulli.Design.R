# ------------------------------------------------------------------------------
#     prior --> list containing the information for prior
#          [[1]] - the prior distribution type:
#                1 - DIP
#                2 - Beta(a,b)
#          [[2]] - a: first parameter of the Beta distribution
#          [[3]] - b: second parameter of the Beta distribution
# ------------------------------------------------------------------------------

#' Two sample Bernoulli model - Trial Design
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
#' \item Beta(a,b), where a = shape, b = scale}
#' The second and third elements of the list are the parameters a and b, respectively.
#' @param nmin The start searching total sample size for two treatment groups.
#' @param nmax The stop searching total sample size for two treatment groups.
#' @param p1 The response rate of the new treatment.
#' @param p2 The response rate of the compared treatment.
#' @param d The target improvement (minimal clinically meaningful difference).
#' @param ps The efficacy boundary (upper boundary).
#' @param pf The futility boundary (lower boundary).
#' @param power The power to achieve.
#' @param t1error The controlled type-I-error.
#' @param alternative less (lower values imply greater efficacy) or greater (larger
#' values imply greater efficacy).
#' @param seed The seed for simulations.
#' @param sim The number of simulations.
#' @return A list of the arguments with method and computed elements
#' @examples
#' \donttest{
#' # with traditional Bayesian prior Beta(1,1)
#' TwoSampleBernoulli.Design(list(2,1,1), nmin = 100, nmax = 120, p1 = 0.5, p2 = 0.3, d = 0,
#'                    ps = 0.90, pf = 0.05, power = 0.8, t1error = 0.05, alternative = "greater",
#'                    seed = 202210, sim = 10)
#' # with DIP
#' TwoSampleBernoulli.Design(list(1,0,0), nmin = 100, nmax = 120, p1 = 0.5, p2 = 0.3, d = 0,
#'                    ps = 0.90, pf = 0.05, power = 0.8, t1error = 0.05, alternative = "greater",
#'                    seed = 202210, sim = 10)
#' }
#' @import stats
#' @export TwoSampleBernoulli.Design


TwoSampleBernoulli.Design <- function(prior, nmin=10, nmax = 200, p1, p2, d = 0,
                               ps = 0.95, pf = 0.05, power = 0.80, t1error = 0.05,
                               alternative = c("less", "greater"), seed = 202209, sim = 500){

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
  if(!is.null(nmax) && (!is.numeric(nmax) || nmax <= nmin || nmax >= 300))
    stop("nmax must greater than 'nmin' and less than 300")
  ## p1 limit
  if(!is.null(p1) && (!is.numeric(p1) || (p1 < 0 | p1 > 1)))
    stop("p1 must be numeric in [0,1]")
  ## p2 limit
  if(!is.null(p2) && (!is.numeric(p2) || (p2 < 0 | p2 > 1)))
    stop("p2 must be numeric in [0,1]")
  ## d limit
  if(!is.null(d) && (!is.numeric(d) || (d < 0 | d > abs(p1-p2))))
    stop("d must be numeric in [0, |p1-p2|]")
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
  ## number of simulation
  if(!is.numeric(sim))
    stop("simulation number  must be numeric")

  set.seed(seed)

  # Functions to calculate the posterior
  Bernoulli <- function(a,b,y){posterior<-rbeta(1000, a+sum(y), b+(length(y)-sum(y)))}
  Bernoulli.DIP <- function(p0, y, N){
    j<-length(y)
    posterior<-rbeta(1000,1+sum(y)+p0*(N-j),1+(j-sum(y))+(1-p0)*(N-j))
  }

  # Simulated Data
  # calculate N that can achieve the power
  N_v <- NULL
  power_v <- NULL
  n1_v <- NULL
  n2_v <- NULL
  std1_v <- NULL
  std2_v <- NULL
for (N in seq(from=nmin, to=nmax, by=1)){
    cat1s <- 0
    cat1f <- 0
    n.enrolled <- NULL
    n1.enrolled <- NULL
    n2.enrolled <- NULL
  for (k in 1:sim) {
    y.data<-NULL
    Group<-NULL
    j<-0
    r<-0.5 # equal allocation
    cat<-0
    cats<-0
    catf<-0
    pp_stop<-0.5
    while(cat == 0)
    {
      j<-j+1
      u<-runif(1,min = 0,max = 1)
      if(u<=r){
        Group=append(Group,1)
        y.data<-append(y.data,rbinom(1,1,p1))}
      if(u>r){
        Group=append(Group,0)
        y.data<-append(y.data,rbinom(1,1,p2))}
      Matd<-as.data.frame(cbind(y.data,Group))
      y1<-Matd$y[which(Matd$Group==1)]
      y2<-Matd$y[which(Matd$Group==0)]
      sn1 <- length(y1[y1==1]) # number of successes
      sn2 <- length(y2[y2==1])

      if(j>=10 & sn1>0 & sn2>0)
      {
        if (prior[[1]] == 2){
          p1s<-Bernoulli(a = prior[[2]], b = prior[[3]], y = y1)
          p2s<-Bernoulli(a = prior[[2]], b = prior[[3]], y = y2)
        }else if (prior[[1]] == 1){
          N1<-ceiling(N/2)
          N2<-ceiling(N/2)
          p0<-rbeta(1000,1,1) # hyper-prior
          p1s<-Bernoulli.DIP(p0, y=y1, N=N1)
          p2s<-Bernoulli.DIP(p0, y=y2, N=N2)
          p1s[is.na(p1s)]<-sum(y1)/length(y1)
          p2s[is.na(p2s)]<-sum(y2)/length(y2)
        }

        if (alternative == "greater"){
          pp_stop<-sum(p1s>p2s+d)/length(p1s)
        }else if (alternative == "less"){
          pp_stop<-sum(p1s<p2s-d)/length(p1s)
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
    n1.enrolled <- append(n1.enrolled, length(y1))
    n2.enrolled <- append(n2.enrolled, length(y2))
  }
  power.cal <- cat1s/sim

  jitter <- 0.01
  if (power.cal >= power-jitter){
    N_v <- append(N_v, N)
    power_v <- append(power_v, power.cal)
    n1_v <- append(n1_v, round(mean(n1.enrolled), 0))
    n2_v <- append(n2_v, round(mean(n2.enrolled), 0))
    std1_v <- append(std1_v, round(sd(n1.enrolled), 1))
    std2_v <- append(std2_v, round(sd(n2.enrolled), 1))
  }
  result1 <- cbind(N_v, power_v, n1_v, std1_v, n2_v, std2_v)
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
    y.data<-NULL
    Group<-NULL
    j<-0
    r<-0.5 # equal allocation
    cat<-0
    cats<-0
    catf<-0
    pp_stop<-0.5
    while(cat == 0)
    {
      j<-j+1
      u<-runif(1,min = 0,max = 1)
      if(u<=r){
        Group=append(Group,1)
        y.data<-append(y.data,rbinom(1,1,p2))}  #under the null hypothesis p1 = p2
      if(u>r){
        Group=append(Group,0)
        y.data<-append(y.data,rbinom(1,1,p2))}
      Matd<-as.data.frame(cbind(y.data,Group))
      y1<-Matd$y[which(Matd$Group==1)]
      y2<-Matd$y[which(Matd$Group==0)]
      sn1 <- length(y1[y1==1]) # number of successes
      sn2 <- length(y2[y2==1])

      if(j>=10 & sn1>0 & sn2>0)
      {
        if (prior[[1]] == 2){
          p1s<-Bernoulli(a = prior[[2]], b = prior[[3]], y = y1)
          p2s<-Bernoulli(a = prior[[2]], b = prior[[3]], y = y2)
        }else if (prior[[1]] == 1){
          N1<-ceiling(N/2)
          N2<-ceiling(N/2)
          p0<-rbeta(1000,1,1) # hyper-prior
          p1s<-Bernoulli.DIP(p0, y=y1, N=N1)
          p2s<-Bernoulli.DIP(p0, y=y2, N=N2)
          p1s[is.na(p1s)]<-sum(y1)/length(y1)
          p2s[is.na(p2s)]<-sum(y2)/length(y2)
        }

        if (alternative == "greater"){
          pp_stop<-sum(p1s>p2s+d)/length(p1s)
        }else if (alternative == "less"){
          pp_stop<-sum(p1s<p2s-d)/length(p1s)
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
    grp1 = paste(ff$n1_v, " (", ff$std1_v, ")", sep="")
    grp2 = paste(ff$n2_v, " (", ff$std2_v, ")", sep="")
    if (prior[[1]] == 1) {method = "DIP"
    } else if (prior[[1]] == 2) {method = paste("Beta(",prior[[2]], ",", prior[[3]], ")", sep="")}


    z <- list(method = method, planned_sample_size = planN,
              efficacy_boundary = ps, futility_boundary = pf,
              exact_power = exact.power, exact_type_I_error = exact.t1,
              expected_sample_size_and_std_for_the_new_treatment_group = grp1,
              expected_sample_size_and_std_for_the_compared_treatment_group = grp2)
    z
  } # End of Outputs
}


