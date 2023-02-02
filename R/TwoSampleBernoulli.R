# ------------------------------------------------------------------------------
#     prior --> list containing the information for prior
#          [[1]] - the prior distribution type:
#                1 - DIP
#                2 - Beta(a,b)
#          [[2]] - a: first parameter of the Beta distribution
#          [[3]] - b: second parameter of the Beta distribution
# ------------------------------------------------------------------------------

#' Two sample Bernoulli model
#'
#' For a given planned sample size, the efficacy and futility boundaries,
#' return the power, the type I error, the expected sample size and its
#' standard deviation, the probability of reaching the efficacy and futility boundaries.
#' Equal allocation between two treatment groups.
#'
#'
#' @param prior A list of length 3 containing the distributional information of the prior.
#' The first element is a number specifying the type of prior. Options are
#' \enumerate{
#' \item DIP ;
#' \item Beta(a,b), where a = shape, b = scale}
#' The second and third elements of the list are the parameters a and b, respectively.
#' @param N The total planned sample size for two treatment groups.
#' @param p1 The response rate of the new treatment.
#' @param p2 The response rate of the compared treatment.
#' @param d The target improvement (minimal clinically meaningful difference).
#' @param ps The efficacy boundary (upper boundary).
#' @param pf The futility boundary (lower boundary).
#' @param alternative less (lower values imply greater efficacy) or greater (larger
#' values imply greater efficacy).
#' @param seed The seed for simulations.
#' @param sim The number of simulations.
#' @return A list of the arguments with method and computed elements
#' @examples
#' # with traditional Bayesian prior Beta(1,1)
#' TwoSampleBernoulli(list(2,1,1), N = 200, p1 = 0.5, p2 = 0.3, d = 0,
#'                    ps = 0.90, pf = 0.05, alternative = "greater",
#'                    seed = 202210, sim = 5)
#' # with DIP
#' TwoSampleBernoulli(list(1,0,0), N = 200, p1 = 0.5, p2 = 0.3, d = 0,
#'                    ps = 0.90, pf = 0.05, alternative = "greater",
#'                    seed = 202210, sim = 5)
#' @import stats
#' @export TwoSampleBernoulli


TwoSampleBernoulli <- function(prior, N = 200, p1, p2, d = 0,
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
  # calculate power
  n.enrolled <- NULL
  n1.enrolled <- NULL
  n2.enrolled <- NULL
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
          y.data<-append(y.data,rbinom(1,1,p1))}
        if(u>r){
          Group=append(Group,0)
          y.data<-append(y.data,rbinom(1,1,p2))}
        Matd<-as.data.frame(cbind(y.data,Group))
        y1<-Matd$y[which(Matd$Group==1)]
        y2<-Matd$y[which(Matd$Group==0)]
        sn1 <- length(y1[y1==1]) # number of successes
        sn2 <- length(y2[y2==1])

        if(sn1>0 & sn2>0)
        {
          if (prior[[1]] == 2){
            p1s<-Bernoulli(a = prior[[2]], b = prior[[3]], y = y1)
            p2s<-Bernoulli(a = prior[[2]], b = prior[[3]], y = y2)
          }else if (prior[[1]] == 1){
            N1<-ceiling(N/2)
            N2<-ceiling(N/2)
            p0<-rbeta(1000,1,1) # hyper-prior
            p1s<-Bernoulli.DIP(p0, y1, N1)
            p2s<-Bernoulli.DIP(p0, y2, N2)
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
ss1 <- round(mean(n1.enrolled), digits = 0)
ss2 <- round(mean(n2.enrolled), digits = 0)
std1 <- round(sd(n1.enrolled), digits = 1)
std2 <- round(sd(n2.enrolled), digits = 1)
fut.rate <- cat1f/sim
power <- cat1s/sim


# calculate type I error
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

        if(sn1>0 & sn2>0)
        {
          if (prior[[1]] == 2){
            p1s<-Bernoulli(a = prior[[2]], b = prior[[3]], y = y1)
            p2s<-Bernoulli(a = prior[[2]], b = prior[[3]], y = y2)
          }else if (prior[[1]] == 1){
            N1<-ceiling(N/2)
            N2<-ceiling(N/2)
            p0<-rbeta(1000,1,1) # hyper-prior
            p1s<-Bernoulli.DIP(p0, y1, N1)
            p2s<-Bernoulli.DIP(p0, y2, N2)
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
t1error <- cat1s/sim

  # Outputs
  if (prior[[1]] == 1) {method = "DIP"
  } else if (prior[[1]] == 2) {method = paste("Beta(",prior[[2]], ",", prior[[3]], ")", sep="")}

  grp1 = paste(ss1, " (", std1, ")", sep="")
  grp2 = paste(ss2, " (", std2, ")", sep="")

  z <- list(method = method, power = power, type_I_error = t1error,
            The_planned_sample_size_for_each_group = N/2,
            expected_sample_size_and_std_for_the_new_treatment = grp1,
            expected_sample_size_and_std_for_the_compared_treatment = grp2,
            the_prob_efficacy = power, the_prob_futility = fut.rate)
  z
}

