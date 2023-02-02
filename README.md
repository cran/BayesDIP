README
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesDIP

Provide early termination phase II trial designs with a decreasingly
informative prior (DIP) or a regular Bayesian prior chosen by the user.
The program can determine the minimum planned sample size necessary to
achieve the user-specified admissible designs. The program can also
perform power and expected sample size calculations for the tests in
early termination Phase II trials.\[1\]

## Installation

You can install from CRAN with:

``` r
install.packages("BayesDIP")
```

Or try the development version from \[GitHub\] with:

``` r
# install.packages("devtools")
devtools::install_github("chenw10/BayesDIP")
```

## Example

``` r
library(BayesDIP)

# Calculate the minimum planned sample size within the range 10<=N<=100,
# under an admissible design which is set as 80% power and 5% type I error here.

# One sample Bernoulli model with the response rate for the new treatment is 0.5, 
# the null response rate is 0.3, and the target improvement to achieve is 0.
# The alternative hypothesis: p1 > p0 + d
# Simulate 10 replicate trials using this design with efficacy boundary 0.98 
# and futility boundary 0.05.

### Designs with traditional Bayesian prior Beta(1,1)
### Designs and operating characteristics based on 100 simulations:
OneSampleBernoulli.Design(list(2,1,1), nmin = 10, nmax=100, p0 = 0.3, p1 = 0.5, d = 0,
                          ps = 0.98, pf = 0.02, power = 0.80, t1error=0.05, alternative = "greater",
                          seed = 202210, sim = 100)
#> 
#> Prior:   Beta(1,1) 
#> Planned Sample Size:   92 
#> Efficacy Boundary:   0.98 
#> Futility Boundary:   0.02 
#> Exact Power:  0.99 
#> Exact Type I error:   0.05 
#> Expected sample size:  24 
#> Expected sample size standard deviation:  17.6


### Designs with DIP
### Designs and operating characteristics based on 10 simulations:
OneSampleBernoulli.Design(list(1,0,0), nmin = 10, nmax=100, p0 = 0.3, p1 = 0.5, d = 0,
                          ps = 0.98, pf = 0.02, power = 0.80, t1error=0.05, alternative = "greater",
                          seed = 202210, sim = 100)
#> 
#> Prior:   DIP 
#> Planned Sample Size:   44 
#> Efficacy Boundary:   0.98 
#> Futility Boundary:   0.02 
#> Exact Power:  0.81 
#> Exact Type I error:   0.05 
#> Expected sample size:  30 
#> Expected sample size standard deviation:  9.7


# Calculate the power, type I error and the expected sample size given a planned sample size

# One sample Bernoulli model with the response rate for the new treatment is 0.5, 
# the null response rate is 0.3, and the target improvement to achieve is 0.05.
# The alternative hypothesis: p1 > p0 + d
# Simulate 100 replicate trials for a given planned sample size 100 using this design
# with efficacy boundary 0.98 and futility boundary 0.05.  

## with traditional Bayesian prior Beta(1,1)
## Operating characteristics based on 100 simulations:
OneSampleBernoulli(list(2,1,1), N = 100, p0 = 0.3, p1 = 0.5, d = 0.05,
                   ps = 0.98, pf = 0.05, alternative = "greater",
                   seed = 202210, sim = 100)
#> 
#> Prior:   Beta(1,1) 
#> Power:  0.89 
#> Type I error:   0.05 
#> Expected sample size:   42.7 
#> Expected sample size standard deviation:   29.06 
#> The probability of reaching the efficacy boundary:   0.89 
#> The probability of reaching the futility boundary:   0.02


## with DIP
## Operating characteristics based on 100 simulations:
OneSampleBernoulli(list(1,0,0), N = 100, p0 = 0.3, p1 = 0.5, d = 0.05,
                   ps = 0.98, pf = 0.05, alternative = "greater",
                   seed = 202210, sim = 100)
#> 
#> Prior:   DIP 
#> Power:  0.86 
#> Type I error:   0.01 
#> Expected sample size:   72.1 
#> Expected sample size standard deviation:   17.28 
#> The probability of reaching the efficacy boundary:   0.86 
#> The probability of reaching the futility boundary:   0
```

## Reference

\[1\] Wang C, Sabo RT, Mukhopadhyay ND, and Perera RA. Early termination
in single-parameter model phase II clinical trial designs using
decreasingly informative priors. , 9(2): April - June 2022.
<https://doi.org/10.18203/2349-3259.ijct20221110>
