# The paper uses X = (X1, X2), tau(X) for the treatment effect, and c(X)
# for confounding. U affects treatment assignment but is never observed.
H <- function(x1, x2) 1 + .35*x1 + .15*x1*x1 + .25*x2*x2

paper_setting <- function(setting = "S4") {
  switch(setting,
    S1 = list(tau = function(x1, x2) (x1-.25)*H(x1,x2),
              c = function(x1,x2) .40*x1+.10*x2,
              eta0 = function(x1) 1, eta1 = 1.25, cutoff = .25, degree = 3L,
              z = function(x1,x2) cbind(x1,x2)),
    S3 = list(tau = function(x1,x2) (x1-.25)*(1+.20*x2),
              c = function(x1,x2) .35+.25*x1-.25*x2,
              eta0 = function(x1) sqrt(.49+2*x1^2), eta1 = .25, cutoff = .25,
              degree = 2L, z = function(x1,x2) polynomial_basis(x1,x2,1L)),
    # Same tau, c, cutoff, bases and target outcomes in S2 and S4.
    # Only the RCT outcome noise changes, not the covariate distribution.
    S2 = , S4 = list(tau = function(x1,x2) x1*(H(x1,x2)*(1+.10*x1)),
                     c = function(x1,x2) -1.65*x1^2,
                     eta0 = function(x1) .04+.65*x1*x1,
                     eta1 = if (setting == "S2") 1.25 else .625,
                     cutoff = 0, degree = 4L, z = function(x1,x2) cbind(x1^2)),
    stop("setting must be S1, S2, S3 or S4"))
}

# Each row changes at most one feature of the baseline.
conditions <- data.frame(
  condition = c("baseline","half_target","double_target","double_trial","contraction","expansion"),
  n0 = c(4000,2000,8000,4000,4000,4000),
  n1 = c(100,100,100,200,100,100),
  sigma1 = c(1,1,1,1,.8,1.25))

generate_data <- function(setting = "S4", condition = "baseline", rep = 1L) {
  d <- paper_setting(setting)
  row <- match(condition, conditions$condition)
  stopifnot(!is.na(row), length(rep)==1L, rep %in% 1:25)
  design <- conditions[row, ]
  seed <- 9730000L + rep  # same streams as the published S1--S3 run

  draw_population <- function(n, s, offset) {
    set.seed(seed + offset)
    sigma <- if (s == 0) 1 else design$sigma1
    draw_x <- function() sigma*qnorm(runif(n,pnorm(-3/sigma),pnorm(3/sigma)))
    x1 <- draw_x()
    x2 <- draw_x()
    U <- sample(c(-1,1), n, replace = TRUE)
    A <- rbinom(n, 1, .5 + (1-s)*U/4)
    eta <- if (s == 0) d$eta0(x1) else d$eta1
    tau <- d$tau(x1,x2)
    c_x <- d$c(x1,x2)
    Y0 <- -.5*tau + (1-s)*U*c_x + eta*rnorm(n)
    Y1 <-  .5*tau + (1-s)*U*c_x + eta*rnorm(n)
    data.frame(x1,x2,A,Y=ifelse(A==1,Y1,Y0),Y0,Y1)
  }
  observed <- c("x1","x2","A","Y")
  list(target = draw_population(design$n0,0,200000L)[observed],
       trial  = draw_population(design$n1,1,300000L)[observed],
       test   = draw_population(50000L,0,100000L),
       setting = setting, condition = condition, rep = rep, seed = seed)
}

# Orthonormal Hermite products span the bivariate polynomial working model.
# They are normalized under N(0,I), not necessarily under the truncated design.
polynomial_basis <- function(x1,x2,degree) {
  hermite <- function(x) {
    a <- matrix(1,length(x),degree+1L)
    if (degree>=1L) a[,2] <- x
    if (degree>=2L) for (j in 2:degree)
      a[,j+1L] <- (x*a[,j]-sqrt(j-1)*a[,j-1L])/sqrt(j)
    a
  }
  a <- hermite(x1); b <- hermite(x2)
  ij <- do.call(rbind,lapply(0:degree,function(k)cbind(k:0,0:k)))
  basis <- vapply(seq_len(nrow(ij)),function(k)
    a[,ij[k,1]+1L]*b[,ij[k,2]+1L],numeric(length(x1)))
  basis <- matrix(basis,nrow=length(x1))
  colnames(basis) <- apply(ij,1,function(k)paste0("H",k[1],"_",k[2]))
  basis
}

# Evaluation only. A threshold can disagree with the optimal rule only on
# intervals separated by its root and the true cutoff. On each interval
# tau has constant sign. Integrate over X2 first (E[X2]=0, E[X2^2]=m2).
population_regret <- function(theta, setting) {
  d <- paper_setting(setting)
  z <- pnorm(3)-pnorm(-3)
  m2 <- 1-6*dnorm(3)/z
  x2 <- if (setting=="S3") 0 else sqrt(m2)
  breaks <- c(-3,d$cutoff,3)
  if (theta[2]!=0) breaks <- c(breaks,min(3,max(-3,-theta[1]/theta[2])))
  breaks <- sort(unique(breaks))
  sum(vapply(seq_len(length(breaks)-1L), function(j) {
    lo <- breaks[j]; hi <- breaks[j+1L]; mid <- (lo+hi)/2
    if ((mid>d$cutoff)==(theta[1]+theta[2]*mid>0)) return(0)
    integrate(function(x) abs(d$tau(x,x2))*dnorm(x)/z,lo,hi,rel.tol=1e-10)$value
  },numeric(1)))
}
