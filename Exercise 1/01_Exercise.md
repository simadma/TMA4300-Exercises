Exercise 1
================
Mads Adrian Simonsen, William Scott Grundeland Olsen
08 februar, 2021

``` r
#remove this------------------------------------------------------------------------------
generate_from_exp <- function(n, rate = 1) {
  Y <- runif(n)
  X <- -(1 / rate) * log(Y)
  return(X)
}
std_normal <- function(n) {
  X1 <- pi * runif(n)   # n samples from Uniform(0, pi)
  X2 <- generate_from_exp(n, 1/2)   # n samples from Exponential(1/2)
  Z <- X2^(1/2) * cos(X1)   # Z ~ Normal(0, 1)
  return(Z)
}
#-----------------------------------------------------------------------------------------
```

# Problem A: Stochastic simulation by the probability integral transform and bivariate techniques

## 1.

Let *X* ∼ Exp (*λ*), with the cdf
*F*<sub>*X*</sub>(*x*) = 1 − *e*<sup> − *λ**x*</sup>,  *x* ≥ 0.
Then the random variable *Y* := *F*<sub>*X*</sub>(*X*) has a
Uniform (0, 1) distribution. The probability integral transform becomes

Thus, we sample *Y* from `runif()` and transform it using , to sample
from the exponential distribution. Figure shows one million samples
drawn from the `generate_from_exp()` function defined in the code chunk
below.

``` r
set.seed(123)

generate_from_exp <- function(n, rate = 1) {
  Y <- runif(n)
  X <- -(1 / rate) * log(1 - Y)
  X
}

# sample
n <- 1000000  # One million samples
lambda <- 4.32
x <- generate_from_exp(n, rate = lambda)

# plot
hist(x,
  breaks      = 80,
  probability = TRUE,
  xlim        = c(0, 2)
)
curve(dexp(x, rate = lambda),
  add = TRUE,
  lwd = 2,
  col = "red"
)
```

![Normalized histogram of one million samples drawn from the exponential
distribution, together with the theoretical pdf, with
*λ* = 4.32.](01_Exercise_files/figure-gfm/generate%20from%20exp-1.png)

## 2.

### (a)

More testing

### (b)

## 3.

### (a)

### (b)

### (c)

## 4.

## 5

# Problem B: The gamma distribution

## 1.

### (a)

Let *f*(*x*) be the target distribution we wish to sample from, and let
*g*(*x*) be the proposal distribution. For the rejection sampling
algorithm, we require that for some constant *c* &gt; 0. Let *X* and *U*
be independent samples where *X* ∼ *g*(*x*) and *U* ∼ Uniform (0, 1).
Then the acceptance probability is
$$
\\begin{split}
  \\Pr\\left(U\\leq\\frac{f(X)}{c\\cdot g(X)}\\right) &= \\int\_{-\\infty}^\\infty\\int\_{0}^{f(x)/(c\\ g(x))}f\_{X,U}(x, u)\\,du\\,dx \\\\
  &= \\int\_{-\\infty}^\\infty\\int\_{0}^{f(x)/(c\\ g(x))}g(x)\\cdot 1\\,du\\,dx \\\\
  &= \\int\_{-\\infty}^\\infty\\frac{f(x)}{c\\ g(x)}g(x)\\,dx \\\\
  &= \\frac{1}{c}\\int\_{-\\infty}^\\infty f(x)\\,dx \\\\
  &= \\frac{1}{c}
\\end{split}
$$

We wish to sample from Gamma (*α*, *β* = 1), using the proposal
distribution *g*(*x*) given in `eqref{??????}`. We want to choose *c*
such that the acceptance probability is maximized while is satisfied. We
must check three cases. The trivial case when *x* ≤ 0, we have
*f*(*x*) = *g*(*x*) = 0 so is satisfied for all *c*. When
0 &lt; *x* &lt; 1 we have

$$
\\begin{split}
  f(x) &\\leq c\\,g(x) \\\\
  \\frac{1}{\\Gamma(\\alpha)}x^{\\alpha - 1}e^{-x} &\\leq c\\,\\frac{1}{\\alpha^{-1} + e^{-1}}x^{\\alpha - 1} \\\\
  c&\\geq \\left(\\frac{1}{\\alpha} + \\frac{1}{e}\\right)\\frac{1}{\\Gamma(\\alpha)}e^{-x} \\\\
  c&\\geq \\left(\\frac{1}{\\alpha} + \\frac{1}{e}\\right)\\frac{1}{\\Gamma(\\alpha)}.
\\end{split}
$$

The last case, when *x* ≥ 1, we have
$$
\\begin{split}
  f(x) &\\leq c\\,g(x) \\\\
  \\frac{1}{\\Gamma(\\alpha)}x^{\\alpha - 1}e^{-x} &\\leq c\\,\\frac{1}{\\alpha^{-1} + e^{-1}}e^{-x} \\\\
  c&\\geq \\left(\\frac{1}{\\alpha} + \\frac{1}{e}\\right)\\frac{1}{\\Gamma(\\alpha)}x^{\\alpha - 1} \\\\
  c&\\geq \\left(\\frac{1}{\\alpha} + \\frac{1}{e}\\right)\\frac{1}{\\Gamma(\\alpha)}.
\\end{split}
$$
That is, we choose
*c* := (*α*<sup> − 1</sup> + *e*<sup> − 1</sup>)/*Γ*(*α*), such that the
acceptance probability becomes
$$
\\Pr\\left(U\\leq\\frac{f(X)}{c\\cdot g(X)}\\right) = \\frac{1}{c} = \\frac{\\Gamma(\\alpha)}{\\alpha^{-1} + e^{-1}},\\quad \\alpha\\in(0, 1).
$$

### (b)

``` r
set.seed(137)

sample_from_gamma_rej <- function(n, shape = 0.5) {
  c <- (1 / shape + 1 / exp(1)) / gamma(shape)    # constant that minimizes the envelope
  x <- vector(mode = "numeric", length = n)
  for (i in 1:n) {
    repeat {
      x[i] <- generate_from_gx(1, alpha = shape)  # draw from proposal
      u <- runif(1)                               # draw from U(0, 1)
      f <- dgamma(x[i], shape = shape)            # target value
      g <- theo_gx(x[i], alpha = shape)           # proposal value
      alpha <- (1 / c) * (f / g)
      if (u <= alpha) {
        break
      }
    }
  }
  return(x)
}

# n <- 1000000
# alpha <- 0.9
# x <- sample_from_gamma_rej(n, shape = alpha)
# hist(x,
#   breaks      = 80,
#   probability = TRUE,
#   xlim        = c(0, 6)
# )
# curve(dgamma(x, shape = alpha),
#   add = TRUE,
#   lwd = 2,
#   col = "red"
# )
```

## 2.

### (a)

We will now use the ratio-of-uniforms method to simulate from
Gamma (*α*, *β* = 1). Additionally we have *α* &gt; 1 this time. Let us
define

and

such that
*C*<sub>*f*</sub> ⊂ \[0, 1\] × \[*b*<sub>−</sub>, *b*<sub>+</sub>\].

First we find sup<sub>*x*</sub>*f*<sup>\*</sup>(*x*). This must be when
*x* &gt; 0. We differentiate *f*<sup>\*</sup>(*x*) and setting the
expression equal to zero to find the stationary point.
$$
\\begin{split}
  0 &= \\frac{d}{dx} f^\*(x) \\\\
  &= \\frac{d}{dx} x^{\\alpha - 1}e^{-x} \\\\
  &= e^{-x}x^{\\alpha - 2}\\left((\\alpha - 1) - x\\right) \\\\
  \\Rightarrow\\quad x &= \\alpha - 1,\\quad \\text{where}\\quad\\alpha &gt; 1.
\\end{split}
$$
Since we have only one stationary point, *f*<sup>\*</sup>(*x*) is
continuous, *f*<sup>\*</sup>(*x*) &gt; 0 ∀*x* &gt; 0 and
lim<sub>*x* → 0+</sub>*f*<sup>\*</sup>(*x*) = lim<sub>*x* → ∞</sub>*f*<sup>\*</sup>(*x*) = 0,
then *x* = *α* − 1 must be the global maximum point. That is

We now wish to find *b*<sub>+</sub>.

$$
\\begin{split}
  0 &= \\frac{d}{dx} x^2f^\*(x) \\\\
  &= \\frac{d}{dx} x^{\\alpha + 1}e^{-x} \\\\
  &= e^{-x}x^{\\alpha}\\left((\\alpha + 1) - x\\right) \\\\
  \\Rightarrow\\quad x &= \\alpha + 1,\\quad \\text{where}\\quad\\alpha &gt; 1.
\\end{split}
$$
Using the same reasoning as for *a*, we have that *x* = *α* + 1 is a
global maximum point for *x*<sup>2</sup>*f*<sup>\*</sup>(*x*). Then

Finally, we have that

### (b)

To avoid producing `NaNs`, we will implement the ratio-of-uniform method
on a log scale. We get the following log-transformations.
$$
\\begin{split}
  X\_1\\sim \\operatorname{Uniform}(0, a)&\\Rightarrow \\log X\_1 = \\log a + \\log U\_1,\\quad U\_1\\sim\\operatorname{Uniform}(0,1); \\\\
  X\_2\\sim \\operatorname{Uniform}(b\_-=0, b\_+ = b) &\\Rightarrow \\log X\_2 = \\log b + \\log U\_2,\\quad U\_2\\sim\\operatorname{Uniform}(0,1); \\\\
  y = \\frac{x\_2}{x\_1} &\\Rightarrow y = \\exp\\{(\\log x\_2) - (\\log x\_1)\\}; \\\\
  0\\leq x\_1\\leq \\sqrt{f^\*(y)} &\\Rightarrow \\log x\_1 \\leq \\frac{1}{2}\\log f^\*(y); \\\\
  f^\*(y) = \\begin{cases}y^{\\alpha - 1}e^{-y},&y&gt;0, \\\\ 0, &\\text{otherwise},\\end{cases} &\\Rightarrow \\log f^\*(y) = \\begin{cases}(\\alpha - 1)\\log y - y,&y&gt;0, \\\\ -\\infty, &\\text{otherwise.}\\end{cases}
\\end{split}
$$

``` r
set.seed(434)

lgamma_core <- function(x, alpha = 2) {
  ifelse(
    test = x <= 0,
    yes  = -Inf,
    no   = (alpha - 1)*log(x) - x
  )
}

sample_from_gamma_rou <- function(n, shape = 2, include_trials = FALSE) {
  log_a <- ((shape - 1) / 2) * (log(shape - 1) - 1)
  log_b <- ((shape + 1) / 2) * (log(shape + 1) - 1)
  trials <- 0
  y <- vector(mode = "numeric", length = n)
  for (i in 1:n) {
    repeat {
      log_x1 <- log_a + log(runif(1))
      log_x2 <- log_b + log(runif(1))
      y[i] <- exp(log_x2 - log_x1)
      log_f <- lgamma_core(y[i], alpha = shape)
      if (log_x1 <= 0.5 * log_f) {
        break
      } else {
        trials <- trials + 1
      }
    }
  }
  if (include_trials) {
    return(list(x = y, trials = trials))
  }
  return(y)
}

# generate 1000 samples for each alpha and record number of trials
n <- 1000
m <- 50
alpha <- seq(2, 2000, length.out = m)
trials <- vector(mode = "integer", length = m)

for (i in 1:m) {
  trials[i] <- sample_from_gamma_rou(n, alpha[i], include_trials = TRUE)$trials
}

# plot trials wrt. alpha
ggplot(mapping = aes(x = alpha, y = trials)) +
  geom_point()
```

![Number of trials before accepting *N* = 1000 simulations for various
shape parameters (*α*) using the ratio-of-uniform
method.](01_Exercise_files/figure-gfm/gamma%20ratio-of-uniform%20method-1.png)

Figure strongly suggests that the acceptance probability decreases with
increasing *α*. That is, the ratio of the area of the square
*a* ⋅ (*b*<sub>+</sub> − *b*<sub>−</sub>) = *a**b* and the region
*C*<sub>*f*</sub> is increasing when *α* increases.

## 3.

### (a)

Let *X*<sub>1</sub> and *X*<sub>2</sub> be independent random variables
where *X*<sub>1</sub> ∼ Gamma (*α*<sub>1</sub>, *β* = 1) and
*X*<sub>2</sub> ∼ Gamma (*α*<sub>2</sub>, *β* = 1). Then
*X*<sub>1</sub> + *X*<sub>2</sub> has the following mgf:
$$
\\begin{split}
  M\_{X\_1 + X\_2}(t)  &= M\_{X\_1}(t)\\cdot M\_{X\_2}(t) \\\\
  &= (1 - t)^{-\\alpha\_1}\\cdot(1 - t)^{-\\alpha\_2} \\\\
  &= (1 - t)^{-(\\alpha\_1 + \\alpha\_2)}, \\quad t &lt; 1.
\\end{split}
$$
That is,
*X*<sub>1</sub> + *X*<sub>2</sub> ∼ Gamma (*α*<sub>1</sub> + *α*<sub>2</sub>, *β* = 1).

### (b)

A Exp (1) distribution is a special case of Gamma (*α*, *β*) with
parameters *α* = *β* = 1. Using the result obtained in **(a)**, we then
have that the sum of a random sample
*X*<sub>1</sub>, …, *X*<sub>*k*</sub> drawn from Exp (1) has a
Gamma (*k*, 1) distribution. We can use this fact to improve our
algorithm. Let *k* ∈ ℕ<sub>0</sub> and let *r* ∈ \[0, 1) and assume
*k* + *r* &gt; 0. Then we can decompose any *α* &gt; 0 as
*α* = *k* + *r*.
Let *X*<sub>1</sub>, …, *X*<sub>*k*</sub> ∼ Exp (1) be a random sample
and *W* ∼ Gamma (*r*, 1), where *W* and *X*<sub>*i*</sub>,
*i* = 1, …, *k* are mutual independent. Then

That is, for any *α* &gt; 0, we will only use the rejection sampling
method for *r* = *α*mod 1 to sample *W* ∼ Gamma (*r*, 1), and for the
remaining *k* = *α* − *r* (possibly zero), we sample
*X*<sub>1</sub>, …, *X*<sub>*k*</sub> ∼ Exp (1), and use to sample from
Gamma (*α*, 1).

``` r
sample_from_gamma_improved <- function(n, shape = 1) {
  r <- shape %% 1
  if (r > 0) {
    w <- sample_from_gamma_rej(n, shape = r)
  } else {
    w <- 0
  }
  k <- shape - r
  if (k >= 1) {
    xk <- matrix(generate_from_exp(n * k, rate = 1),
                 nrow = n,
                 ncol = k
    )
    x <- rowSums(xk)
  } else {
    x <- 0
  }
  return(x + w)
}
```

## 4.

Since *β* is an inverse scale parameter, we can simply draw samples from
Gamma (*α*, 1) and divide every sample by *β*. This can be shown by
looking at the mgf of *X*/*β* where *X* ∼ Gamma (*α*, 1).

$$
  M\_{X/\\beta}(t) = M\_X\\left(\\frac{t}{\\beta}\\right) = \\left(1 - \\frac{t}{\\beta}\\right)^{-\\alpha} \\sim \\operatorname{Gamma}(\\alpha,\\beta).
$$

``` r
sample_from_gamma_final <- function(n, shape = 1, rate = 1) {
  (1 / rate) * sample_from_gamma_improved(n, shape = shape)
}
```

## 5.

### (a)

Let *X* and *Y* be independent random variables where
*X* ∼ Gamma (*α*, 1) and *Y* ∼ Gamma (*β*, 1). Let
$$
z = g\_1(x, y) = \\frac{x}{x + y} \\quad \\text{and}\\quad w = g\_2(x, y) = x + y
$$
Then *Z* = *g*<sub>1</sub>(*X*, *Y*) ∈ (0, 1) and
*W* = *g*<sub>2</sub>(*X*, *Y*) &gt; 0. This gives us
$$
\\begin{split}
  &x = g\_1^{-1}(z, w) = zw\\quad\\text{and}\\quad y = g\_2^{-1}(z, w) = w(1 - z), \\\\
  &\|\\det(J)\| = \\left\|\\det\\left(\\begin{matrix}
    \\partial\_zg\_1^{-1}(z, w) & \\partial\_wg\_1^{-1}(z, w) \\\\
    \\partial\_zg\_2^{-1}(z, w) & \\partial\_wg\_2^{-1}(z, w)
  \\end{matrix}\\right)\\right\| = \\left\|\\det\\left(\\begin{matrix}
    w & z \\\\
    -w & 1 - z
  \\end{matrix}\\right)\\right\| = \|w\| = w.
\\end{split}
$$
The marginal distribution *f*<sub>*Z*</sub>(*z*) is then found as
follows.
$$
\\begin{split}
  f\_Z(z) &= \\int\_{0}^\\infty f\_{Z,W}(z, w)\\,dw \\\\
  &=\\int\_0^\\infty f\_{X,Y}\\left(g\_1^{-1}(z, w), g\_2^{-1}(z, w)\\right)\|\\det(J)\|\\,dw \\\\
  &= \\int\_0^\\infty f\_X\\left(zw\\right)\\cdot f\_Y\\left(w(1 - z)\\right)\\cdot w\\,dw \\\\
  &= \\int\_0^\\infty \\frac{1}{\\Gamma(\\alpha)}(zw)^{\\alpha - 1}e^{-zw}\\cdot\\frac{1}{\\Gamma(\\beta)}(w(1 - z))^{\\beta - 1}e^{-w(1 - z)}\\cdot w\\,dw \\\\
  &= \\frac{1}{\\Gamma(\\alpha)\\Gamma(\\beta)}z^{\\alpha - 1}(1 - z)^{\\beta - 1} \\Gamma(\\alpha + \\beta) \\int\_{0}^\\infty \\frac{1}{\\Gamma(\\alpha + \\beta)} w^{(\\alpha + \\beta) - 1}e^{-w}\\,dw \\\\
  &= \\frac{\\Gamma(\\alpha + \\beta)}{\\Gamma(\\alpha)\\Gamma(\\beta)}z^{\\alpha - 1}(1 - z)^{\\beta - 1},\\quad z\\in(0, 1).\\\\
\\end{split}
$$
That is, *f*<sub>*Z*</sub>(*z*) ∼ Beta (*α*, *β*).

### (b)

``` r
sample_from_beta <- function(n, alpha, beta) {
  x <- sample_from_gamma_final(n, shape = alpha)
  y <- sample_from_gamma_final(n, shape = beta)
  return(x / (x + y))
}
```

# Problem C: Monte Carlo integration and variance reduction

## 1.

Let *X* ∼ N (0, 1), and
*θ* = Pr (*X* &gt; 4) ≈ 3.1671242 × 10<sup> − 5</sup>. Let also
*h*(*x*) = *I*(*x* &gt; 4), where *I*( ⋅ ) is the indicator function.
Then
E \[*h*(*X*)\] = ∫<sub> − ∞</sub><sup>∞</sup>*h*(*x*)*f*<sub>*X*</sub>(*x*) *d**x* = ∫<sub> − ∞</sub><sup>∞</sup>*I*(*x* &gt; 4)*f*<sub>*X*</sub>(*x*) *d**x* = Pr (*X* &gt; 4) = *θ*.

Let *X*<sub>1</sub>, …*X*<sub>*n*</sub> ∼ N (0, 1) be a sample. Then the
simple Monte Carlo estimator of *θ* is
$$
  \\hat\\theta\_{\\mathrm{MC}} = \\frac{1}{n}\\sum\_{i=1}^n h(X\_i),
$$
with expectation
$$
  \\operatorname E\\left\[\\hat\\theta\_\\mathrm{MC}\\right\] = \\frac{1}{n}\\sum\_{i=1}^n\\operatorname E\\left\[h(X\_i)\\right\] = \\frac{1}{n}\\sum\_{i=1}^n\\theta = \\theta,
$$

and sampling variance
$$
  \\widehat{\\operatorname{Var}}\\left\[\\hat\\theta\_\\mathrm{MC}\\right\] = \\frac{1}{n^2}\\sum\_{i=1}^n\\widehat{\\operatorname{Var}}\\left\[h(X\_i)\\right\]= \\frac{1}{n}\\widehat{\\operatorname{Var}}\[h(X)\]=\\frac{1}{n(n-1)}\\sum\_{i=1}^n\\left(h(X\_i) - \\hat\\theta\_\\mathrm{MC}\\right)^2.
$$

Then the statistic
$$
  T = \\frac{\\hat\\theta\_\\mathrm{MC} - \\theta}{\\sqrt{\\widehat{\\operatorname{Var}}\\left\[\\hat\\theta\_\\mathrm{MC}\\right\]}}\\sim\\mathrm t\_{n - 1},
$$
and
*t*<sub>*α*/2, *n* − 1</sub> = *F*<sub>*T*</sub><sup> − 1</sup>(1 − *α*/2),
where *F*<sub>*T*</sub><sup> − 1</sup>( ⋅ ) is the quantile function of
the *t*<sub>*n* − 1</sub> distribution.

``` r
set.seed(321)
theta <- pnorm(4, lower.tail = FALSE)
n <- 100000
x <- std_normal(n)
h <- function(x) {
  1 * (x > 4)
}

hh <- h(x)                # I(X > 4) vector of ones and zeros
theta_MC <- mean(hh)      # Monte Carlo estimate of Pr(X > 4)
sample_var_MC <- var(hh)  # Sampling variance

t <- qt(0.05/2, df = n - 1, lower.tail = FALSE)  # quantile with 5% significance level
ci_MC <- theta_MC + t * sqrt(sample_var_MC / n) * c(-1, 1)  # Confidence Interval

# Result
list(
  theta_MC      = theta_MC,
  sample_var_MC = sample_var_MC,
  confint       = ci_MC,
  error         = abs(theta_MC - theta)
)
```

    ## $theta_MC
    ## [1] 6e-05
    ## 
    ## $sample_var_MC
    ## [1] 5.9997e-05
    ## 
    ## $confint
    ## [1] 0.0000119915 0.0001080085
    ## 
    ## $error
    ## [1] 2.832876e-05

## 2.

We will sample from the proposal distribution
$$
  g\_X(x) = \\begin{cases}cx e^{-\\frac{1}{2}x^2}, & x &gt; 4 \\\\ 0, & \\text{otherwise}.\\end{cases}
$$
but first we must find the normalizing constant *c*.

$$
\\begin{split}
  c &= \\left(\\int\_{4}^\\infty x e^{-\\frac{1}{2}x^2}\\,dx\\right)^{-1} 
  = \\left(\\int\_{\\frac{1}{2}4^2}^\\infty e^{-u}\\,du\\right)^{-1} 
  = \\left(e^{-\\frac{1}{2}4^2} - 0\\right)^{-1} 
  = e^{\\frac{1}{2}4^2}, \\\\
  \\Rightarrow g\_X(x) &= \\begin{cases}x e^{-\\frac{1}{2}(x^2 - 4^2)}, & x &gt; 4, \\\\ 0, & \\text{otherwise}.\\end{cases}
\\end{split}
$$

We can easily sample from the proposal distribution using inversion
sampling. The cdf for *x* &gt; 4 is found by integrating.
$$
  G\_X(x) = \\int\_4^x y e^{-\\frac{1}{2}(y^2 - 4^2)}\\,dy = \\int\_0^{\\frac{1}{2}(x^2 - 4^2)}e^{-u}\\,du = 1 - e^{-\\frac{1}{2}(x^2 - 4^2)},\\quad x&gt; 4,
$$
and *G*<sub>*X*</sub>(*x*) = 0 for *x* ≤ 4. Let
*U* = *G*<sub>*X*</sub>(*X*) ∼ Uniform (0, 1). Then we solve for *X*.
$$
\\begin{split}
  U &= 1 - e^{-\\frac{1}{2}(X^2 - 4^2)} \\\\
  -\\frac{1}{2}(X^2 - 4^2) &= \\log(1 - U) \\\\
  X &= \\sqrt{4^2 -2\\log(1 - U)}, \\quad U\\sim \\operatorname{Uniform}(0, 1).
\\end{split}
$$

Let *X*<sub>1</sub>, …, *X*<sub>*n*</sub> be a sample drawn from the
proposal distribution *g*<sub>*X*</sub>(*x*). Then the importance
sampling estimator of *θ* is given by
$$
  \\hat\\theta\_\\mathrm{IS} = \\frac{1}{n}\\sum\_{i=1}^n h(X\_i)w(X\_i),
$$
where *w*(*x*) = *f*<sub>*X*</sub>(*x*)/*g*<sub>*X*</sub>(*x*), with
expectation

$$
\\begin{split}
\\operatorname E\\left\[\\hat\\theta\_\\mathrm{IS}\\right\] &= \\frac{1}{n}\\sum\_{i=1}^n\\int\_{0}^\\infty h(x\_i)w(x\_i)g\_X(x\_i)\\,dx\_i \\\\
&= \\frac{1}{n}\\sum\_{i=1}^n\\int\_{0}^\\infty h(x\_i)f\_X(x\_i)\\,dx\_i \\\\
&= \\frac{1}{n}\\sum\_{i=1}^n\\operatorname E\\left\[h(X\_i)\\mid X\_i\\sim \\operatorname{N}(0, 1)\\right\] \\\\
&= \\frac{1}{n}\\sum\_{i=1}^n\\theta \\\\
&= \\theta,
\\end{split}
$$

and sampling variance

$$
  \\widehat{\\operatorname{Var}}\\left\[\\hat\\theta\_\\mathrm{IS}\\right\] = \\frac{1}{n^2}\\sum\_{i=1}^n\\widehat{\\operatorname{Var}}\\left\[h(X\_i)w(X\_i)\\right\]= \\frac{1}{n}\\widehat{\\operatorname{Var}}\[h(X)w(X)\]=\\frac{1}{n(n-1)}\\sum\_{i=1}^n\\left(h(X\_i)w(X\_i) - \\hat\\theta\_\\mathrm{IS}\\right)^2.
$$

``` r
set.seed(321)

sample_from_proposal <- function(n) {
  u <- runif(n)
  sqrt(4^2 - 2 * log(1 - u))
}

n <- 100000
x <- sample_from_proposal(n)

w <- function(x) {
  f <- dnorm(x)                              # target density
  g <- ifelse(                               # proposal density
         test = x > 4,
         yes  = x * exp(-0.5 * (x^2 - 16)),
         no   = 0
       )
  return(f / g)
}

hw <- h(x) * w(x)

theta_IS <- mean(hw)      # Importance sampling estimate of Pr(X > 4)
sample_var_IS <- var(hw)  # Sampling variance

t <- qt(0.05/2, df = n - 1, lower.tail = FALSE)  # quantile with 5% significance level
ci_IS <- theta_IS + t * sqrt(sample_var_IS / n) * c(-1, 1)  # Confidence Interval

# Result
list(
  theta_IS      = theta_IS,
  sample_var_IS = sample_var_IS,
  confint       = ci_IS,
  error         = abs(theta_IS - theta)
)
```

    ## $theta_IS
    ## [1] 3.167611e-05
    ## 
    ## $sample_var_IS
    ## [1] 2.410122e-12
    ## 
    ## $confint
    ## [1] 3.166649e-05 3.168573e-05
    ## 
    ## $error
    ## [1] 4.866683e-09

The number of samples *m* needed for the simple Monte Carlo estimator to
achieve the same precision as the importance sampling approach, we would
need

$$
  m = n\\frac{\\widehat{\\operatorname{Var}}\[h(X)\]}{\\widehat{\\operatorname{Var}}\[h(X)w(X)\]} = 10^{5}\\frac{5.9997\\times 10^{-5}}{2.4101218\\times 10^{-12}} = 2.4893763\\times 10^{12},
$$
samples. That is, we need about 10 million times more samples.

## 3.

### (a)

We modify `sample_from_proposal()` to return a pair of samples, where
one takes *U* ∼ Uniform (0, 1) as argument and the other 1 − *U* as
argument.

``` r
sample_from_proposal_mod <- function(n) {
  u <- runif(n)
  list(
    x_1 = sqrt(4^2 - 2 * log(1 - u)),
    x_2 = sqrt(4^2 - 2 * log(u))
  )
}
```

### (b)

``` r
set.seed(53)
n <- 50000
sample_pair <- sample_from_proposal_mod(n)

hw1 <- h(sample_pair$x_1) * w(sample_pair$x_1)
hw2 <- h(sample_pair$x_2) * w(sample_pair$x_2)
hw_AS <- 0.5 * (hw1 + hw2)  # Antithetic sample

theta_AS <- mean(hw_AS)      # Antithetic sampling estimate of Pr(X > 4)
sample_var_AS <- var(hw_AS)  # Sampling variance

t <- qt(0.05/2, df = n - 1, lower.tail = FALSE)  # quantile with 5% significance level
ci_AS <- theta_AS + t * sqrt(sample_var_AS / n) * c(-1, 1)  # Confidence Interval

# Result
list(
  theta_AS      = theta_AS,
  sample_var_AS = sample_var_AS,
  confint       = ci_AS,
  error         = abs(theta_AS - theta)
)
```

    ## $theta_AS
    ## [1] 3.167307e-05
    ## 
    ## $sample_var_AS
    ## [1] 2.849882e-13
    ## 
    ## $confint
    ## [1] 3.166839e-05 3.167774e-05
    ## 
    ## $error
    ## [1] 1.823745e-09

# Problem D: Rejection sampling and importance sampling

## 1.

## 2.

## 3.

## 4.
