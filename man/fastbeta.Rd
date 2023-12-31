\name{fastbeta}
\alias{fastbeta}
\title{Estimate a Time-Varying Infectious Disease Transmission Rate}
\description{
Generates a discrete time estimate of a transmission rate
\eqn{\beta(t)} from an equally spaced incidence time series
and other data.
}
\details{
\code{fastbeta} works by discretizing the system of ordinary
differential equations
\deqn{
  \begin{aligned}
  \frac{\text{d} S}{\text{d} t} &= \nu(t) - \beta(t) S I + \delta R - \mu(t) S \\
  \frac{\text{d} I}{\text{d} t} &= \beta(t) S I - \gamma I - \mu(t) I \\
  \frac{\text{d} R}{\text{d} t} &= \gamma I - \delta R - \mu(t) R
  \end{aligned}
}{
  dS/dt = nu(t) - beta(t) S I + delta R - mu(t) S
  dI/dt = beta(t) S I - gamma I - mu(t) I
  dR/dt = gamma I - delta R - mu(t) R
}
where \eqn{t} is understood to be a unitless measure of time
relative to the duration of an observation interval, then computing
the iteration
\deqn{
  \begin{aligned}
  I_{t + 1} &= \frac{(1 - \frac{1}{2} (\gamma + \mu_{t})) I_{t} + Z_{t + 1}}{
    1 + \frac{1}{2} (\gamma + \mu_{t + 1})} \\
  R_{t + 1} &= \frac{(1 - \frac{1}{2} (\delta + \mu_{t})) R_{t} + \frac{1}{2} \gamma (I_{t} + I_{t + 1})}{
    1 + \frac{1}{2} (\delta + \mu_{t + 1})} \\
  S_{t + 1} &= \frac{(1 - \frac{1}{2} \mu_{t}) S_{t} - Z_{t + 1} + B_{t + 1} + \frac{1}{2} \delta (R_{t} + R_{t + 1})}{
    1 + \frac{1}{2} \mu_{t + 1}} \\
  \beta_{t} &= \frac{Z_{t} + Z_{t + 1}}{2 S_{t} I_{t}}
  \end{aligned}
}{
  I[t + 1] = ((1 - 0.5 * (gamma + mu[t])) * I[t] + Z[t + 1]) /
      (1 + 0.5 * (gamma + mu[t + 1]))
  R[t + 1] = ((1 - 0.5 * (delta + mu[t])) * R[t] + 0.5 * gamma * (I[t] + I[t + 1]))/
      (1 + 0.5 * (delta + mu[t + 1]))
  S[t + 1] = ((1 - 0.5 *          mu[t] ) * S[t] - Z[t + 1] + B[t + 1] + 0.5 * delta * (R[t] + R[t + 1])) /
      (1 + 0.5 *          mu[t + 1] )
  beta[t] = (Z[t] + Z[t + 1]) / (2 * S[t] * I[t])
}
where
\deqn{
  \begin{aligned}
  X_{t} &\sim X(t) \quad [X = \beta, \mu, S, I, R] \\
  Z_{t} &\sim \int_{t - 1}^{t} \beta(s) S(s) I(s) \, \text{d} s \\
  B_{t} &\sim \int_{t - 1}^{t} \nu(s) \, \text{d} s
  \end{aligned}
}{
  X[t] ~ X(t)    [X = S, I, R, beta, mu]
  Z[t] ~ integrate(beta(s) S(s) I(s), t - 1, t)
  B[t] ~ integrate(            nu(s), t - 1, t)
}
and it is understood that indexing starts at \eqn{t = 0}.
\eqn{Z_{t}}{Z[t]}, \eqn{B_{t}}{B[t]}, and \eqn{\mu_{t}}{mu[t]} denote
incidence and births between times \eqn{t - 1} and \eqn{t} and the
per capita natural mortality rate at time \eqn{t}; they are supplied
together as argument \code{series}.
}
\usage{
fastbeta(series, constants, \dots)
}
\arguments{
\item{series}{a \dQuote{multiple time series} object, inheriting from
  class \code{\link[=ts]{mts}}, with three columns storing
  (\dQuote{parallel}, equally spaced) time series of incidence, births,
  and the per capita natural mortality rate, in that order.}
\item{constants}{a numeric vector of the form
  \code{c(S0, I0, R0, gamma, delta)}, specifying an initial state
  and rates of removal and loss of immunity, in that order.}
\item{\dots}{optional arguments passed to \code{\link{deconvolve}},
  if the first column of \code{series} represents \emph{reported}
  incidence or mortality rather than incidence.}
}
\value{
A \dQuote{multiple time series} object, inheriting from class
\code{\link[=ts]{mts}}, with four columns (named \code{S}, \code{I},
\code{R}, and \code{beta}) storing the result of the iteration
described in \sQuote{Details}.
It is completely parallel to argument \code{series}, having the same
\code{\link{tsp}} attribute.
}
%\seealso{}
\references{
Jagan, M., deJonge, M. S., Krylova, O., & Earn, D. J. D. (2020).
Fast estimation of time-varying infectious disease transmission rates.
\emph{PLOS Computational Biology},
\emph{16}(9), Article e1008124, 1-39.
\doi{10.1371/journal.pcbi.1008124}
}
\examples{
\dontshow{
## for R_DEFAULT_PACKAGES=NULL
library(   stats, pos = "package:base", verbose = FALSE)
library(graphics, pos = "package:base", verbose = FALSE)
library(   utils, pos = "package:base", verbose = FALSE)
}
data(sir.e02, package = "fastbeta")
a <- attributes(sir.e02)
str(sir.e02)
plot(sir.e02)

## We suppose that we have perfect knowledge of incidence,
## births, and the data-generating parameters
series <- cbind(sir.e02[, c("Z", "B")], mu = a[["mu"]](0))
colnames(series) <- c("Z", "B", "mu") # FIXME: stats:::cbind.ts mangles dimnames
constants <- c(S0 = sir.e02[[1L, "S"]],
               I0 = sir.e02[[1L, "I"]],
               R0 = sir.e02[[1L, "R"]],
               gamma = a[["gamma"]],
               delta = a[["delta"]])

X <- fastbeta(series, constants)
str(X)
plot(X)

plot(X[, "beta"], ylab = "Transmission rate")
lines(a[["beta"]](time(X)), col = "red") # the "truth"
}
