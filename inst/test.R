rungekutta <-
function (f, t, y, theta, h, length.out, a, b, c)
{
	stopifnot(length(h) == 1L, h > 0,
	          length(b) > 0L, length(c) == length(b),
	          length(a) == (length(b) * (length(b) - 1L)) %/% 2L)
	F <- vector("list", length(b))
	ans <- vector("list", length.out)
	ans[[1L]] <- y
	for (k in seq_len(length(ans) - 1L)) {
		t. <- t + (k - 1L) * h
		y. <- y
		pos <- 0L
		for (j in seq_along(F)) {
			tmp <- 0
			for (i in seq_len(j - 1L))
				tmp <- tmp + a[pos <- pos + 1L] * F[[i]]
			F[[j]] <- f(t. + c[j] * h, y. + tmp * h, theta)
		}
		tmp <- 0
		for (j in seq_along(F))
			tmp <- tmp + b[j] * F[[j]]
		ans[[k + 1L]] <- y <- y + tmp * h
	}
	ans
}

dormandprince <- # ode45 in MATLAB
function (f, t, y, theta, h, length.out)
{
	a <- .fmpq(num = c(1L,
	                   3L, 9L,
	                   44L, -56L, 32L,
	                   19372L, -25360L, 64448L, -212L,
	                   9017L, -355L, 46732L, 49L, -5103L,
	                   35L, 0L, 500L, 125L, -2187L, 11L),
	           den = c(5L,
	                   40L, 40L,
	                   45L, 15L, 9L,
	                   6561L, 2187L, 6561L, 729L,
	                   3168L, 33L, 5247L, 176L, 18656L,
	                   384L, 1L, 1113L, 192L, 6784L, 84L))
	b <- .fmpq(num = c(5179L, 0L, 7571L, 393L, -92097L, 187L, 1L),
	           den = c(57600L, 1L, 16695L, 640L, 339200L, 2100L, 40L))
	c <- .fmpq(num = c(0L, 1L, 3L, 4L, 8L, 1L, 1L),
	           den = c(1L, 5L, 10L, 5L, 9L, 1L, 1L))
	rungekutta(f = f, t = t, y = y, theta = theta,
	           h = h, length.out = length.out,
	           a = a, b = b, c = c)
}

library(flint)
## flintPrec(prec = .Machine[["double.digits"]])
## flintPrec(prec = .Machine[["double.max.exp"]])
flintPrec(prec = .Machine[["double.max.exp"]] * 0x1p+4L)


f <- function (t, y, theta) c(.arb(x = -1), exp(y[1L]))
t <- .arb(x = 0)
y <- .arb(x = c(0, 0))
theta <- NULL
h <- .arb(x = 0x1p-12)
length.out <- 1L + as.integer(0x1p+6/h)

L <- dormandprince(f = f, t = t, y = y, theta = theta,
                   h = h, length.out = length.out)

last <- L[[length(L)]]
last[2L] - (1 - exp(.arb(x = -0x1p+7)))


## d/dt a(t) = -a(t)    =>   d/dt log(a(t)) = -1
## d/dt b(t) =  a(t)
x <- deSolve::lsoda(y = c(0, 0),
                    times = seq.int(from = 0, to = 100, by = 0.01),
                    func = function (t, x, theta) list(c(-1, exp(x[1L]))),
                    parms = NULL,
                    ynames = FALSE)

## log(a(t)) = log(exp(-t)) = -t
x[, 2L] - (-x[, 1L])
## b(t) = integrate(exp(-s), 0, t) = 1 - exp(-t)
x[, 3L] - (1 - exp(-x[, 1L]))

plot(1 - x[, 3L], log = "y")
x[, 3L]











R0 <- 16
ell <- 0.5

init.infected <- 1e-6
init.infectious <- init.infected * 1e-3

init <- c(1 - init.infected, # S
          log(init.infected - init.infectious), # log(E)
          log(init.infectious), # log(I)
          log(init.infected - init.infectious), # log(Y_E)
          log(init.infectious)) # log(Y_I)

a <- 1/ell

F <-
function (t, x, theta)
{
	    S <- x[1L]
	log.E <- x[2L]
	log.I <- x[3L]
	log.YE <- x[4L]
	log.YI <- x[5L]
    list(c(-R0 * exp(log.I) * S, # d/dt S
	       R0 * exp(log.I - log.E) * S - 1/ell, # d/dt log(E)
	       1/ell * exp(log.E - log.I) - 1, # d/dt log(I)
           R0 * exp(log.I - log.YE) * S - a, # d/dt log(Y_E)
           a * exp(log.YE - log.YI) - exp(log.I - log.YI))) # d/dt log(Y_I)
}

## Y_E' = F * S - a * Y_E
## Y_I' = a * Y_E - F * (1/R0)

## Y' = (Y_E + Y_I)' = F * (S - 1/R0)


from <- 0; to <- 100; by <- 0.01
x. <- deSolve::lsoda(y = init,
                     times = seq.int(from = from, to = to, by = by),
                     func = F,
                     parms = NULL,
                     ynames = FALSE)
plot(exp(x.[, 5L]) + exp(x.[, 6L]), log = "y")
x <- ts(cbind(x.[, 2L], exp(x.[, 3L:4L]), ),
        start = from, end = to, deltat = by,
        names = c("S", "E", "I", "Y"))
plot(x, log = "y")
tail(x)




