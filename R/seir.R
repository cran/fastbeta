seir <-
function (length.out = 1L,
          beta, nu, mu, sigma = gamma, gamma = 1, delta = 0,
          init, m = length(init) - n - 2L, n = 1L,
          stochastic = TRUE, prob = 1, delay = 1,
          useCompiled = TRUE, ...)
{
	stopifnot(requireNamespace(if (stochastic) "adaptivetau" else "deSolve"),
	          is.integer(length.out) && length(length.out) == 1L && length.out >= 1L,
	          is.function(beta) && !is.null(formals(beta)),
	          is.function(nu  ) && !is.null(formals(nu  )),
	          is.function(mu  ) && !is.null(formals(mu  )),
	          is.double(sigma) && length(sigma) == 1L && sigma >= 0,
	          is.double(gamma) && length(gamma) == 1L && gamma >= 0,
	          is.double(delta) && length(delta) == 1L && sigma >= 0,
	          is.integer(m) && length(m) == 1L && m >= 0L && m < 4096L,
	          is.integer(n) && length(n) == 1L && n >= 1L && n < 4096L,
	          is.double(init),
	          length(init) == m + n + 2L,
	          all(is.finite(init)),
	          min(init) >= 0,
	          is.double(prob),
	          any(length(prob) == c(1L, length.out - 1L)),
	          min(prob) >= 0,
	          max(prob) <= 1,
	          is.double(delay),
	          length(delay) >= 1L,
	          min(delay) >= 0,
	          sum(delay) > 0)

	p <- m + n + 2L

	init. <- c(init, 0, 0)
	names(init.) <- nms <-
		c("S",
		  sprintf("E.%03x", seq_len(m)),
		  sprintf("I.%03x", seq_len(n)),
		  "R",
		  "Z", # incidence
		  "B") # births

	if (!useCompiled) {
		sigma <- sigma * m
		gamma <- gamma * n
		delta <- delta * 1
		i.S <- 1L
		i.E <- seq.int(    2L, length.out = m)
		i.I <- seq.int(m + 2L, length.out = n)
		i.R <- p
	}

	if (stochastic) {

		if (any(init. != (tmp <- trunc(init.)))) {
			warning(gettextf("truncating fractional part of '%s'",
			                 "init"),
			        domain = NA)
			init. <- tmp
		}

		tl.params <-
			(function (maxtau = 1, ...) list(maxtau = maxtau, ...))(...)

		tran <-
			c(list(`names<-`(c(-1, 1, 1), c("S", nms[2L], "Z"))),
			  list(`names<-`(c(1, 1), c("S", "B"))),
			  lapply(1L:p,
			         function(i) `names<-`(-1, nms[i])),
			  lapply(2L:(p - 1L),
			         function(i) `names<-`(c(-1, 1), nms[c(i, i + 1L)])),
			  list(`names<-`(c(1, -1), c("S", "R"))))
		## infection, birth, natural mortality, removal, loss of immunity

		if (useCompiled) {
			.Call(R_adseir_initialize, beta, nu, mu, sigma, gamma, delta, m, n)
			on.exit(.Call(R_adseir_finalize))
			ff <- function (x, theta, t) .Call(R_adseir_dot, t, x)
			Df <- function (x, theta, t) .Call(R_adseir_jac, t, x)
		}
		else {
			## D[i, j] = d(rate of transition j)/d(state i)
			##
			##     [S   ]  x  0  mu    0   0   0      0      0      0
			##     [E[i]]  0  0   0   mu   0   0  sigma      0      0
			##     [I[j]]  y  0   0    0  mu   0      0  gamma      0
			##     [R   ]  0  0   0    0   0  mu      0      0  delta
			##     [Z   ]  0  0   0    0   0   0      0      0      0
			##     [B   ]  0  0   0    0   0   0      0      0      0
			##
			## where x = beta * sum(I), y = beta * S
			D <- matrix(0, p + 2L, p + p + 1L)
			k.0 <- seq.int(from = m + 2L, by = 1L,
			               length.out = n)
			k.1 <- seq.int(from = (p + 2L) * 2L + 1L, by = p + 3L,
			               length.out = p)
			k.2 <- seq.int(from = (p + 2L) * (p + 2L) + 2L, by = p + 3L,
			               length.out = p - 1L)
			D[k.2] <- rep.int(c(sigma, gamma, delta), c(m, n, 1L))

			ff <-
			function (x, theta, t)
			{
				x.S <- x[i.S]
				x.E <- x[i.E]
				x.I <- x[i.I]
				x.R <- x[i.R]
				beta <- beta(t)
				nu   <- nu  (t)
				mu   <- mu  (t)
				ok <- sum(x.E) + (s <- sum(x.I)) > 1
				c(        beta * s * x.S               ,
				                nu                     ,
				                mu * x.S               ,
				  if (ok)       mu * x.E else double(m),
				  if (ok)       mu * x.I else double(n),
				                mu * x.R               ,
				             sigma * x.E               ,
				  if (ok)    gamma * x.I else double(n),
				             delta * x.R               )
			}
			Df <-
			function (x, theta, t)
			{
				x.S <- x[i.S]
				x.I <- x[i.I]
				beta <- beta(t)
				mu   <- mu  (t)
				D[ 1L] <<- beta * sum(x.I)
				D[k.0] <<- beta * x.S
				D[k.1] <<- mu
				D
			}
		}

		X. <- adaptivetau::ssa.adaptivetau(
			init.values  = init.,
			transitions  = tran,
			rateFunc     = ff,
			params       = NULL,
			tf           = length.out - 1L,
			jacobianFunc = Df,
			tl.params    = tl.params)

		D. <- dim(X.)
		i <- D.[1L] - match(seq.int(0L, length.out = length.out),
		                    as.integer(ceiling(X.[D.[1L]:1L, 1L]))) + 1L
		if (anyNA(i)) {
			## tl.params[["maxtau"]] constrains leaps but not steps => LOCF
			i[is.na(i)] <- 0L
			k <- c(which(i[-1L] != i[-length(i)]), length(i)) # run stops
			ik <- i[k]
			w <- which(ik == 0L)
			ik[w] <- ik[w - 1L]
			i <- rep.int(ik, k - c(0L, k[-length(k)]))
		}
		X <- X.[i, -1L, drop = FALSE] # discarding time

	}
	else {

		## E[i], I[j], R: logarithmic scale
		stopifnot(min(init[-1L]) > 0)
		init.[2L:p] <- log(init.[2L:p])

		if (useCompiled) {
			.Call(R_deseir_initialize, beta, nu, mu, sigma, gamma, delta, m, n)
			on.exit(.Call(R_deseir_finalize))
			gg <- "R_deseir_dot"
			Dg <- "R_deseir_jac"
		}
		else {
			## D[i, j] = d(rate of change in state i)/d(state j)
			##
			## nonzero pattern in m = 4, n = 4 case :
			##
			##       S E E E E I I I I R Z B
			##     S | . . . . | | | | | . .
			##     E | | . . . | | | | . . .
			##     E . | | . . . . . . . . .
			##     E . . | | . . . . . . . .
			##     E . . . | | . . . . . . .
			##     I . . . . | | . . . . . .
			##     I . . . . . | | . . . . .
			##     I . . . . . . | | . . . .
			##     I . . . . . . . | | . . .
			##     R . . . . . . . . | | . .
			##     Z | . . . . | | | | . . .
			##     B . . . . . . . . . . . .
			##
			## nonzero pattern in m = 0, n = 4 case :
			##
			##       S I I I I R Z B
			##     S | | | | | | . .
			##     I | . | | | . . .
			##     I . | | . . . . .
			##     I . . | | . . . .
			##     I . . . | | . . .
			##     R . . . . | | . .
			##     Z | | | | | . . .
			##     B . . . . . . . .
			##
			## where log(.) is suppressed only for pretty printing
			D <- matrix(0, p + 2L, p + 2L)
			k.0 <- seq.int(from = (p + 2L) * 2L + 3L, by = p + 3L,
			               length.out = p - 2L)
			k.1 <- k.0 - p - 2L

			i.1 <- 2L:(p - 1L)
			i.0 <- 3L:p
			a.1 <- rep.int(c(sigma, gamma), c(m, n)); a.11 <- a.1[1L]
			a.0 <- c(a.1[-1L], delta)

			gg <-
			function (t, x, theta)
			{
				x.S <- x[i.S]
				x.E <- x[i.E]
				x.I <- x[i.I]
				x.R <- x[i.R]
				beta <- beta(t)
				nu   <- nu  (t)
				mu   <- mu  (t)
				s.1 <- sum(exp(x.I        ))
				s.2 <- sum(exp(x.I - x[2L]))
				list(c(nu + delta * exp(x.R) - (beta * s.1 + mu) * x.S,
				       beta * s.2 * x.S - (a.11 + mu),
				       a.1 * exp(x[i.1] - x[i.0]) - (a.0 + mu),
				       beta * s.1 * x.S,
				       nu))
			}
			Dg <-
			function (t, x, theta)
			{
				x.S <- x[i.S]
				x.I <- x[i.I]
				x.R <- x[i.R]
				beta <- beta(t)
				mu   <- mu  (t)
				s.1 <- sum(u.1 <- exp(x.I        ))
				s.2 <- sum(u.2 <- exp(x.I - x[2L]))
				D[i.S, i.S] <<- -(D[p + 1L, i.S] <<- beta * s.1) - mu
				D[i.S, i.I] <<- -(D[p + 1L, i.I] <<- beta * u.1 * x.S)
				D[i.S, i.R] <<- delta * exp(x.R)
				D[ 2L, i.S] <<- beta * s.2
				D[ 2L, i.I] <<- beta * u.2 * x.S
				D[ 2L,  2L] <<- if (m) -beta * s.2 else 0
				D[k.0] <<- -(D[k.1] <<- a.1 * exp(x[i.1] - x[i.0]))
				D
			}
		}

		X. <- deSolve::lsoda(
			y        = init.,
			times    = seq.int(0, length.out = length.out),
			func     = gg,
			parms    = NULL,
			jacfunc  = Dg,
			jactype  = "fullusr",
			hmax     = 1,
			ynames   = FALSE,
			dllname  = if (useCompiled) "fastbeta",
			initfunc = NULL,
			initpar  = NULL,
			...)

		X. <- X.[, -1L, drop = FALSE]
		X.[, 2L:p] <- exp(X.[, 2L:p])
		D. <- dim(X.)
		X <-
			if (D.[1L] < length.out)
				## not seen, but help("lsoda") says that it is possible
				rbind(X., matrix(NaN, length.out - D.[1L], D.[2L]),
				      deparse.level = 0L)
			else X.

	}

	if (length.out > 1L) {
	head <- 1L:(length.out - 1L)
	tail <- 2L:length.out
	X[tail, p + 1L:2L] <- X[tail, p + 1L:2L] - X[head, p + 1L:2L]
	}
	X[  1L, p + 1L:2L] <- NA_real_

	m.p <- missing(prob)
	m.d <- missing(delay)
	if (doObs <- !(m.p && m.d)) {
		X <- cbind(X, NA_real_, deparse.level = 0L)
		if (length.out > 1L) {
		Z <- X[tail, p + 1L]
		if (stochastic) {
			Z <- as.integer(Z)
			if (!m.p)
				Z <- rbinom(length.out - 1L, Z, prob)
			if (!m.d)
				## FIXME? 'rmultinom' is more efficient, but not vectorized ...
				Z <- tabulate(rep.int(seq_len(length.out - 1L), Z) +
				              sample(seq.int(0L, length.out = length(delay)),
				                     size = sum(Z),
				                     replace = TRUE,
				                     prob = delay),
				              length.out - 1L)
		}
		else {
			if (!m.p)
				Z <- Z * prob
			if (!m.d) {
				d <- length(delay) - 1L
				Z <- filter(c(double(d), Z), delay / sum(delay),
				            sides = 1)[(d + 1L):(d + length.out - 1L)]
			}
		}
		X[tail, p + 3L] <- Z
		}
	}

	nms <- rep.int(c("S", "E", "I", "R", "Z", "B", "Z.obs"),
	               c(1L, m, n, 1L, 1L, 1L, if (doObs) 1L else 0L))
	ts(X, start = 0, names = nms)
}

seir.R0 <-
function (beta, nu, mu, sigma, gamma, delta, m = 0L, n = 1L)
{
	sigma <- sigma * m
	gamma <- gamma * n
	delta <- delta * 1
	(nu / mu) * (sigma / (sigma + mu))^m * (beta / (gamma + mu)) *
		sum(cumprod(rep.int(c(1, gamma / (gamma + mu)), c(1L, n - 1L))))
}

seir.ee <-
function (beta, nu, mu, sigma, gamma, delta, m = 0L, n = 1L)
{
	sigma <- sigma * m
	gamma <- gamma * n
	delta <- delta * 1
	y <- cumprod(rep.int((c(delta, gamma) + mu) / gamma, c(1L, n - 1L)))
	if (m == 0L) {
		S. <- (tmp <- (gamma + mu) * y[n]) / beta / sum(y)
		R. <- (nu - mu * S.) / (tmp - delta)
		ans <- c(S., y[n:1L] * R., R.)
	}
	else {
		x <- cumprod(rep.int((c(gamma, sigma) + mu) / sigma, c(1L, m - 1L))) *
			y[n]
		S. <- (tmp <- (sigma + mu) * x[m]) / beta / sum(y)
		R. <- (nu - mu * S.) / (tmp - delta)
		ans <- c(S., x[m:1L] * R., y[n:1L] * R., R.)
	}
	names(ans) <- rep.int(c("S", "E", "I", "R"), c(1L, m, n, 1L))
	ans
}

sir <-
function (length.out = 1L,
          beta, nu, mu, gamma = 1, delta = 0,
          init, n = 1L,
          stochastic = TRUE, prob = 1, delay = 1,
          useCompiled = TRUE, ...)
{
	if (any(...names() == "m"))
		stop(gettextf("call '%s', not '%s', when setting '%s'",
		              "seir", "sir", "m"),
		     domain = NA)
	call <- match.call(expand.dots = TRUE)
	call[[1L]] <- quote(seir)
	call[["m"]] <- 0L
	eval(call, parent.frame())
}
