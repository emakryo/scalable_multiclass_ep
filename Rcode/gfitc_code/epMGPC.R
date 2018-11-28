#########################################################################################################
#
# Implements a multi-class Gaussian Process classification using parallel EP updates for the likelihood.
# Uses the FITC approximation to approximate the covariance matrix of the prior. 
#

t0 <- NULL

##
# Function which computes the cholesky decomposition of the inverse
# of a particular matrix.
#
# @param	M	m x m positive definite matrix.
#
# @return	L	m x m upper triangular matrix such that
#			M^-1 = L %*% t(L)
#

cholInverse <- function(M) { rot180(forwardsolve(t(chol(rot180(M))), diag(nrow(M)))) }

##
# Function which rotates a matrix 180 degreees.
#

rot180 <- function(M) { matrix(rev(as.double(M)), nrow(M), ncol(M)) }

##
# This function computes the covariance matrix for the GP
#

kernel <- function(X, l, sigma0, sigma) {
	X <- X / matrix(sqrt(l), nrow(X), ncol(X), byrow = TRUE)
	distance <- as.matrix(dist(X))^2
	sigma * exp(-0.5 * distance) + diag(sigma0, nrow(X)) + diag(rep(1e-10, nrow(X)))
}

##
# Function which computes the kernel matrix between the observed data and the test data
#

kernel_nm <- function(X, Xnew, l, sigma) {
	X <- X / matrix(sqrt(l), nrow(X), ncol(X), byrow = TRUE)	
	Xnew <- Xnew / matrix(sqrt(l), nrow(Xnew), ncol(Xnew), byrow = TRUE)
	n <- nrow(X)
	m <- nrow(Xnew)
	Q <- matrix(apply(X^2, 1, sum), n, m)
	Qbar <- matrix(apply(Xnew^2, 1, sum), n, m, byrow = T)
	distance <- Qbar + Q - 2 * X %*% t(Xnew)
	sigma * exp(-0.5 * distance)
}

##
# Function which computes the diagonal of the kernel matrix for the data
# points.
#
# @param	X 		n x d matrix with the n data points.
# @param	sigma		scalar with the amplitude of the GP.
# @param	sigma0		scalar with the noise level in the GP.
#
# @return	diagKnn		n-dimensional vector with the diagonal of the
#				kernel matrix for the data points.
#

computeDiagKernel <- function(X, sigma, sigma0) { rep(sigma, nrow(X)) + 1e-10 + sigma0 }

##
# Function that initializes the struture with the problem information.
#
# @param	X	n x d matrix with the data points.
# @param	Xbar	m x d matrix with the pseudo inputs.
# @param	sigma	scalar with the log-amplitude of the GP.
# @param	sigma0	scalar with the log-noise level in the GP.
# @param	l	d-dimensional vector with the log-lengthscales.
# 
# @return	gFITCinfo	List with the problem information 
#

initialize_kernel_FITC <- function(X, Xbar, sigma, sigma0, l) {

	# We initialize the structure with the data and the kernel
	# hyper-parameters

	gFITCinfo <- list()
	gFITCinfo$X <- X
	gFITCinfo$Xbar <- Xbar
	gFITCinfo$m <- nrow(Xbar)
	gFITCinfo$d <- ncol(Xbar)
	gFITCinfo$n <- nrow(X)
	gFITCinfo$sigma <- sigma
	gFITCinfo$sigma0 <- sigma0
	gFITCinfo$l <- l

	# We compute the kernel matrices
	#browser()
	gFITCinfo$Kmm <- kernel(Xbar, gFITCinfo$l, gFITCinfo$sigma0, gFITCinfo$sigma)
	gFITCinfo$Knm <- kernel_nm(X, Xbar, gFITCinfo$l, gFITCinfo$sigma)
	gFITCinfo$P <- gFITCinfo$Knm
	gFITCinfo$R <- cholInverse(gFITCinfo$Kmm)
	gFITCinfo$PRt <- gFITCinfo$P %*% t(gFITCinfo$R)

	# We compute the diagonal matrices

	gFITCinfo$diagKnn <- computeDiagKernel(X, gFITCinfo$sigma, gFITCinfo$sigma0)
	gFITCinfo$D <- gFITCinfo$diagKnn - (gFITCinfo$PRt)^2 %*% rep(1, gFITCinfo$m)

	gFITCinfo
}

##
# Function that computes the marginal means and variances of the product
# of the FITC prior and a multivariate Gaussian density with diagonal
# correlation matrix.
#
# @param	gFITCinfo	List with the problem information
#				(see initializegFITCinfo).
# @param	f1Hat		list with n-dimensional vector with the inverse
#				variances times the mean of the Gaussian and the invese variances
#				approximation to the likelihood factor.
#
# @return	ret		A list with the marginal means and variances.
#

reconstructPosterior <- function(gFITCinfo, f1Hat) {

	Dnew <- gFITCinfo$D / (1 + gFITCinfo$D * f1Hat$eta2)
	Pnew <- matrix(1 / (1 + gFITCinfo$D * f1Hat$eta2), gFITCinfo$n,
		gFITCinfo$m) * gFITCinfo$P

	Rnew <- backsolve(rot180(t(chol(rot180(diag(gFITCinfo$m) +
		t(gFITCinfo$PRt) %*% (matrix(f1Hat$eta2 / (1 + gFITCinfo$D *
		f1Hat$eta2), gFITCinfo$n, gFITCinfo$m) * gFITCinfo$PRt))))),
		gFITCinfo$R)

	aNew <- Dnew * f1Hat$eta1
	gammaNew <- t(Rnew) %*% (Rnew %*% (t(Pnew) %*% f1Hat$eta1))

	# We obtain the new marginal means (this is emplained in Miguel Lazaro Gredilla thesis. See the appendix)

	vNew <- as.double(Dnew + (Pnew %*% t(Rnew))^2 %*% rep(1, gFITCinfo$m))
	mNew <- as.double(aNew) + as.double(Pnew %*% gammaNew)

	list(mNew = mNew, vNew = vNew, Rnew = Rnew)
}

#########################################################################################################
#
# The main function that should be used in this module is "epMGPC_internal". You have to call it with the arguments:
#
# 	X  -> Design matrix for the classification problem.
# 	Xbar  -> Matrix of Pseudo inputs.
#       Y  -> Target vector for the classification problem. (must be an R factor)
#       sigma -> variance of the latent noise  (value)
#       sigma0 -> parameter for the gaussian kernel (value)
#       l -> length scale parameter for the gaussian kernel (value)
#       start -> initial EP approximation (NULL if not used) It should be "a".
#
# returns the approximate distribution for the posterior as a list with several components.
#

epMGPCInternal <- function(X, Xbar, Y, sigma = rep(1.0, length(levels(Y))), 
	sigma0 = rep(1.0, length(levels(Y))), l = matrix(1.0, length(levels(Y)), ncol(X)), start = NULL) {

	# This initializes the paramters (otherwise the initialization is applied when the argument is used
	# for the first time) wierd behavior of R

	log_l <- log(l) * 1.0
	log_sigma0 <- log(sigma0) * 1.0
	log_sigma <- log(sigma) * 1.0

	log_l[ log_l < -50 ] <- -50
	log_sigma0[ log_sigma0 < -50 ] <- -50
	log_sigma[ log_sigma < -50 ] <- -50

	log_l[ log_l > 50 ] <- 50
	log_sigma0[ log_sigma0 > 50 ] <- 50
	log_sigma[ log_sigma > 50 ] <- 50

	nK <- length(levels(Y))
	n <- nrow(X)
	d <- ncol(X)
	levelsY <- levels(Y)
	Y <- as.integer(Y)

	# We initialize the approximate factors  reusing any previous solution

	if (!is.null(start))  {
		f1Hat <- start$f1Hat
	} else { 
		f1Hat <- list(eta1k = matrix(0, n, nK), eta1kyi = matrix(0, n, nK), 
			eta2k = matrix(1e-2, n, nK), eta2kyi = matrix(1e-2, n, nK))
	}

	# We initialize the FITC approximations (one for each class)

	gFITCinfo  <- list()

	for (i in 1 : nK) 
		gFITCinfo[[ i ]] <- initialize_kernel_FITC(X, Xbar[[ i ]], sigma[ i ], sigma0[ i ], l[ i, ])

	a <- list(gFITCinfo = gFITCinfo, f1Hat = f1Hat, nK = nK, n = n)

	aOld <- a

	# Main loop of EP

	i <- 1
	damping <- .5
	convergence <- FALSE

	while (! convergence && i < 1e3) {

		update_correct <- FALSE
		damping_inner <- damping
		second_update <- fail <- FALSE

		while(update_correct != TRUE) {

			error <- FALSE

			tryCatch(aNew <<- process_likelihood_factors(a, Y, damping_inner), error = function(x) error <<- TRUE)

			if (error == FALSE) {

				if (fail == TRUE && second_update == FALSE) {
					a <- aNew
					second_update <- TRUE
				} else 
					update_correct <- TRUE

			} else {

				if (i == 1)
					stop("Failure in first EP iteration!")

				cat("Reducing damping factor to guarantee EP update! Damping:", damping_inner, "\n")

				a <- aOld
				damping_inner <- damping_inner * 0.5
				fail <- TRUE
				second_update <- FALSE
			}
		}

		aOld <- a
		a <- aNew

		# We check for convergence

		change <- max(abs(aOld$f1Hat$eta1k - a$f1Hat$eta1k))
		change <- max(change, abs(aOld$f1Hat$eta2k - a$f1Hat$eta2k))
		change <- max(change, abs(aOld$f1Hat$eta1kyi - a$f1Hat$eta1kyi))
		change <- max(change, abs(aOld$f1Hat$eta2kyi - a$f1Hat$eta2kyi))
		
		if (change < 1e-3)
			convergence <- T

		evidence <- computeEvidence(a, Y)
		
		cat("\tIteration",  i, change, "Evidence:", evidence, "\n")
		
		if (REPORT == TRUE) {
		  
		  t_before <- proc.time()
		  ret <- 	list(f1Hat = a$f1Hat, gFITCinfo = a$gFITCinfo, n = n, nK = nK, Y = Y, a = a, levelsY = levelsY)
		  prediction <- predictMGPC(ret, X_test)
		  ll <- -mean(log(prediction$prob[ cbind(1:nrow(X_test), as.integer(Y_test))]))
		  ee <- mean(as.integer(prediction$labels) != Y_test) 
		  t_after <- proc.time()
		  
		  t0 <<- t0 + (t_after - t_before)
		  
		  write.table(t(c(ee, ll, proc.time() - t0)), 
		              file = paste("./results/time_outter_", CONT, ".txt", sep = ""), row.names = F, col.names = F, append = TRUE)
		}

		# Annealed damping scheme

		damping <- damping * 0.99

		i <- i + 1
	}

	# We compute the evidence and its gradient

	gradientLogZ <- list()

	for (i in 1 : nK) {

		eta1 <- a$f1Hat$eta1k[ , i ]
		eta2 <- a$f1Hat$eta2k[ , i ]
		eta1[ Y == i ] <- eta1[ Y == i ] + rowSums(a$f1Hat$eta1kyi[ Y == i ,, drop = FALSE])
		eta2[ Y == i ] <- eta2[ Y == i ] + rowSums(a$f1Hat$eta2kyi[ Y == i ,, drop = FALSE])

		gradientLogZ[[ i ]] <- computeGrads(a$gFITCinfo[[ i ]], l[ i, ], sigma0[ i ], sigma[ i ], eta2, eta1) 
	}

	# We compute the evidence

	logZ <- computeEvidence(a, Y)

	# We are done!

	list(f1Hat = a$f1Hat, logZ = logZ, gradientLogZ = gradientLogZ, gFITCinfo = a$gFITCinfo, 
		n = n, nK = nK, Y = Y, a = a, levelsY = levelsY)
}

###
# Function which processes the likelihood
#
# @param	a 		The approximation 
# @param	Y		The class labels.
# @param	damping		The damping factor 
#
# @return	a 		The updated approximation
#
#

process_likelihood_factors <- function(a, Y, damping) {

	# We get the marginals of the posterior approximation

	v_kyi <- v_k <- matrix(0, a$n, a$nK)
	m_kyi <- m_k <- matrix(0, a$n, a$nK)

	for (i in 1 : a$nK) {

		eta1 <- a$f1Hat$eta1k[ , i ]
		eta2 <- a$f1Hat$eta2k[ , i ]
		eta1[ Y == i ] <- eta1[ Y == i ] + rowSums(a$f1Hat$eta1kyi[ Y == i ,, drop = FALSE])
		eta2[ Y == i ] <- eta2[ Y == i ] + rowSums(a$f1Hat$eta2kyi[ Y == i ,, drop = FALSE])

		f1Hat <- list(eta1 = eta1, eta2 = eta2)

		ret <- reconstructPosterior(a$gFITCinfo[[ i ]], f1Hat)

		v_k[ , i ] <- ret$vNew
		m_k[ , i ] <- ret$mNew
		v_kyi[ Y == i, ] <- matrix(ret$vNew[ Y == i ], sum(Y == i), a$nK)
		m_kyi[ Y == i, ] <- matrix(ret$mNew[ Y == i ], sum(Y == i), a$nK)
	}

	v_kyi_old <- (v_kyi^-1 - a$f1Hat$eta2kyi)^-1
	m_kyi_old <- v_kyi_old * (m_kyi / v_kyi - a$f1Hat$eta1kyi)
	v_k_old <- (v_k^-1 - a$f1Hat$eta2k)^-1
	m_k_old <- v_k_old * (m_k / v_k - a$f1Hat$eta1k)


	if (any(v_kyi_old < 0) || any(v_k_old < 0)) {
		stop("Negative variances!")
	}

	logZ <- pnorm((m_kyi_old - m_k_old) / sqrt(v_k_old + v_kyi_old), log.p = TRUE)

	ratio <- exp(-logZ + dnorm((m_kyi_old - m_k_old) / sqrt(v_kyi_old + v_k_old), log = T))

	# Alpha and Beta are the first and second derivatives of logZ with respect to the mean of the cavity

	alpha_k <- -ratio / sqrt(v_k_old + v_kyi_old)
	alpha_kyi <- ratio / sqrt(v_k_old + v_kyi_old)
	beta_k <- - ratio * ((m_kyi_old - m_k_old) / sqrt(v_k_old + v_kyi_old) + ratio) / (v_k_old + v_kyi_old)
	beta_kyi <- - ratio * ((m_kyi_old - m_k_old) / sqrt(v_k_old + v_kyi_old) + ratio) / (v_k_old + v_kyi_old)

	eta2HatNew_k <- -beta_k / (1 + beta_k * v_k_old)
	eta1HatNew_k <- (alpha_k - m_k_old * beta_k) / (1 + beta_k * v_k_old)
	eta2HatNew_kyi <- -beta_kyi / (1 + beta_kyi * v_kyi_old)
	eta1HatNew_kyi <- (alpha_kyi - m_kyi_old * beta_kyi) / (1 + beta_kyi * v_kyi_old)

	# We set to uniform those factors that do not actually appear in the likelihood

	eta1HatNew_k[ cbind(1 : a$n, Y) ] <- 0
	eta2HatNew_k[ cbind(1 : a$n, Y) ] <- 0
	eta1HatNew_kyi[ cbind(1 : a$n, Y) ] <- 0
	eta2HatNew_kyi[ cbind(1 : a$n, Y) ] <- 0

	# We do the dammped updates

	a$f1Hat$eta2k <- damping * eta2HatNew_k + (1 - damping) * a$f1Hat$eta2k
	a$f1Hat$eta1k <- damping * eta1HatNew_k + (1 - damping) * a$f1Hat$eta1k

	a$f1Hat$eta2kyi <- damping * eta2HatNew_kyi + (1 - damping) * a$f1Hat$eta2kyi
	a$f1Hat$eta1kyi <- damping * eta1HatNew_kyi + (1 - damping) * a$f1Hat$eta1kyi

	a
}

###
# Function which computes the EP approximation of the log evidence.
#
# @param	f1Hat		The approximation for the first factor.
# @param	gFITCinfo	The list with the problem information.
# @param	Y		The class labels.
#
# @return	logZ		The log evidence.
#

computeEvidence <- function(a, Y) {

	# We get the marginals of the posterior approximation

	v_kyi <- v_k <- matrix(0, a$n, a$nK)
	m_kyi <- m_k <- matrix(0, a$n, a$nK)
	Rnew <- list()

	for (i in 1 : a$nK) {

		eta1 <- a$f1Hat$eta1k[ , i ]
		eta2 <- a$f1Hat$eta2k[ , i ]
		eta1[ Y == i ] <- eta1[ Y == i ] + rowSums(a$f1Hat$eta1kyi[ Y == i ,, drop = FALSE])
		eta2[ Y == i ] <- eta2[ Y == i ] + rowSums(a$f1Hat$eta2kyi[ Y == i ,, drop = FALSE])

		f1Hat <- list(eta1 = eta1, eta2 = eta2)

		ret <- reconstructPosterior(a$gFITCinfo[[ i ]], f1Hat)

		Rnew[[ i ]] <- ret$Rnew
		v_k[ , i ] <- ret$vNew
		m_k[ , i ] <- ret$mNew
		v_kyi[ Y == i, ] <- matrix(ret$vNew[ Y == i ], sum(Y == i), a$nK)
		m_kyi[ Y == i, ] <- matrix(ret$mNew[ Y == i ], sum(Y == i), a$nK)
	}

	v_kyi_old <- (v_kyi^-1 - a$f1Hat$eta2kyi)^-1
	m_kyi_old <- v_kyi_old * (m_kyi / v_kyi - a$f1Hat$eta1kyi)
	v_k_old <- (v_k^-1 - a$f1Hat$eta2k)^-1
	m_k_old <- v_k_old * (m_k / v_k  - a$f1Hat$eta1k)

	if (any(v_kyi_old < 0) || any(v_k_old < 0)) {
		stop("Negative variances in the computation of the evidence!")
	}

	# We assume approximate factors of the form exp(eta1 * x - 0.5 * eta2 * x^2)

	logZ <- pnorm((m_kyi_old - m_k_old) / sqrt(v_k_old + v_kyi_old), log.p = TRUE)

	to_add <- logZ - (0.5 * log(v_k) - 0.5 * log(v_k_old) - 0.5 * m_k_old^2 / v_k_old + 0.5 * m_k^2 / v_k) - 
		(0.5 * log(v_kyi) - 0.5 * log(v_kyi_old) - 0.5 * m_kyi_old^2 / v_kyi_old + 0.5 * m_kyi^2 / v_kyi)

	to_add[ cbind(1 : a$n, Y) ] <- 0

	logZret <- sum(to_add)

	for (i in 1 : a$nK) {

		eta1 <- a$f1Hat$eta1k[ , i ]
		eta2 <- a$f1Hat$eta2k[ , i ]
		eta1[ Y == i ] <- eta1[ Y == i ] + rowSums(a$f1Hat$eta1kyi[ Y == i ,, drop = FALSE])
		eta2[ Y == i ] <- eta2[ Y == i ] + rowSums(a$f1Hat$eta2kyi[ Y == i ,, drop = FALSE])

		# We add the difference between the log normalizer of the posterior and the prior 
		# This uses the matrix determinant lemma and the particular form of Rnew to cancel some of the
		# terms that appear

		logZret <- logZret + (sum(log(diag(Rnew[[ i ]]))) - sum(log(diag(a$gFITCinfo[[ i ]]$R)))) - 
			0.5 * sum(log(1 + eta2 * a$gFITCinfo[[ i ]]$D)) + 0.5 * sum(eta1 * m_k[, i ]) 
	}

	logZret
}

##
# This function computes the gradients of the ML approximation provided by EP once it has converged
#

computeGrads <- function(gFITCinfo, l, sigma0, sigma, tauTilde, muTilde) {

	# The gradient is given by Eq. (34) in http://arxiv.org/pdf/1507.04513.pdf

	# We have to take into account the form of the FITC prior, that is, D + P R^T R P^T

	# M = Kprior^-1 - Kprior^-1 Kpost Kprior^-1 - Kprior^-1 meanPost meanPost^T Kprior^-1
	# We compute first the representation of M = (K + PI^-1)^-1 - (PI^-1 + K)^-1 PI^-1 eta1 eta1^T PI^-1 (PI^-1 + K)^-1
	# (K + PI^-1)^-1 comes from the woodbury formula
	# and (K + PI^-1)^-1 = Da - PaRatRaPat

	# We use the chain rule of matrix derivatives which is the trace of (M * dKdparam) 
	# It also uses the fact that the gradient of the kernel can be written as tr(M * dKdtheta) and
	# dKdtheta has some structure. For example, dKdtheta = K o xy^T where o is the Hadamard product.
	# The trace satisfies then that tr(M (K o x y^T)) = x^T (M^T o K) y. This allows to write all gradients
	# using matrix multiplications.

	Da <- 1 / (gFITCinfo$D + tauTilde^-1)
	Pa <- matrix(Da, gFITCinfo$n, gFITCinfo$m) * gFITCinfo$PRt
	Ra <- cholInverse(diag(gFITCinfo$m) + t(gFITCinfo$PRt) %*% Pa)
	PaRat <- Pa %*% t(Ra)
	RtRP <- t(gFITCinfo$PRt %*% gFITCinfo$R)

	upsilon <- Da * muTilde / tauTilde - PaRat %*% (t.default(PaRat) %*% (muTilde / tauTilde))

	diagM <- Da - rowSums(PaRat^2) - upsilon^2
	RtRPM <- matrix(Da, gFITCinfo$m, gFITCinfo$n, byrow = TRUE) * RtRP - (RtRP %*% PaRat) %*% t(PaRat) - (RtRP %*% upsilon) %*% t(upsilon)
	RtRPMPtRRt <- RtRPM %*% t(RtRP)
	RtRPdiagM <- RtRP * matrix(diagM, gFITCinfo$m, gFITCinfo$n, byrow = TRUE)
	RtRPdiagMPtRRt <- RtRPdiagM %*% t(RtRP)

	# We compute the derivatives of the kernel with respect to log_sigma

	dKnn_dlog_sigma0 <- rep(sigma0, gFITCinfo$n)
	dKmm_dlog_sigma0 <- diag(gFITCinfo$m) * sigma0
	dP_dlog_sigma0 <- matrix(0, gFITCinfo$n, gFITCinfo$m)

	gr_log_sigma0 <- -0.5 * sum(diagM * dKnn_dlog_sigma0) + sum(RtRPdiagM * t(dP_dlog_sigma0)) - 
		0.5 * sum(RtRPdiagMPtRRt * dKmm_dlog_sigma0) - sum(RtRPM * t(dP_dlog_sigma0)) + 0.5 * sum(RtRPMPtRRt * dKmm_dlog_sigma0)

	# We compute the derivatives of the kernel with respect to log_sigma0

	dKnn_dlog_sigma <- rep(sigma, gFITCinfo$n)
	dKmm_dlog_sigma <- gFITCinfo$Kmm - diag(rep(gFITCinfo$sigma0, gFITCinfo$m)) - diag(1e-10, gFITCinfo$m)
	dP_dlog_sigma <- gFITCinfo$P 

	gr_log_sigma <- -0.5 * sum(diagM * dKnn_dlog_sigma) + sum(RtRPdiagM * t(dP_dlog_sigma)) - 
		0.5 * sum(RtRPdiagMPtRRt * dKmm_dlog_sigma) - sum(RtRPM * t(dP_dlog_sigma)) + 0.5 * sum(RtRPMPtRRt * dKmm_dlog_sigma)

	# We compute the derivatives of the kernel with respect to l

	gr_log_l <- rep(0, length(l))

#	for (i in 1 : length(l)) {
#
#		dKnn_dlog_l <- rep(0, gFITCinfo$n)
#
#		distance <- as.matrix(dist(gFITCinfo$Xbar[, i, drop = FALSE ]))^2
#		dKmm_dlog_l <- gFITCinfo$Kmm * 0.5 * distance / l[ i ]
#
#		Q <- matrix(gFITCinfo$X[ , i ]^2, gFITCinfo$n, gFITCinfo$m)
#		Qbar <- matrix(gFITCinfo$Xbar[ , i ]^2, gFITCinfo$n, gFITCinfo$m, byrow = T)
#		distance <- Qbar + Q - 2 * gFITCinfo$X[ , i ] %*% t(gFITCinfo$Xbar[ , i ])
#
#                dP_dlog_l <- gFITCinfo$P * 0.5 * distance / l[ i ]

##		gr_log_l[ i ] <- -0.5 * sum(diagM * dKnn_dlog_l) + sum(RtRPdiagM * t(dP_dlog_l)) - 
##			0.5 * sum(RtRPdiagMPtRRt * dKmm_dlog_l) - sum(RtRPM * t(dP_dlog_l)) + 0.5 * sum(RtRPMPtRRt * dKmm_dlog_l)

#		gr_log_l[ i ] <- sum(RtRPdiagM * t(dP_dlog_l)) - 0.5 * sum(RtRPdiagMPtRRt * dKmm_dlog_l) - 
#			sum(RtRPM * t(dP_dlog_l)) + 0.5 * sum(RtRPMPtRRt * dKmm_dlog_l)
#	}

	Ml <- 0.5 * (t(RtRPdiagM) - t(RtRPM)) * gFITCinfo$P
	Xbarl <-  (gFITCinfo$Xbar / matrix(sqrt(l), nrow(gFITCinfo$Xbar), ncol(gFITCinfo$Xbar), byrow = TRUE))
	Xl <-  (gFITCinfo$X / matrix(sqrt(l), nrow(gFITCinfo$X), ncol(gFITCinfo$X), byrow = TRUE))
	Ml2 <- (-0.5 * RtRPdiagMPtRRt + 0.5 * RtRPMPtRRt) * gFITCinfo$Kmm * 0.5

	gr_log_l <- ((rep(1, ncol(Ml)) %*% t(Ml)) %*% Xl^2) - 2  * colSums(Xl * (Ml %*% Xbarl)) +  ((rep(1, nrow(Ml)) %*% Ml) %*% Xbarl^2) +
		((rep(1, ncol(Ml2)) %*% t(Ml2)) %*% Xbarl^2) - 2 * colSums(Xbarl * (Ml2 %*% Xbarl)) + ((rep(1, nrow(Ml2)) %*% Ml2) %*% Xbarl^2)

#	gr_log_l <- colSums(t(Ml) %*% Xl^2) - 2  * colSums(Xl * (Ml %*% Xbarl)) +  colSums(Ml %*% Xbarl^2) +
#		colSums(t(Ml2) %*% Xbarl^2) - 2 * colSums(Xbarl * (Ml2 %*% Xbarl)) + colSums(Ml2 %*% Xbarl^2)

	gr_xbar <- matrix(0, gFITCinfo$m, length(l)) 

#	for (i in 1 : length(l)) {
#
#		distance <- matrix(gFITCinfo$X[ , i ], gFITCinfo$n, gFITCinfo$m) - matrix(gFITCinfo$Xbar[ , i ], 
#			gFITCinfo$n, gFITCinfo$m, byrow = T)
#		dP_dXbar <- gFITCinfo$P * distance / gFITCinfo$l[ i ]
#
#		distance <- (matrix(gFITCinfo$Xbar[ , i ], gFITCinfo$m, gFITCinfo$m) - matrix(gFITCinfo$Xbar[ , i ],
#			gFITCinfo$m, gFITCinfo$m, byrow = T))
#		dKmm_dXbar <- - 2 * (gFITCinfo$Kmm - diag(gFITCinfo$sigma0, gFITCinfo$m) - diag(1e-10, gFITCinfo$m)) * distance / gFITCinfo$l[ i ]
#
#		dKnn_dXbar <- rep(0, gFITCinfo$n)
#
##		gr_xbar[ ,i ] <- -0.5 * rep(sum(diagM * dKnn_dXbar), gFITCinfo$m) + rowSums(RtRPdiagM * t(dP_dXbar)) - 
##			0.5 * rowSums(RtRPdiagMPtRRt * dKmm_dXbar) - rowSums(RtRPM * t(dP_dXbar)) + 0.5 * rowSums(RtRPMPtRRt * dKmm_dXbar)
#
#		gr_xbar[ ,i ] <- rowSums(RtRPdiagM * t(dP_dXbar)) - 0.5 * rowSums(RtRPdiagMPtRRt * dKmm_dXbar) - 
#			rowSums(RtRPM * t(dP_dXbar)) + 0.5 * rowSums(RtRPMPtRRt * dKmm_dXbar)
#	}

	Xbar <- (gFITCinfo$Xbar / matrix(l, nrow(gFITCinfo$Xbar), ncol(gFITCinfo$Xbar), byrow = TRUE))
        X <- (gFITCinfo$X / matrix(l, nrow(gFITCinfo$X), ncol(gFITCinfo$X), byrow = TRUE))
	M <- - 0.5 * ((RtRPdiagMPtRRt - RtRPMPtRRt) * - 2 * (gFITCinfo$Kmm - diag(gFITCinfo$sigma0, gFITCinfo$m) - diag(1e-10, gFITCinfo$m)))
	P <- (t(RtRPdiagM) - t(RtRPM)) * gFITCinfo$P
	gr_xbar <- (Xbar * matrix(rep(1, gFITCinfo$m) %*% M, gFITCinfo$m, length(l))) - M %*% Xbar + t(P) %*% X - 
		(Xbar * matrix(rep(1, gFITCinfo$n) %*% P, gFITCinfo$m, length(l)))

	list(gr_log_l = gr_log_l, gr_log_sigma0 = gr_log_sigma0, gr_log_sigma = gr_log_sigma, gr_xbar = gr_xbar)
}


###
# Function which computes the probability of class 1 on new data.
#
# @param	ret	The list returned by epGPCExternal
# @param	Xtest	The n x d matrix with the new data points.
#
# @return	pOne	The probability of class 1 on the new data.
#

predictMGPC <- function(ret, Xtest) {

	n <- nrow(Xtest)
	v_k <- v_k <- matrix(0, n, ret$nK)
	m_k <- m_k <- matrix(0, n, ret$nK)

	# To make prediction we need the means and the variances of the process at the new point.
	# These are obtained by computing int p(f*|f) p(f|D) d f, where p(f|D) is the posterior
	# approximation computed. Importantly Kf*f (the covariance between the observed points and the
	# new points) is the approxmate covariance used by the FITC approximation.

	# The integrals needed for making predictions are given by Eq. (19) in 
	# http://papers.nips.cc/paper/4241-robust-multi-class-gaussian-process-classification.pdf

	for (i in 1 : ret$nK) {

		eta1 <- ret$f1Hat$eta1k[ , i ]
		eta2 <- ret$f1Hat$eta2k[ , i ]
		eta1[ ret$Y == i ] <- eta1[ ret$Y == i ] + rowSums(ret$f1Hat$eta1kyi[ ret$Y == i ,, drop = FALSE])
		eta2[ ret$Y == i ] <- eta2[ ret$Y == i ] + rowSums(ret$f1Hat$eta2kyi[ ret$Y == i ,, drop = FALSE])

		muTilde <- eta1
		tauTilde <- eta2

		gFITCinfo <- ret$gFITCinfo[[ i ]]

		pStar <- kernel_nm(Xtest, gFITCinfo$Xbar, gFITCinfo$l, gFITCinfo$sigma)
		dStar <- computeDiagKernel(Xtest, gFITCinfo$sigma, gFITCinfo$sigma0) - (pStar %*%t(gFITCinfo$R))^2 %*% rep(1, gFITCinfo$m)
	
		Dnew <- gFITCinfo$D / (1 + gFITCinfo$D * tauTilde)
		Pnew <- matrix(1 / (1 + gFITCinfo$D * tauTilde), gFITCinfo$n, gFITCinfo$m) * gFITCinfo$P
	
		Rnew <- backsolve(rot180(t(chol(rot180(diag(gFITCinfo$m) + t(gFITCinfo$PRt) %*% (matrix(tauTilde / (1 + gFITCinfo$D *
			tauTilde), gFITCinfo$n, gFITCinfo$m) * gFITCinfo$PRt))))), gFITCinfo$R)
	
		gammaNew <- t(Rnew) %*% (Rnew %*% (t(Pnew) %*% muTilde))
	
		# We obtain the predictive means and variances
	
		mPrediction <- pStar %*% gammaNew 
		vPrediction <- dStar + (pStar %*% t(Rnew))^2 %*% rep(1, gFITCinfo$m)

		m_k[ , i ] <- mPrediction
		v_k[ , i ] <- vPrediction
	}

	# We compute the label of each instance using quadrature

	labels <- rep(ret$Y[ 1 ], n)
	prob <- matrix(0, n, length(ret$levelsY))
	colnames(prob) <- ret$levelsY

	for (i in 1 : n)  {

		maxProb <- 0

		for (k in 1 : length(ret$levelsY)) {

#			pTmp <- integrate(function(u) 
#				apply(as.matrix(u), 1, function(u_tmp) {
#					exp(sum(pnorm((u_tmp - m_k[ i, -k ]) / sqrt(v_k[ i, -k ]), log.p = TRUE)) + 
#					dnorm(u_tmp, mean = m_k[ i, k ], sd = sqrt(v_k[ i, k ]), log = TRUE))
#				}) , -Inf, +Inf)$value

			f <- function(u) { apply(as.matrix(u), 1, function(u_tmp) {
                                       exp(sum(pnorm((u_tmp - m_k[ i, -k ]) / sqrt(v_k[ i, -k ]), log.p = TRUE)) + 
                                       dnorm(u_tmp, mean = m_k[ i, k ], sd = sqrt(v_k[ i, k ]), log = TRUE))
                               })}

			grid <- seq(-6, 6, length = 100) * sqrt(v_k[ i, k ]) + m_k[ i, k ]
			values <- f(grid)
			pTmp <- sum(values * (grid[ 2 ] - grid[ 1 ]))

			prob[ i , k ] <- pTmp

			if (maxProb < pTmp) {
				labels[ i ] <- ret$levelsY[ k ]
				maxProb <- pTmp
			}
		}
	}

	list(labels = labels, prob = prob, maxprob = apply(prob, 1, max))
}

##
# Function which adjusts the kernel parameters by gradient descent.
#
# @param	X	n x d matrix with the data points.
# @param	Y	n-dimensional vector with the class labels.
# @param	m	The number of pseudo-inputs to use.
#
# @return	ret	A list with the following elements:
#
#			a		The posterior approximation.
#			f1Hat		The approximation to the first factor.
#			f2Har		The approximation to the second factor.
#			logZ		The approximation to the log evidence.
#			gradientLogZ	The gradient of the log evidence.
#
#			optimize_flags: booleans corresponding to sigma, sigma0, l, pesudoinputs
#					indicates whether that parameter should be optimized or not
#

epMGPCExternal_gradient_ascent <- function(X, Y, m, sigma = rep(1, length(levels(Y))), Xbar_ini = NULL,
	sigma0 = rep(1, length(levels(Y))), l = matrix(1, length(levels(Y)), ncol(X)), max_iters = 250) {

	# We initialize the hyper-parameters (We take the log)

	sigma <- log(sigma)
	sigma0 <- log(sigma0)
	l <- log(l)

	Xbar <- Xbar_ini

	if (is.null(Xbar)) {
		Xbar <- list()
		for (i in 1 : length(levels(Y))) 
			Xbar[[ i ]] <- X[ sample(1 : nrow(X),  m), ]
	}

	# We initialize the gradient optimization process

	ret <- epMGPCInternal(X, Xbar, Y, exp(sigma), exp(sigma0), exp(l))

	best <- ret

	# cat(0, "New evidence:", ret$logZ, "\n")

	eps <- 0.001
	convergence <- F
	iteration <- 1
	
	while (!convergence && iteration < max_iters) {

		# We update the parameters

		for (i in 1 : length(levels(Y))) {
			sigma[ i ] <- sigma[ i ] + eps * ret$gradientLogZ[[ i ]]$gr_log_sigma
			sigma0[ i ] <- sigma0[ i ] + eps * ret$gradientLogZ[[ i ]]$gr_log_sigma0
			l[ i, ] <- l[ i, ] + eps * ret$gradientLogZ[[ i ]]$gr_log_l
			Xbar[[ i ]] <- Xbar[[ i ]] + eps * ret$gradientLogZ[[ i ]]$gr_xbar
		}

		# We train the model using the previous solution as the starting point. If that fails (the starting point
		# is unfeasible we start from scracht)

		tryCatch( retNew <<- epMGPCInternal(X, Xbar, Y, exp(sigma), exp(sigma0), exp(l), ret$a) , error = function(x)  {
			tryCatch( retNew <<- epMGPCInternal(X, Xbar, Y, exp(sigma), exp(sigma0), exp(l))
                	, error = function(x)  convergence <<- TRUE)
			}
                )

		if (is.nan(retNew$logZ))
			return(ret)

		if (abs(retNew$logZ - ret$logZ) < 1e-4)
			convergence <- T

		if (retNew$logZ < ret$logZ)
			eps <- eps * 0.5
		else
			eps <- eps * 1.1

		#cat(iteration, "New evidence:", retNew$logZ, "eps:", eps, "Change:", abs(retNew$logZ - ret$logZ), "\n")

		ret <- retNew

		if (ret$logZ > best$logZ)
			best <- ret

		iteration <- iteration + 1
	}

	best
}

##
# This trains the model using LBFGS. It reuses the last EP solutions.
#

epMGPCExternal_lBFGS <- function(X, Y, m, sigma = rep(1, length(levels(Y))), Xbar_ini = NULL,
	sigma0 = rep(1, length(levels(Y))), l = matrix(1, length(levels(Y)), ncol(X)), max_iters = 250) {
  
        t0 <<- proc.time()
  
	# We initialize the hyper-parameters (We take the log)

	sigma <- log(sigma)
	sigma0 <- log(sigma0)
	l <- log(l)

	Xbar <- Xbar_ini

	if (! is.null(Xbar)) {
		Xbar <- list()
		for (i in 1 : length(levels(Y))) 
			Xbar[[ i ]] <- X[ sample(1 : nrow(X),  m), ]
	}

	# We initialize the gradient optimization process

	epSolution <- epMGPCInternal(X, Xbar, Y, exp(sigma), exp(sigma0), exp(l))

	# These functions codify and decofify the hyper-parameters as an array of doubles

	codify_params <- function(Xbar, sigma, sigma0, l) {
		
		params <- NULL

		params <- c(params, sigma, sigma0, l)

		for (i in 1 : length(Xbar)) {
			params <- c(params, Xbar[[ i ]])
		}

		params
	}

	decodify_params <- function(params) {

		sigma <- params[ 1 : length(sigma) ]
		params <- params[ -(1 : length(sigma)) ]

		sigma0 <- params[ 1 : length(sigma0) ]
		params <- params[ -(1 : length(sigma0)) ]

		l <- matrix(params[ 1 : length(c(l)) ], nrow(l), ncol(l))
		params <- params[ -(1 : length(c(l))) ]

		retXbar <- list()

		for (i in 1 : length(Xbar)) {
			retXbar[[ i ]] <- matrix(params[ 1 : length(c(Xbar[[ i ]])) ], nrow(Xbar[[ i ]]), ncol(Xbar[[ i ]]))
			params <- params[ -(1 : length(c(Xbar[[ i ]]))) ]
		}

		list(sigma = sigma, sigma0 = sigma0, l = l, Xbar = retXbar)
	}

	# The next functions evaluate the objective and the gradient

	evaluate_objective <- function(params) {
		
		result <- decodify_params(params)

		# We train the model using the previous solution as the starting point. If that fails (the starting point
		# is unfeasible we start from scracht)
		
		tryCatch( epSolution <<- epMGPCInternal(X, result$Xbar, Y, exp(result$sigma), exp(result$sigma0), exp(result$l), 
			epSolution$a) , error = function(x)  {
			tryCatch(epSolution <<- epMGPCInternal(X, result$Xbar, Y, exp(result$sigma), exp(result$sigma0), exp(result$l))
	                	, error = function(x) {browser(); stop("Failure inside LBFGS when evaluating the objective: ",x)})
			}
                )

		- 1.0 * epSolution$logZ
	}

	evaluate_gradient <- function(params) {

		result <- decodify_params(params)
		sigma <- result$sigma
		sigma0 <- result$sigma0
		l <- result$l
		Xbar <- result$Xbar

		# We train the model using the previous solution as the starting point. If that fails (the starting point
		# is unfeasible we start from scracht)

		tryCatch(epSolution <<- epMGPCInternal(X, result$Xbar, Y, exp(result$sigma), exp(result$sigma0), 
			exp(result$l), epSolution$a) , error = function(x)  {
			tryCatch(epSolution <<- epMGPCInternal(X, result$Xbar, Y, exp(result$sigma), exp(result$sigma0), exp(result$l))
	                	, error = function(x) {browser(); stop("Failure inside LBFGS when evaluating the gradient") })
		})

		# We obtain the gradient and codify it

		gradient <- NULL

		for (i in 1 : length(epSolution$gradientLogZ))
			gradient <- c(gradient, epSolution$gradientLogZ[[ i ]]$gr_log_sigma)

		for (i in 1 : length(epSolution$gradientLogZ))
			gradient <- c(gradient, epSolution$gradientLogZ[[ i ]]$gr_log_sigma0)

		gr_l <- matrix(0, nrow(l), ncol(l))
		for (i in 1 : length(epSolution$gradientLogZ))
			gr_l[ i, ] <- epSolution$gradientLogZ[[ i ]]$gr_log_l

		gradient <- c(gradient, gr_l)

		for (i in 1 : length(epSolution$gradientLogZ))
			gradient <- c(gradient, epSolution$gradientLogZ[[ i ]]$gr_xbar)


		- 1.0 * gradient
	}

	# Now we optimize using LBFGS

	starting_point <- codify_params(Xbar, sigma, sigma0, l)

	result <- optim(starting_point, evaluate_objective, evaluate_gradient, method = "L-BFGS-B", 
		control = list(trace = 1, REPORT = 1, maxit = max_iters))

	epSolution
}



