#' multivariate ecdfs
#' 
#' The function mecdf, maps a data matrix to an object of class "mecdf",
#' representing multivariate ECDF models. In the mecdf package, such models are
#' also functions. In principle, the mecdf function's first argument is a
#' numeric matrix, representing multivariate data. Each row, represents one
#' multivariate realisation (we can think of this as either, the realised value
#' of a vector random variable or the realised values of multiple random
#' variables). Each column, represents one (random) variable. The mecdf
#' function's second argument determines whether or not the model is a step
#' function or a continuous function. As of mecdf 0.6.0, all models are step
#' functions by default.
#' 
#' 
#' @aliases mecdf print.mecdf plot.mecdf
#' @param m An object of class "mecdf".
#' @param x In principle, a numeric matrix representing multivariate
#' realisations, where each row represents one realisation and each column
#' represents one variable. A vector can also be used, in which case its
#' converted to a matrix, using cbind.
#' @param continuous Logical (defaults to false), if true a continuous
#' function, otherwise a step function.
#' @param validate Logical (defaults to true), whether or not to validate the
#' mecdf arguments. In general, this should be true.
#' @param expand Logical (default, true if continuous, false if step), whether
#' or not to correct the model, by adding two extra points.
#' @param project Logical (defaults to false), if true, the matrix is
#' transformed, such that the function's values fall within a unit cube with
#' uniform marginals.
#' @param expandf Numeric (defaults to 0.1), giving the expansion factor.
#' Ignored, if expand is false.
#' @param \dots .
#' @return The model that's returned (reiterating the model's a function), maps
#' a data matrix to a vector of cumulative probabilities. In principle, the
#' function's only argument is a numeric matrix, with the same row and column
#' conventions as above. It maps each row of the matrix to one value in the
#' vector. Note that a regular vector can be used as an argument, however it's
#' meaning is ambiguous. If the model is univariate, vector arguments are
#' equivalent to matrices with one column. If the model is multivariate, vector
#' arguments, are equivalent to a matrices with one row.
#' @export mecdf
mecdf = function (x, continuous=FALSE, ...,
	validate=TRUE, expand=continuous, project=FALSE, expandf=0.1)
{	x = cbind (x)
	nraw = nr = nrow (x)
	nc = ncol (x)
	if (validate)
	{	if (length (list (...) ) > 0)
			stop ("invalid constructor argument")
		if (!is.numeric (x) ) stop ("x must be numeric")
		if (!all (is.finite (x) ) ) stop ("all x must be finite")
		for (j in 1:nc) if (length (unique (x [,j]) ) < 2)
			stop ("each variable requires at least 2 distinct realisations")
		if (nc == 1) x [] = sort (x)
		if (is.null (colnames (x) ) ) colnames (x) = paste ("x", 1:ncol (x), sep="")
		if (is.null (rownames (x) ) ) rownames (x) = 1:nr
	}
	if (expand)
	{	nr = nr + 2
		a = b = numeric (nc)
		for (j in 1:nc)
		{	xrng = range (x [,j])
			xf = expandf * diff (xrng)
			a [j] = xrng [1] - xf
			b [j] = xrng [2] + xf
		}
		x = rbind (a, x, b)
	}
	if (project)
		for (j in 1:nc) x [,j] = (order (order (x [,j]) ) - 1) / (nr - 1)
	Fh = Fst = NULL
	if (nc > 1)
	{	if (continuous)
		{	Fh = .mecdf.continuous
			Fst = .mecdf.vertex
		}
		else Fh = FUNCTION (.mecdf.step)
	}
	else
	{	if (continuous) Fh = .uecdf.continuous
		else Fh =.uecdf.step
	}
	f <- extend(FUNCTION(.mecdf.main), "mecdf")
	f$continuous = continuous
	f$Fh = Fh
	f$Fst = Fst
	f$nraw = nraw
	f$nr = nr
	f$nc = nc
	f$x = x
	f
}

.mecdf.main = function (u)
{	if (.$nc > 1)
	{	if (!is.matrix (u) ) u = rbind (u)
		if (.$nc != ncol (u) )
			stop ("k-variate mecdf requires k-column matrix")
		.mecdf.interpolate (.$Fh, .$Fst, .$nr, .$nc, .$x, u)
	}
	else
	{	if (is.matrix (u) && ncol (u) > 1)
			stop ("univariate mecdf doesn't accept multicolumn matrix")
		.uecdf.interpolate (.$Fh, .$nr, .$x, u)
	}
}

print.mecdf = function (m, ...)
{	variate = if (m$nc == 1) "univariate"
	else if (m$nc == 2) "bivariate"
	else paste (m$nc, "-variate", sep="")
	type = if (m$continuous) "continuous" else "step"
	cat ("mecdf_{", variate, ", ", type, "}\n", sep="")
	#print ( (m$x) )
}

plot.mecdf = function (m, ...)
{	p = m (m$x)
	if (m$nc == 1) .uecdf.plot (m, p, m$continuous, ...)
	else if (m$nc == 2) .becdf.plot (m, p, ...)
	else stop ("s3x_plot.mecdf only supports univariate and bivariate models")
}

.uecdf.plot = function (e, p, continuous, ...)
{	xlab = colnames (e$x)
	ylab = "Fh(x)"
	if (continuous)
		plot (e$x, p, ylim=c (0, 1), yaxs="i", type="l", xlab=xlab, ylab=ylab, ...)
	else
	{	plot (e$x, p, ylim=c (0, 1), yaxs="i", xlab=xlab, ylab=ylab, pch=NA, ...)
		x1 = e$x [-e$nr]
		x2 = e$x [-1]
		p0 = p [-e$nr]
		segments (x1, p0, x2, p0)
		segments (e$x, c (0, p), e$x, c (p, 1) )
	}
}

.becdf.plot = function (e, p, lines=TRUE, lty=1, col=rgb (0.975, 0.7, 0), ...)
{	labs = colnames (e$x)
	x1 = e$x [,1]; x2 = e$x [,2]
	plot (x1, x2, xlab=labs [1], ylab=labs [2], pch=NA, ...)
	if (lines)
	{	segments (x1, x2, x1 - 2 * diff (range (x1) ), x2, lty=lty, col=col)
		segments (x1, x2, x1, x2 - 2 * diff (range (x2) ), lty=lty, col=col)
	}
	text (x1, x2, round (p, 2) )
}



