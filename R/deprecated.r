bcdfplot = function (...) plotbcdf (...)
bcdf.plot = function (...) plotbcdf (...)

mecdf.theme = function (theme)
{	#for .First.lib
	rgb = grDevices::rgb

	if (theme == "simple") .mecdfcols (fline="black", ffill="white")
	else if (theme == "emerald")
		.mecdfcols (rgb (0.08, 0.6, 0.4), c (0, 0.4, 0.05), c (0.2, 1, 0.1) )
	else if (theme == "blue-sapphire")
		.mecdfcols ("blue", c (0.3, 0, 0.6), c (0, 0, 1) )
	else if (theme == "ruby")
		.mecdfcols ("black", c (0.4, 0, 0), c (1, 0.2, 0.2) )
	else if (theme == "silver")
		.mecdfcols ("grey30", c (0.60, 0.60, 0.62), c (0.95, 0.95, 0.94) )
	else if (theme == "gold")
		.mecdfcols (rgb (0.98, 1, 0.15), c (0.5, 0.55, 0), c (0.98, 1, 0.15) )
	else stop ("unsupported theme")
}

.mecdfcols = function (sline="black", sfill1=c (0.92, 0.92, 0.92), sfill2=sfill1,
	fline="grey70", ffill="grey95")
{	options (mecdf.frame=list (line=fline, fill=ffill) )
	options (mecdf.surface = list (line=sline, fill1=sfill1, fill2=sfill2) )
}



#' univariate continuous ecdf inverse
#' 
#' This function computes the inverse of a univariate continuous ECDF, with no
#' model correction. This function is deprecated. Unlike the mecdf function, it
#' doesn't return a function, rather it returns a scalar value.
#' 
#' 
#' @param x A numeric vector.
#' @param p Probability.
#' @param sort Whether or not to sort the data. Unless the data is already
#' sorted, this must be true.
#' @export ecdfinv
ecdfinv = function (x, p, sort=TRUE)
{	n = length (x)
	if (sort) x = sort (x)
	k = (n - 1) * p + 1
	if (k < 1) x [1]
	else if (k > n) x [n]
	else
	{	kf = floor (k)
		while (kf > 1 && x [kf] == x [kf + 1]) kf = kf - 1
		if (k == kf)
			x [kf]
		else
		{	kc = ceiling (k)
			while (kc < n && x [kc] == x [kc + 1]) kc = kc + 1
			a = (k - kf) / (kc - kf)
			b = 1 - a
			b * x [kf] + a * x [kc]
		}
	}
}


