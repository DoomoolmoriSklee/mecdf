#' surface plots of bivariate cdfs
#' 
#' Prototype functions for plotting bivariate CDFs using three dimensional
#' graphics. Note that setting regular to false, can be very slow. Also note
#' bcdfplot and bcdf.plot are deprecated.
#' 
#' 
#' @aliases plotbcdf plotbcdf.mecdf plotbcdf.matrix bcdfplot bcdf.plot
#' @param m Either an mecdf object or a square matrix.
#' @param regular Logical (defaults to true), whether or not to use a regularly
#' spaced grid. If false, the plot uses the actual data values.
#' @param res Ignored unless regular is true. The number of points in both
#' directions of the grid.
#' @param ulim Ignored unless regular is true. A vector giving the minimum and
#' maximum values for the first variable on the grid. Defaults to the minimum
#' and maximum observed values.
#' @param vlim Same as ulim, except for the second variable.
#' @param mmin Minimum value.
#' @param mmax Maximum value.
#' @param labs .
#' @param \dots .
#' @export plotbcdf
plotbcdf = function (m, ...) UseMethod ("plotbcdf")

plotbcdf.mecdf = function (m, regular=TRUE, res=16, ulim, vlim, ...)
{	x = m$x
	if (missing (ulim) ) ulim = range (x [,1])
	if (missing (vlim) ) vlim = range (x [,2])
	if (regular)
	{	u = seq (ulim [1], ulim [2], length=res)
		v = seq (vlim [1], vlim [2], length=res)
		mst = matrix (numeric (), nr=res, nc=res)
		for (i in 1:res) for (j in 1:res) mst [i, j] = m (c (u [i], v [j]) )
		plotbcdf (mst, ..., labs=colnames (x) )
	}
	else
	{	u = sort (unique (x [,1]) )
		v = sort (unique (x [,2]) )
		if (!m$continuous)
		{	uo = 0.1 * diff (range (u) )
			vo = 0.1 * diff (range (v) )
			u = c (u [1] - uo, u, u [length (u)] + uo)
			v = c (v [1] - vo, v, v [length (v)] + vo)
		}
		nu = length (u)
		nv = length (v)
		mst = matrix (numeric (), nr=nu, nc=nv)
		for (i in 1:nu) for (j in 1:nv) mst [i, j] = m (c (u [i], v [j]) )
		urng = range (u)
		vrng = range (v)
		u = (u - urng [1]) / diff (urng)
		v = (v - vrng [1]) / diff (vrng)
		plotbcdf.matrix (NULL, ..., labs=colnames (x) )
		if (m$continuous) .plotbcdf.irregulargrid (mst, u, v)
		else .plotbcdf.bsf (mst, u, v)
	}
}

plotbcdf.matrix = function (m, mmin=0, mmax=1, ..., labs=NULL)
{	p0 = par (mar=c (1.75, 0.175, 0.9, 0.175) )
	plot.new ()
	plot.window (c (-0.75, 0.75), c (0, 1.5) )
	.plotbcdf.plane.uv ()
	.plotbcdf.plane.u2 ()
	.plotbcdf.plane.v2 ()
	if (ifo (labs) )
		.plotbcdf.labs (rev (labs) )
	if (ifo (m) )
	{	n = nrow (m)
		if (n != ncol (m) ) stop ("square matrix required")
		m = (m - mmin) / (mmax - mmin)
		.plotbcdf.regulargrid (m, n)
	}
	par (p0)
}

.plotbcdf.regulargrid = function (m, n)
{	n2 = n - 1
	incr = 1 / n2
	for (i in 1:n2) for (j in 1:n2)
	{	u1 = (i - 1) * incr
		u2 = u1 + incr
		v1 = (j - 1) * incr
		v2 = v1 + incr
		u = c (u1, u1, u2, u2)
		v = c (v1, v2, v2, v1)
		w = c (m [i, j], m [i, j + 1], m [i + 1, j + 1], m [i + 1, j])
		dir = (w [3] - w [1]) / sqrt (2 * incr^2)
		if (dir < 0) dir = 0
		#scaling to interval (0, 1), tangent=1 -> dir=0.5
		dir = 2 * atan (dir) / pi
		.plotbcdf.poly (u, v, w, getOption ("mecdf.surface")$line,
			.colinterp (dir) )
	}
}

.plotbcdf.irregulargrid = function (m, u, v)
{	nu = length (u) - 1
	nv = length (v) - 1
	for (i in 1:nu) for (j in 1:nv)
	{	u1 = u [i]
		u2 = u [i + 1]
		v1 = v [j]
		v2 = v [j + 1]
		up = c (u1, u1, u2, u2)
		vp = c (v1, v2, v2, v1)
		w = c (m [i, j], m [i, j + 1], m [i + 1, j + 1], m [i + 1, j])
		dir = (w [3] - w [1]) / sqrt ( (u2 - u1)^2 + (v2 - v1)^2)
		if (dir < 0) dir = 0
		dir = 2 * atan (dir) / pi
		.plotbcdf.poly (up, vp, w, getOption ("mecdf.surface")$line,
			col=.colinterp (dir) )
	}
}

.plotbcdf.bsf = function (m, u, v)
{	nu = length (u) - 1
	nv = length (v) - 1
	bc = rgb (0.08, 0.5, 0.5)
	fc = rgb (0, 0.8, 0.1)
	ec = rgb (0.7, 0.8, 0.9)
	for (i in nu:1) for (j in nv:1)
	{	u1 = u [i]
		u2 = u [i + 1]
		v1 = v [j]
		v2 = v [j + 1]
		up = c (u1, u1, u2, u2)
		vp = c (v1, v2, v2, v1)
		w0 = m [i, j]
		w = c (w0, w0, w0, w0)
		up1 = c (u1, u1, u1, u1)
		vp1 = c (v1, v1, v2, v2)
		w1 = c (0, w0, w0, 0)
		up2 = c (u1, u1, u2, u2)
		vp2 = c (v1, v1, v1, v1)
		dir = (w [3] - w [1]) / sqrt ( (u2 - u1)^2 + (v2 - v1)^2)
		if (i == 1 || i == nu || j == 1 || j == nv)
		{	.plotbcdf.poly (up1, vp1, w1, border=bc, col=ec)
			.plotbcdf.poly (up2, vp2, w1, border=bc, col=ec)
			.plotbcdf.poly (up, vp, w, border=bc, col=ec)
		}
		else
		{	.plotbcdf.poly (up1, vp1, w1, border=bc, col=fc)
			.plotbcdf.poly (up2, vp2, w1, border=bc, col=fc)
			.plotbcdf.poly (up, vp, w, border=bc, col="darkgreen")
		}
	}
}

.plotbcdf.plane.uv = function ()
{	.plotbcdf.poly (c (0, 0, 1, 1), c (0, 1, 1, 0), 0)
	p1 = .project (0, 0, 0)
	p2 = .project (0, 1, 0)
	p3 = .project (1, 0, 0)
	z = 0.03
	col = getOption ("mecdf.frame")$line
	arrows (p1 [1] - z, p1 [2] - z, p2 [1] - z, p2 [2] - z, col=col, length=0.12)
	arrows (p1 [1] + z, p1 [2] - z, p3 [1] + z, p3 [2] - z, col=col, length=0.12)
}

.plotbcdf.plane.u2 = function ()
{	.plotbcdf.poly (1, c (0, 0, 1, 1), c (0, 1, 1, 0) )
	.plotbcdf.lines (1, 0:1, 0.25)
	.plotbcdf.lines (1, 0:1, 0.5)
	.plotbcdf.lines (1, 0:1, 0.75)
	.plotbcdf.lines (1, 0.25, 0:1)
	.plotbcdf.lines (1, 0.5, 0:1)
	.plotbcdf.lines (1, 0.75, 0:1)
}

.plotbcdf.plane.v2 = function ()
{	.plotbcdf.poly (c (0, 0, 1, 1), 1, c (0, 1, 1, 0) )
	.plotbcdf.lines (0:1, 1,  0.25)
	.plotbcdf.lines (0:1, 1, 0.5)
	.plotbcdf.lines (0:1, 1, 0.75)
	.plotbcdf.lines (0.25, 1,  0:1)
	.plotbcdf.lines (0.5, 1, 0:1)
	.plotbcdf.lines (0.75, 1, 0:1)
}

.plotbcdf.labs = function (labs)
{	x = c (-0.525, 0.525)
	y = 0.16
	text (x, y, labs)
}

.plotbcdf.poly = function (u, v, w, border=getOption ("mecdf.frame")$line,
	col=getOption ("mecdf.frame")$fill)
{	m = .project (u, v, w)
	polygon (m [,1], m [,2], border=border, col=col)
}

.plotbcdf.lines = function (u, v, w, col=getOption ("mecdf.frame")$line)
{	m = .project (u, v, w)
	lines (m [,1], m [,2], col=col)
}

#todo, fix this...
.project = function (u, v, w)
{	x = u * cos (pi / 4) + v * cos (pi * 3 / 4)
	y = u * sin (pi / 4) + v * sin (pi * 3 / 4)
	y = 0.71 * y + 0.5 * w
	cbind (x, y)
}

.colinterp = function (x)
{	col1 = getOption ("mecdf.surface")$fill1
	col2 = getOption ("mecdf.surface")$fill2
	col = col1 + x * (col2 - col1)
	rgb (col [1], col [2], col [3])
}




