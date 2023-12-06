\name{hsv_palette}
\alias{hsv_palette}
\title{Sequential color palette based on HSV colors}
\description{
Computes a sequential color palette based on HSV colors by varying the
saturation, given hue and value.
}
\usage{
hsv_palette(h = 2/3, from = 0.7, to = 0.2, v = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{hue}
  \item{from}{lower bound for saturation}
  \item{to}{upper bound for saturation}
  \item{v}{value}
}
\value{
  A function with one argument: the size of the palette, i.e., the number
  of colors.
}
\author{David Meyer \email{David.Meyer@R-project.org}}
\seealso{\code{\link[grDevices]{hsv}}}
\examples{
pie(rep(1, 10), col = hsv_palette()(10))
pie(rep(1, 10), col = hsv_palette(h = 0)(10))
}
\keyword{hplot}% __ONLY ONE__ keyword per line
