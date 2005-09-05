\name{matchControls}
\alias{matchControls}
\title{Find Matched Control Group}
\usage{
matchControls(formula, data = list(), subset, contlabel = "con",
               caselabel = NULL, dogrep = TRUE, replace = FALSE)
}
\arguments{
    \item{formula}{A formula indicating cases, controls and the
	variables to be matched. Details are described below.}
    \item{data}{an optional data frame containing the variables in the
	model.  By default the variables are taken from the environment
	which \code{matchControls} is called from.}
    \item{subset}{an optional vector specifying a subset of observations
 	to be used in the matching process.}
    \item{contlabel}{A string giving the label of the control group.}
    \item{caselabel}{A string giving the labels of the cases.}
    \item{dogrep}{If \code{TRUE}, then \code{contlabel} and
	\code{contlabel} are matched using \code{\link{grep}}, else
	string comparison (exact equality) is used.}
    \item{replace}{If \code{FALSE}, then every control is used only
	once.}
}
\description{
    Finds controls matching the cases as good as possible.
}
\details{
    The left hand side of the \code{formula} must be a factor
    determining whether an observation belongs to the case or the
    control group.  By default, all observations where a grep of
    \code{contlabel} matches, are used as possible controls, the rest is
    taken as cases.  If \code{caselabel} is given, then only those
    observations are taken as cases.  If \code{dogrep = TRUE}, then both
    \code{contlabel} and \code{caselabel} can be regular expressions.

    The right hand side of the \code{formula} gives the variables that
    should be matched.  The matching is done using the
    \code{\link{daisy}} distance from the \code{cluster} package, i.e.,
    a model frame is built from the formula and used as input for
    \code{\link{daisy}}. For each case, the nearest control is
    selected. If \code{replace = FALSE}, each control is used only
    once.
}
\value{
    Returns a list with components
    \item{cases}{Row names of cases.}
    \item{controls}{Row names of matched controls.}
    \item{factor}{A factor with 2 levels indicating cases and controls
	(the rest is set to \code{NA}.}
}
\author{Friedrich Leisch}
\examples{
Age.case <- 40 + 5 * rnorm(50)
Age.cont <- 45 + 10 * rnorm(150)
Age <- c(Age.case, Age.cont)

Sex.case <- sample(c("M", "F"), 50, prob = c(.4, .6), replace = TRUE)
Sex.cont <- sample(c("M", "F"), 150, prob = c(.6, .4), replace = TRUE)
Sex <- as.factor(c(Sex.case, Sex.cont))

casecont <- as.factor(c(rep("case", 50), rep("cont", 150)))

## now look at the group properties:
boxplot(Age ~ casecont)
barplot(table(Sex, casecont), beside = TRUE)

m <- matchControls(casecont ~ Sex + Age)

## properties of the new groups:
boxplot(Age ~ m$factor)
barplot(table(Sex, m$factor))
}
\keyword{manip}
