% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/printrcure.R
\name{printrcure}
\alias{printrcure}
\title{Prints rcure object}
\usage{
printrcure(x, Var = TRUE, ROC = TRUE, ...)
}
\arguments{
\item{x}{an object of rcure}

\item{Var}{if it is TRUE, the funciton returns standard error estimated by bootstrap method. If it is False, the function only returns estimators of coefficients. By default, Var = TRUE}

\item{ROC}{if TRUE, the function plots ROC. By default, ROC = TRUE}

\item{...}{further arguments to be passed to the printrcure function}
}
\description{
Prints rcure object
}
\references{
Cai, C., Zou, Y., Peng, Y., & Zhang, J. (2012). smcure: An R-Package for estimating semiparametric mixture cure models. Computer methods and programs in biomedicine, 108(3), 1255-1260
}

