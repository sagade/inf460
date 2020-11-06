\name{Multiplot}
\alias{Multiplot}
\title{muliplot}
\usage{
  Multiplot(..., plotlist = NULL, cols)
}
\arguments{
  \item{...}{plots passed one by one}

  \item{plotlist}{plots passed in a list (appended to ...)}

  \item{cols}{number of cols in the plots table}
}
\value{
  none
}
\description{
  Grid of Plots in ggplot2
}
\details{
  Arrange multiple plots in a simple table layout in
  ggplot2.
}
\examples{
tmpp1 <- ggheatmap(matrix(rnorm(50*10), 50, 10))
tmpp2 <- ggheatmap(matrix(rnorm(50*10), 50, 10))
multiplot(tmpp1, tmpp2, cols=1)
}
\author{
  Winston Chang
  \url{http://wiki.stdout.org/rcookbook/Graphs/Multiple
  graphs on one page (ggplot2)/}
}

