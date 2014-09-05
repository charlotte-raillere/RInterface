\name{transmart.getGeneSignature}
\alias{transmart.getGeneSignature}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function returns Gene signature data
}
\description{
This function will retrieve a gene signature based on its name.
}
\usage{
transmart.getGeneSignature(name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
	The name of a gene signature.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame is returned with the Gene symbol, the probe identifier (can be NA) and the Fold change (can be NA)
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#Assumes a connection to the tranSMART DB has been made.

geneSig <- transmart.getGeneSignature(name = 'gene signature test')


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
