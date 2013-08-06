\name{writeROI}
\alias{writeROI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write ROI file.
}
\description{
Write data.frame to a ROI file.
}
\usage{
writeROI(roi, file.name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{roi}{
    The data.frame to write to file.
    }
  \item{file.name}{
    The path to the file to write to.
    } 
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
%%\value{  
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}
%%\references{
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (roi, fpname) 
{
    fpConn <- file(fpname, "w")
    colCount = 1
    write("# ROI written from R", fpConn, append = F)
    colNames = names(roi)
    colNames = colNames[!colNames == "counts"]
    outStr = ""
    outStr = paste(c(outStr, "##", colNames[1]), collapse = "")
    for (colIds in c(2:length(colNames))) {
        outStr = paste(c(outStr, colNames[colIds]), collapse = "\t")
    }
    outStr = paste(c(outStr, "counts"), collapse = "\t")
    write(outStr, file = fpConn, append = TRUE)
    apply(roi, 1, writeRoiLine, fpConn)
    close(fpConn)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
