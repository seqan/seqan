\name{readROI}
\alias{readROI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read ROI file.
}
\description{
Read ROI file into data.frame.
}
\usage{
readROI(file.name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.name}{
    The path to the ROI file to read.
    }
}
\value{
  The function returns a data.frame with the data from the ROI file.
}
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
function (fpname) 
{
    fi <- function(x, i) {
        x[i]
    }
    fni <- function(x, i) {
        as.numeric(x[[i]])
    }
    getVec <- function(x) {
        mylen = as.numeric(x[5])
        veclen = length(x)
        as.numeric(unlist(strsplit(x[veclen], ",")))
    }
    getVals <- function(y, x, columnNames) {
        veclen = length(columnNames)
        t1 = veclen - 1
        for (colN in c(7:t1)) {
            y[, columnNames[colN]] = unlist(lapply(x, fni, colN))
        }
        return(y)
    }
    con = gzfile(fpname)
    rLines = readLines(con)
    close(con)
    values = rLines[substr(rLines, 1, 1) != "#"]
    values = strsplit(values, "\t")
    if (length(rLines[substr(rLines, 1, 2) == "##"]) == 0) {
        columnNames = unlist(list("##ref", "begin_pos", "end_pos", 
            "region_name", "length", "strand", "max_count", "cg_content", 
            "counts"))
    }
    else {
        columnNames = unlist(strsplit(rLines[substr(rLines, 1, 
            2) == "##"], "\t"))
    }
    df = data.frame(ref = unlist(lapply(values, fi, 1)), begin_pos = as.integer(unlist(lapply(values, 
        fni, 2))), end_pos = as.integer(unlist(lapply(values, 
        fni, 3))), region_name = unlist(lapply(values, fi, 4)), 
        length = as.integer(unlist(lapply(values, fni, 5))), 
        strand = unlist(lapply(values, fi, 6)))
    df = getVals(df, values, columnNames)
    df$counts = lapply(values, getVec)
    df$counts = lapply(df$counts, unlist)
    roiNames = names(df)
    return(df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
