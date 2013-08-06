\name{ngsroi-package}
\alias{ngsroi-package}
\alias{ngsroi}
\docType{package}
\title{
NGS Regions of Interest Analysis
}
\description{
Routines for I/O of NGS ROI files and manipulation thereof.
}
\details{
\tabular{ll}{
Package: \tab ngsroi\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-04-12\cr
License: \tab BSD_3_clause\cr
}
}
\author{
Bernd Jagla <bernd.jagla@pasteur.fr>\cr
Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>

Maintainer: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
}
\references{
Jagla, B, Holtgrewe, M, Reinert, K: NGS ROI. To appear.
}
\keyword{ package }
\examples{
library(ngsroi)

# Load ROI file into data.frame.
#roi = readROI("dmel.bowtie.sam.roi")

# Compute some metrics.
#roi$min = as.numeric(lapply(roi$counts ,min))
#roi$median =  as.numeric(lapply(roi$counts ,median))
#roi$mean =  as.numeric(lapply(roi$counts ,mean))
#roi$quantile75 =  as.numeric(lapply(roi$counts ,quantile, probs=0.75))
#roi$quantile95 =  as.numeric(lapply(roi$counts ,quantile, probs=0.95))

# Write data.frame into ROI file again.
#writeROI(roi, "dmel.bowtie.sam.trans.roi");
}
