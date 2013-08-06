
writeROI <-
function(roi, file.name)
{
  writeRoiLine <- function(t1, fpname)
  {
    outStr="x"
    colNames = names(t1)
    colNames = colNames[!colNames == "counts"]
    outStr = paste(  t1[colNames],collapse="\t")
    outStr = paste(c(outStr, paste(unlist(t1$counts),collapse=",")), collapse="\t")
    write(outStr, file=fpname, append=TRUE)
  }

   fpConn <- file(file.name, "w")
   colCount=1
   write("# ROI written from R", fpConn, append=F)
   colNames = names(roi)
   colNames = colNames[!colNames == "counts"]

   outStr=""
   outStr=paste(c(outStr, "##",colNames[1]), collapse="")
   for(colIds in c(2:length(colNames)))
   {
      outStr = paste(c(outStr, colNames[colIds]), collapse="\t")
   }
   outStr = paste(c(outStr, "counts"), collapse="\t")
   write(outStr, file=fpConn, append=TRUE)

   apply(roi,1,writeRoiLine,fpConn)
   close(fpConn)
}
