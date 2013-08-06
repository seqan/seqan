readROI <-
function(file.name){
   # return element from string array
   fi <- function(x,i){
      x[i]
   }
   # return element from string array as number
   fni <- function(x,i){
      as.numeric(x[[i]])
   }
   # return vector of counts for last column
   getVec <- function(x){
      mylen=as.numeric(x[5])
      veclen=length(x)
      as.numeric(unlist(strsplit(x[veclen],",")))
   }
   #get non-fixed values (i.e. all columns between the cg_count and the count vector
   getVals <- function(y, x, columnNames){
      #mylen=as.numeric(x[[1]][5])
      veclen=length(columnNames)
      t1=veclen - 1
      for ( colN in c(7:t1)) {
   	y[,columnNames[colN]]=unlist(lapply(x,fni,colN))
      }
      return(y)
   }

   con=gzfile(file.name)
   rLines = readLines(con)
   close(con)
   values = rLines[substr(rLines,1,1)!="#"]
   values = strsplit(values,"\t")
   if (length(rLines[substr(rLines,1,2)=="##"])==0){
      columnNames=unlist(list("##ref", 
                              "begin_pos",
			      "end_pos",
			      "region_name",
			      "length", 
			      "strand",
			      "max_count",
			      "gc_content",
			      "counts"))
   }else{
      columnNames = unlist(strsplit(rLines[substr(rLines,1,2)=="##"],"\t"))
   }



   df=data.frame(ref=unlist(lapply(values,fi,1)),
                 begin_pos=as.integer(unlist(lapply(values,fni,2))),
                 end_pos=as.integer(unlist(lapply(values,fni,3))),
                 region_name=unlist(lapply(values,fi,4)),
                 length=as.integer(unlist(lapply(values,fni,5))),
                 strand=unlist(lapply(values,fi,6))
		 )

   df=getVals(df, values, columnNames)
   df$counts = lapply(values,getVec)
   df$counts=lapply(df$counts,unlist)
   roiNames=names(df)
   return (df)
}
