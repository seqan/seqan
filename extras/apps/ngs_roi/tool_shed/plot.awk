###
#
# generate gnuplot input for generating ROI plots (9 per page)
# USAGE:   gawk -v fileName=${IN_FILE}.tmp.ps \
#               -f plot.awk \
#               ${USED_INFILE} |gnuplot  2> /dev/null
#  ${IN_FILE}.tmp.ps is the file gnuplot write to 
#  ${USED_INFILE} is the uncompressed ROI file.
#


BEGIN{
#   if(fileName==""){
#       print "fileName not specified: use -v fileName=..." 
#       exit
#   }
   print "set terminal postscript color blacktext solid \"Helvetica\" 6" 
   print "set output \""fileName"\"" 
   print "set multiplot" 
   print "set size 0.9,0.9" 
   print "set origin 0.1,0.1" 
   pos[1]="0.,0.6" 
   pos[2]="0.3,0.6" 
   pos[3]="0.6,0.6" 
   pos[4]="0.,0.3" 
   pos[5]="0.3,0.3" 
   pos[6]="0.6,0.3" 
   pos[7]="0.,0." 
   pos[8]="0.3,0." 
   pos[9]="0.6,0." 
   idx=0
}

!/^#/{
   #if($5>100){
   if(1==1){
   idx++
   print "set size 0.3,0.3" 
   print "set origin " pos[idx]
   print "set title \"" $1,$2,$3,$4 "\\n" $5, $6, $7"\"" 
   print "plot \"-\"  with lines" 
   split($NF, counts, ",")
   for(i=1; i<=length(counts);i++){
      print counts[i]
   }
   print "e" 
   if(idx==9){
      print "unset multiplot" 
      print "set size 0.9,0.9" 
      print "set origin 0.1,0.1" 
      print "set multiplot" 
      idx=0
   }
   }
}
