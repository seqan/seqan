BEGIN{
longstring="\
[ /Rect [124 87 240 257]\n\
/Border [0 0 0]\n\
/Color [0 0 1]\n\
/Action << /Subtype /URI /URI (__URI1__) >>\n\
/Subtype /Link\n\
/ANN pdfmark\n\
[ /Rect [275 87 391 257]\n\
/Border [0 0 0]\n\
/Color [0 0 1]\n\
/Action << /Subtype /URI /URI (__URI4__) >>\n\
/Subtype /Link\n\
/ANN pdfmark\n\
[ /Rect [426 87 542 257]\n\
/Border [0 0 0]\n\
/Color [0 0 1]\n\
/Action << /Subtype /URI /URI (__URI7__) >>\n\
/Subtype /Link\n\
/ANN pdfmark\n\
[ /Rect [124 297 240 473]\n\
/Border [0 0 0]\n\
/Color [0 0 1]\n\
/Action << /Subtype /URI /URI (__URI2__) >>\n\
/Subtype /Link\n\
/ANN pdfmark\n\
[ /Rect [275 297 391 473]\n\
/Border [0 0 0]\n\
/Color [0 0 1]\n\
/Action << /Subtype /URI /URI (__URI5__) >>\n\
/Subtype /Link\n\
/ANN pdfmark\n\
[ /Rect [426 297 542 473]\n\
/Border [0 0 0]\n\
/Color [0 0 1]\n\
/Action << /Subtype /URI /URI (__URI8__) >>\n\
/Subtype /Link\n\
/ANN pdfmark\n\
[ /Rect [124 507 240 689]\n\
/Border [0 0 0]\n\
/Color [0 0 1]\n\
/Action << /Subtype /URI /URI (__URI3__) >>\n\
/Subtype /Link\n\
/ANN pdfmark\n\
[ /Rect [275 507 391 689]\n\
/Border [0 0 0]\n\
/Color [0 0 1]\n\
/Action << /Subtype /URI /URI (__URI6__) >>\n\
/Subtype /Link\n\
/ANN pdfmark\n\
[ /Rect [426 507 542 689]\n\
/Border [0 0 0]\n\
/Color [0 0 1]\n\
/Action << /Subtype /URI /URI (__URI9__) >>\n\
/Subtype /Link\n\
/ANN pdfmark\n\
showpage\n\
"
if(roiFile==""){
   print "Error: I need the variable roiFile to be set (-v roiFile=filename)" > "/dev/stderr"
   exit -1
}
}
/showpage/{
   outstring = longstring
   for(i=1;i<=9;i++){
      roiLine=""
      while(roiLine==""){
         if( (getline roiLine < roiFile) <= 0){
	    break
	 }
	 if(substr(roiLine,1,1)=="#"){
	    roiLine=""
	 }
      }
      split(roiLine,roiArr, "\t")
      gsub("__URI"i"__","http://localhost:60151/goto?locus="roiArr[1]":"roiArr[2]"-"roiArr[3],outstring)
   }
   print outstring
   next
}
{print}
