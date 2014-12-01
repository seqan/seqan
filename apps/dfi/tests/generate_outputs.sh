#!/bin/bash
#
# Currently, we run DFI only on

# We use the current trunk version of 2011-10-13 (r10463) for building the
# reference.
DFI=../../../../build/Release/apps/dfi/dfi
DFI=../../../../build/make.linux/apps/dfi/dfi

set1=( 25.H_sapiens.fa  human_ucsc.fa  bbe.fa CompWindows.fa )
set2=( 59.M_musculus.fa dros_promos.fa kjv.fa CompNonWindows.fa )
out=( Human_Mouse_Proteome human_dros_promoter BasicEnglish_KingJames_Bible Windows_Other )
flags=( -p -d -m -m )
rm -f params.txt results.zip *.res

for (( i = 0 ; i < ${#out[@]} ; i++ ))
do
	inputPos=${set1[$i]}
	inputNeg=${set2[$i]}
	result=${out[$i]}
	flag=${flags[$i]}
	posis=$(grep -c ">" $inputPos)
	negis=$(grep -c ">" $inputNeg)
	newmin=$(echo "scale=0; $posis*0.005/1"| bc)
	newmax=$(echo "scale=0; $negis*0.5/1" | bc)

	# Mine emerging substrings
	RESFILE=${result}_emerging_0.002_5.res
	${DFI} ${inputPos} ${inputNeg} --support 0.002 --growth 5 ${flag} > ${RESFILE}
	echo ${RESFILE} ${inputPos} ${inputNeg} --support 0.002 --growth 5 ${flag} >> params.txt

	# Mine minmax substrings
	RESFILE=${result}_minmax_${newmin}_${newmax}.res
	${DFI} ${inputPos} ${inputNeg} --minmax ${newmin} ${posis} --minmax 1 ${newmax} ${flag} > ${RESFILE}
	echo ${RESFILE} ${inputPos} ${inputNeg} --minmax ${newmin} ${posis} --minmax 1 ${newmax} ${flag} >> params.txt
done

inputPos=CompWindows.fa
inputNeg=CompNonWindows.fa
result=Windows_Other
flag=-m

# Mine emerging substrings with different minimal supports
growth=5
for supp in 0.0005 0.001 0.005 0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1
do
	RESFILE=${result}_emerging_${supp}_${growth}.res
	${DFI} ${inputPos} ${inputNeg} --support ${supp} --growth ${growth} ${flag} > ${RESFILE}
	echo ${RESFILE} ${inputPos} ${inputNeg} --support ${supp} --growth ${growth} ${flag} >> params.txt
done

# Mine emerging substrings with different minimal growth rates
supp=0.002
for growth in 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10
do
	RESFILE=${result}_emerging_${supp}_${growth}.res
	${DFI} ${inputPos} ${inputNeg} --support ${supp} --growth ${growth} ${flag} > ${RESFILE}
	echo ${RESFILE} ${inputPos} ${inputNeg} --support ${supp} --growth ${growth} ${flag} >> params.txt
done

# run dummy example (see dfi README)
supp=1
growth=2
inputPos=fasta1.fa
inputNeg=fasta2.fa
result=Example
flag=

	RESFILE=${result}_emerging_${supp}_${growth}.res
    ${DFI} ${inputPos} ${inputNeg} --support ${supp} --growth ${growth} ${flag} > ${RESFILE}
    echo ${RESFILE} ${inputPos} ${inputNeg} --support ${supp} --growth ${growth} ${flag} >> params.txt


zip results.zip *.res
