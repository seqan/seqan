#!/bin/sh

RAZERS=../razers

CMD="${RAZERS} 01_seq.fasta 01_reads.fasta -id -rr 100 -i 100 -of 4 -o 01_golden.sam"
echo $CMD
${CMD}

CMD="${RAZERS} 02_seq.fasta 02_reads.fasta -id -rr 100 -i 95 -of 4 -o 02_golden.sam"
echo $CMD
${CMD}

CMD="${RAZERS} 03_seq.fasta 03_reads.fasta -id -rr 100 -i 100 -of 4 -o 03_golden.sam"
echo $CMD
${CMD}

CMD="${RAZERS} 04_seq.fasta 04_reads.fasta -id -rr 100 -i 95 -of 4 -o 04_golden.sam"
echo $CMD
${CMD}

CMD="${RAZERS} 05_seq.fasta 05_reads.fasta -id -rr 100 -i 100 -of 4 -o 05_golden.sam"
echo $CMD
${CMD}

CMD="${RAZERS} 06_seq.fasta 06_reads.fasta -id -rr 100 -i 91.9849241194 -of 4 -o 06_golden.sam"
echo $CMD
${CMD}

CMD="${RAZERS} 07_seq.fasta 07_reads.fasta -id -rr 100 -i 95.9924620597 -of 4 -o 07_golden.sam"
echo $CMD
${CMD}

CMD="${RAZERS} 09_seq.fasta 09_reads.fasta -id -i 100 -of 4 -o 09_golden.sam"
echo $CMD
${CMD}

CMD="${RAZERS} 10_seq.fasta 10_reads.fasta -id -i 95 -of 4 -o 10_golden.sam"
echo $CMD
${CMD}

CMD="${RAZERS} 11_seq.fasta 11_reads.fasta -id -i 100 -of 4 -o 11_golden.sam"
echo $CMD
${CMD}

CMD="${RAZERS} 12_seq.fasta 12_reads.fasta -id -i 95 -of 4 -o 12_golden.sam"
echo $CMD
${CMD}
