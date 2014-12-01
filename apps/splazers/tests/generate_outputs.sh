#!/bin/sh
#
# Output generation script for splazers

SPLAZERS="../../../../build/Debug/apps/splazers/splazers" # " -pd ../gapped_params/"
#SPLAZERS="../../../../build/linuxRelease/apps/splazers/splazers -pd ../gapped_params/"

# ============================================================
# First Section
# ============================================================


for rl in 100; do
  # -pd ../gapped_paramsefault run  
  ${SPLAZERS} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1_default.out > se-adeno-reads${rl}_1_default.stdout
  for mml in 16 17 18 19 20 21 22 23 24 25 ; do
 # for mml in 23 ; do
    for ep in 0 1  ; do
      for es in 0 1 2 ; do
        # Run with defaults for everything.
	${SPLAZERS}  -sm $mml -ep $ep -es $es adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1_mml${mml}_ep${ep}_es${es}.out > se-adeno-reads${rl}_1_mml${mml}_ep${ep}_es${es}.stdout
    	# Allow indels.
#    	${SPLAZERS}  -sm $mml -ep $ep -es $es -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1_mml${mml}-id_ep${ep}_es${es}.out > se-adeno-reads${rl}_1_mml${mml}-id_ep${ep}_es${es}.stdout

     done
    done
  done

  # Compute forward/reverse maches only.
  ${SPLAZERS} -id -sm 20 -ep 1 -es 1 -f adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1_mml20-f_ep1_es1.out > se-adeno-reads${rl}_1_mml20-f_ep1_es1.stdout
  ${SPLAZERS} -id -sm 20 -ep 1 -es 1 -r adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1_mml20-r_ep1_es1.out > se-adeno-reads${rl}_1_mml20-r_ep1_es1.stdout


  # Compute with different identity rates.
  for i in  90 91 92 93 94 95 96 97 98 99 100; do
      ${SPLAZERS} -sm 20 -ep 1 -es 1 -i ${i} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1_mml20-i${i}_ep1_es1.out > se-adeno-reads${rl}_1_mml20-i${i}_ep1_es1.stdout
    done

    # Run with different output formats.
    for of in 3 4 ; do
        ${SPLAZERS} -sm 20 -ep 1 -es 1 -of ${of}  adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1_mml20-of${of}_ep1_es1.out > se-adeno-reads${rl}_1_mml20-of${of}_ep1_es1.stdout
    done

    # Run with different match ordering.
    for so in 0 1; do
        ${SPLAZERS} -sm 20 -ep 1 -es 1 -so ${so} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1_mml20-so${so}_ep1_es1.out > se-adeno-reads${rl}_1_mml20-so${so}_ep1_es1.stdout
    done
	
    # run paired end anchored mode
#    ${SPLAZERS} -an adeno-genome.fa adeno-reads${rl}_1.sam -o pe-adeno-reads${rl}_2_default.out > pe-adeno-reads${rl}_2_default.stdout
done

