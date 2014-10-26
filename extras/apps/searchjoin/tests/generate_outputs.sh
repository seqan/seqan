#!/opt/local/bin/bash
# REQUIRES bash version 4.x
#
# Search/Join reference output generation.

# Paths to binaries
BIN=~/Documents/Code/SeqAn/build/Release/bin
SEARCH=$BIN/s4_search
JOIN=$BIN/s4_join

# Threads
THREADS="4"

# Errors per dataset
declare -A K
K[geo]="0 1 2 3"
K[dna]="0 4 8 12 16"

# Seed lengths per dataset
SL[geo]="4 5 6"
SL[dna]="11 13 15"

# ============================================================
# Run Join App
# ============================================================

# Run online join
for alphabet in geo dna; do
    for k in ${K[${alphabet}]}; do
        ${JOIN} -i ${alphabet} -t ${THREADS} -o join_${alphabet}\_${k}\_online.out ${alphabet}\_database.csv ${k} --online # &> join_${alphabet}\_${k}\_online.stdout
        sort --unique join_${alphabet}\_${k}\_online.out > join_${alphabet}\_${k}.out
        rm join_${alphabet}\_${k}\_online.out
    done
done

# Run offline join
#for alphabet in geo dna; do
#    for k in ${K[${alphabet}]}; do
#        for sl in ${SL[${alphabet}]}; do
#            ${JOIN} -i ${alphabet} -t ${THREAds} -o join_${alphabet}\_${k}\_${sl}.out ${alphabet}\_database.csv ${k} --seed-length ${sl} # &> join_${alphabet}\_${k}\_${sl}.stdout
#        done
#    done
#done

# ============================================================
# Run Search App
# ============================================================

# Run online search
for alphabet in geo dna; do
    ${SEARCH} --no-wait -i ${alphabet} -t ${THREADS} -o search_${alphabet}\_online.out ${alphabet}\_database.csv ${alphabet}\_queries.csv --online # &> search_${alphabet}\_online.stdout
    sort --unique search_${alphabet}\_online.out > search_${alphabet}.out
    rm search_${alphabet}\_online.out
done

# Run offline search
#for alphabet in geo dna; do
#    for sl in ${SL[${alphabet}]}; do
#        ${SEARCH} --no-wait -i ${alphabet} -t ${THREADS} -o search_${alphabet}\_${sl}.out ${alphabet}\_database.csv ${alphabet}\_queries.csv --seed-length ${sl} # &> search_${alphabet}\_${sl}.stdout
#    done
#done