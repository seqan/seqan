#!/bin/bash
# ----------------------------------------------------------------------------
#                  SeqAn - The Library for Sequence Analysis
# ----------------------------------------------------------------------------
# Copyright (c) 2006-2016, Knut Reinert, FU Berlin
# All rights reserved.
#
# License: BSD 3-clause
# ----------------------------------------------------------------------------
# Author: Sabrina Krakau <sabrina.krakau@fu-berlin.de>
# Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
# ----------------------------------------------------------------------------
# Binary test helper script.
#
# This script downloads binary packages of SeqAn apps on remote Linux machines
# and calls the programs to test whether the applications can be executed at
# all and thus check whether the system's libraries are compatible.
#
# The following packages have to be installed in order to run on the operating
# systems:
#
# SUSE 64 Bit:
#   zypper install libgomp1-32bit
#
# Fedora 64 Bit:
#   yum install -y glibc.i686 libstdc++.i686 zlib.i686 libgomp.i686
# ----------------------------------------------------------------------------

# Global variables with user name and target host.
USER_NAME=$USER
SERVER_NAMES=

# Parse option values.
while getopts ":u:h:" opt
do
    case $opt in
        u)
            USER_NAME=$OPTARG
            ;;
        h)
            SERVER_NAMES=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done


# Check that ${USER_NAME} and ${SERVER_NAMES} are set.
if [[ "$USER_NAME" == "" || "$SERVER_NAMES" == "" ]]; then
    echo "Missing option -h" >&2
    exit 1
fi

# List of binary packages to tests.
BINARY_URLS=(
    "http://packages.seqan.de/razers3/razers3-3.1.1-Linux-i686.tar.bz2"
    "http://packages.seqan.de/razers3/razers3-3.1.1-Linux-x86_64.tar.bz2"
    "http://packages.seqan.de/stellar/stellar-1.4.1-Linux-i686.tar.bz2"
    "http://packages.seqan.de/stellar/stellar-1.4.1-Linux-x86_64.tar.bz2"
    "http://packages.seqan.de/breakpoint_calculator/breakpoint_calculator-0.2.1-Linux-i686.tar.bz2"
    "http://packages.seqan.de/breakpoint_calculator/breakpoint_calculator-0.2.1-Linux-x86_64.tar.bz2"
    "http://packages.seqan.de/snp_store/snp_store-1.0.1-Linux-i686.tar.bz2"
    "http://packages.seqan.de/snp_store/snp_store-1.0.1-Linux-x86_64.tar.bz2")

# Array with test calls for each package.
APP_CALLS=(
    "./razers3 -o output.razers ../example/genome.fa ../example/reads.fa"
    "./razers3 -o output.razers ../example/genome.fa ../example/reads.fa"
    "./stellar -o ../example/mapped_reads.gff --minLength 30 ../example/NC_001477.fasta ../example/reads.fasta" 
    "./stellar -o ../example/mapped_reads.gff --minLength 30 ../example/NC_001477.fasta ../example/reads.fasta"
    "./breakpoint_calculator -d2 ../example/alignment.maf"
    "./breakpoint_calculator -d2 ../example/alignment.maf"
    "./snp_store ../example/exampleGenome.fa ../example/exampleReads.gff -o output_snp.txt -id output_indel1.txt"
    "./snp_store ../example/exampleGenome.fa ../example/exampleReads.gff -o output_snp.txt -id output_indel1.txt")

# Create tmp directory, download binaries, call program and return success/error message
#
# Arguments:
#
#  $1 server to connect to
#  $2 call to application
function test_app ()
{
    USER_SERVER=$1
    APP_CALL=$2
    # Create temporary directory.
    echo "testing ${F}" >&2
   	TMPDIR=$(ssh $USER_SERVER "mktemp -d")

    # Download package.
    echo -n "  downloading ..." >&2
	ssh $USER_SERVER "cd $TMPDIR && curl $F -o binaryName.tar.bz2 &>/dev/null"
    if [[ "$?" = "0" ]]; then
        echo " OK" >&2
    else
        echo " FAILED" >&2
        echo "FAILED"  # for caller
        return
    fi

    # Extracting package.
    echo -n "  unpacking ..." >&2
	ssh $USER_SERVER "cd $TMPDIR && tar -xjf binaryName.tar.bz2"
    if [[ "$?" = "0" ]]; then
        echo " OK" >&2
    else
        echo " FAILED" >&2
        echo "FAILED"  # for caller
        return
    fi

    # Calling binary.
    echo -n "  calling test ..." >&2
    FILE=$(mktemp)
    MSG=`ssh $USER_SERVER "cd $TMPDIR && cd */bin && (${APP_CALL}) > file.stdout" &>${FILE}`
    RET=$?
    if [[ "${RET}" = "0" ]]; then
        echo " OK" >&2
        echo "OK"  # for caller
        return
    else
        echo " FAILED" >&2
        echo " The program call returns: ${RET}." >&2
        echo " Program output is >>>" >&2
        cat ${FILE} >&2
        rm -f ${FILE}
        echo ${MSG}
        echo "<<<" >&2
        echo "FAILED"  # for caller
        return
    fi

    # Remove temporary directory.
	ssh $USER_SERVER "rm -r $TMPDIR"
    
}

echo "Testing SeqAn Linux Binaries..." >&2
echo >&2

RESULTS=()
IFS=";"; 

# Execute tests for each server.
for SERVER_NAME in ${SERVER_NAMES[@]}; do
    echo "testing on ${SERVER_NAME}"
	USER_SERVER="$USER_NAME@$SERVER_NAME"
	COUNTER=0
	for F in ${BINARY_URLS[@]}; do
        #echo "${APP_CALLS[$COUNTER]}" >&2
        RET=`test_app $USER_SERVER "${APP_CALLS[$COUNTER]}"`
	    RESULTS+=("${RET}")
		let COUNTER=COUNTER+1
	done
done

# For each server and for each binary-url: Output test result
echo -e '\t'
for (( I = 0 ; I < ${#SERVER_NAMES[@]} ; I++ )) do
    echo "Binary tests on ${SERVER_NAMES[$I]}:"
    echo ""
    for (( J = 0 ; J < ${#BINARY_URLS[@]} ; J++ )) do
        K=$(($J+$I*${#BINARY_URLS[@]}))
        echo -e "${BINARY_URLS[$J]}\t${RESULTS[$K]}"
    done
done


