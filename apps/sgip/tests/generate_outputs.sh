#!/bin/sh
#
# Output generation script for sgip

SGIP=../../../../build/Debug/apps/sgip/sgip
DIRR01=../example/r01
DIRSRG=../example/srg
DIRISO=../example/Iso_Data/

${SGIP} -o ${DIRR01}/iso_r01_m200.A00 -c ${DIRR01}/iso_r01_m200.B00 -v 2 -i > iso_r01_m200.A00_B00.stdout
${SGIP} -o ${DIRR01}/iso_r01_m200.A01 -c ${DIRR01}/iso_r01_m200.B01 -v 2 -i > iso_r01_m200.A01_B01.stdout
${SGIP} -o ${DIRR01}/iso_r01_m200.A00 -c ${DIRR01}/iso_r01_m200.B01 -v 2 -i > iso_r01_m200.A00_B01.stdout
${SGIP} -o ${DIRR01}/iso_r01_m200.A00 -v 2 > iso_r01_m200.A00.stdout
${SGIP} -o ${DIRISO}/iso_m2D_m196.A00 -c ${DIRISO}/iso_m2D_m196.B00 -v 2 -i > iso_m2D_m196.A00_B00.stdout
${SGIP} -o ${DIRSRG}/latin-4 -v 2 > srg_latin-4.stdout
${SGIP} -o ${DIRSRG}/lattice-4 -v 2 > srg_lattice-4.stdout
${SGIP} -o ${DIRSRG}/paley-5 -v 2 > srg_paley-5.stdout
${SGIP} -o ${DIRSRG}/sts-7 -v 2 > srg_sts7.stdout
${SGIP} -o ${DIRSRG}/triang-5 -v 2 > srg_triang-5.stdout
