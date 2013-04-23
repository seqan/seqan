#!/bin/bash

SEQAN_ROOT=`basename $0/..`
SNAPSHOT_NAME=seqan-1.3.1

rm -rf ${SNAPSHOT_NAME}
rm -f ${SNAPSHOT_NAME}.zip ${SNAPSHOT_NAME}.tar.gz

pushd ${SEQAN_ROOT}/docs
./make.sh
popd

mkdir ${SNAPSHOT_NAME}
cp -r ${SEQAN_ROOT}/demos ${SEQAN_ROOT}/seqan ${SEQAN_ROOT}/lib ${SEQAN_ROOT}/apps ${SEQAN_ROOT}/cmake ${SEQAN_ROOT}/CHANGELOG ${SEQAN_ROOT}/README ${SNAPSHOT_NAME}
cp -r ${SEQAN_ROOT}/docs/html ${SNAPSHOT_NAME}/docs
find ${SNAPSHOT_NAME} -name '.svn' | xargs rm -rf
find ${SNAPSHOT_NAME} -name '._\*' | xargs rm -rf
${SEQAN_ROOT}/misc/build_forwards.py ${SNAPSHOT_NAME}/seqan
tar czf ${SNAPSHOT_NAME}.tar.gz ${SNAPSHOT_NAME}
zip -qr ${SNAPSHOT_NAME}.zip ${SNAPSHOT_NAME}
