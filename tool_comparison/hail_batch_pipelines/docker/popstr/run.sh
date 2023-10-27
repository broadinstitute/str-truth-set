#!/bin/bash
BAMLIST=$1
REFERENCE=$2
N_CPUS=$3
CODE_DIR=`dirname $0`
CURR_DIR=`pwd`


#call runPerChrom for all autosomes sequentially  {1..22}
for i in {1..22}
do
    ${CODE_DIR}/runPerChrom.sh $BAMLIST $REFERENCE chr${i} $N_CPUS
done
