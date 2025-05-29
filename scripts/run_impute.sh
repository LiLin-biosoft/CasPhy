#!/bin/bash

SCRIPT_PATH=$(dirname $0)
while getopts "p:n:b:t:c:" opt
do
   case $opt in
     p) PREFIX=$OPTARG ;;
     n) NEIGHBOR=$OPTARG ;;
     b) BATCHSIZE=$OPTARG ;;
     t) NTHREAD=$OPTARG ;;
     c) CONTINUE=$OPTARG ;;
   esac
done

if [ $CONTINUE == 0 ];then
  echo "run_mode	2" >${PREFIX}_status.txt
  echo "n_iter	1" >>${PREFIX}_status.txt
  echo "n_neighbor	$NEIGHBOR" >>${PREFIX}_status.txt
  echo "p_mutation	0.5" >>${PREFIX}_status.txt
  echo "batch_size	$BATCHSIZE" >>${PREFIX}_status.txt

  Rscript --vanilla ${SCRIPT_PATH}/run_impute.R \
	  ${PREFIX}.csv \
	  ${PREFIX}_iter1.csv \
	  ${PREFIX}_status.txt \
	  ${NTHREAD} \
	  ${PREFIX}_range
  if [ $? != 0 ];then
    exit 1
  fi

fi

while [ $(grep n_neighbor ${PREFIX}_status.txt|gawk '{print $2}') -gt 1 ];do
  NITER=$(grep n_iter ${PREFIX}_status.txt|gawk '{print $2}')
  Rscript --vanilla ${SCRIPT_PATH}/run_impute.R \
	  ${PREFIX}_iter$[$NITER-1].csv \
	  ${PREFIX}_iter${NITER}.csv \
	  ${PREFIX}_status.txt \
	  ${NTHREAD} \
	  ${PREFIX}_range
  if [ $? != 0 ];then
    exit 1
  fi
  rm ${PREFIX}_iter$[$NITER-1].csv
done
mv ${PREFIX}_iter$NITER.csv ${PREFIX}_imputed.csv
rm ${PREFIX}_range.bk ${PREFIX}_range.desc
rm ${PREFIX}_status.txt
