#!/bin/bash

START=10; END=300; PART=class_a
while getopts ":m:s:e:p:" opt; do
   case "$opt" in
      m) MODE=$OPTARG;;
      s) START=$OPTARG;;
      e) END=$OPTARG;;
      p) PART=$OPTARG;;
   esac
done

# Initialise variables
bindir="/home/hpce17/hpce17314/.bin"
scratchdir="/gpfs/scratch/hpce17/hpce17314"
NUM=`date +%T | sed -e "s/://g"`
NAME=$(ls *parm7 | cut -d "." -f 1)
D=$(date +%F)
depend=0

# Copy files to a new scratch folder
if [[ $PWD != *scratch* ]]; then
   mkdir "$scratchdir/MD_$NUM"
   cp ./* $scratchdir/MD_$NUM
   printf "\nCopying files to ~/scratch/MD_$NUM\n\n"
   printf "Files copied from ~${PWD#~} on $D\n" > $scratchdir/MD_$NUM/sub_history
   cd "$scratchdir/MD_$NUM"
fi

# Minimization
if [ $MODE = 'min' ]; then
   cp -n $NAME.rst7 min0.rst7
   cp $bindir/min* .
   START=1; END=6; c=1
   for (( i=$START; i<=$END; i+=1 )); do sed -e "s/orig/$NAME/g" min$i.sh > run$i.sh; done
fi

# Heating 
if [ $MODE = 'heat' ]; then
   cp -n min6.rst7 heat0.rst7
   cp $bindir/heat* .
   START=1; END=2; c=1
   for (( i=$START; i<=$END; i+=1 )); do sed -e "s/orig/$NAME/g" heat$i.sh > run$i.sh; done
fi

# Production
if [ $MODE = 'prod' ]; then
   if [ $START == 10 ]; then cp -n heat2.rst7 prod0ns.rst7; fi
   cp -n $bindir/prod* .
   c=10
   for (( i=$START; i<=$END; i+=10 )); do
      sed -e "s/prod_n/prod${i}ns/g;s/prod_m/prod$((i-10))ns/g" prod.sh > run$i.sh
      sed -i "s/orig/$NAME/g" run$i.sh
   done
fi

for (( i=$START; i<=$END; i+=c )); do
   if [ $i == $START ]; then
      S=$( sbatch -q $PART run$i.sh )
   else
      S=$( sbatch -q $PART --dependency=${JOB[-1]} run$i.sh )
   fi
   JOB+=($( echo $S | awk '{print $NF}' ))
done

#rm run*
sleep 3

for j in ${JOB[@]}; do
   sacct --jobs $j --format=jobid,jobname,qos,nnodes,alloccpus | sed -n 3p >> sub_history
done
