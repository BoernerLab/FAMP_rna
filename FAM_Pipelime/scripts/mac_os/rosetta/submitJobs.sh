#!/bin/bash

# cmd parsing functions
usage() { echo "submit a Rosetta FARFAR job
Usage: submitJobs -i <FARFAR input script> [-d <directory>] [-p <number of processors>]" 1>&2; exit 1; }
invalidOpt() { echo "Invalid option: -$OPTARG" 1>&2; exit 1; }
missingArg() { echo "Option -$OPTARG requires an argument" 1>&2; exit 1; }

# default cmd arguments
dir=out
proc=1


#------------
# cmd parsing
#------------

while getopts ":i:d:p:h" opt; do
    case $opt in
        i)
            farfar=$OPTARG
            ;;
        d)
            dir=$OPTARG
            ;;
        p)
            proc=$OPTARG
            ;;
        h)
            usage
            ;;
        \?)
            invalidOpt
            ;;
        :)
            missingArg
            ;;
        *)
            usage
            ;;
    esac
done

# no cmd line arguments given
if [ -z "$farfar" ]; then
    usage
fi

# folder exists
if [ -d "$dir" ]; then
    echo "the specified output folder already exists"
    #exit 1
fi


#-----------
# Submission
#-----------

# create output directories
mkdir ./"$dir"/out
for i in `seq 1 $proc`; do
    mkdir ./"$dir"/out/"$i"
done
echo "directories for $proc processes created in ./$dir/out."
sleep 0.5

# create submit script for each processor
for i in `seq 1 $proc`; do
    cat ./$farfar
    sleep 2
    file=`cat ./$farfar | sed "s/-silent /-silent .\/$dir\/out\/$i\//g"`
    if [ ! -d "$dir/submit" ]; then
        mkdir ./$dir/submit
    fi
    echo "$file" > ./$dir/submit/job"$i".sh
    #echo "-run:use_time_as_seed" >> ./$dir/submit/job"$i".sh
    chmod +x ./$dir/submit/job"$i".sh
done

# submit jobs
for i in `seq 1 $proc`; do
    ./$dir/submit/job"$i".sh &
    echo "job $i started..."
    sleep 1.5
done


# kill all running jobs by issuing
#pkill rna_*
