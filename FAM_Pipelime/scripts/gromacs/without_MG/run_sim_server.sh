#!/bin/bash

usage() { echo "Run a single replicate MD simulation
Usage: single_run.sh -c <structure file (.gro)> -d <mdp directory> -a <max. warnings for grompp> (optional) -p <plumed file (.dat)>" 1>&2; exit 1; }
invalidOpt() { echo "Invalid option: -$OPTARG" 1>&2; exit 1; }
missingArg() { echo "Option -$OPTARG requires an argument" 1>&2; exit 1; }
cleanup() { if ls -f $1/\#* 1> /dev/null 2>&1 ; then rm $1/\#* ; fi ; }

while getopts ":c:d:a:p:h" opt; do
    case $opt in
        c) 
            structureFile=$OPTARG
            ;;
        d)
            mdp_dir=$OPTARG
            ;;
        p)  
            plumedFile=$OPTARG
            ;;
        h)
            usage
            ;;
        a)
            maxwarn=$OPTARG
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

structureName=`echo $structureFile | rev | cut -f1 -d"/" | rev | cut -f1 -d"."`
    gmx mdrun -v -s md0/"$structureName".tpr -c md0/"$structureName".gro -x md0/"$structureName".xtc -cpo md0/"$structureName".cpt -e md0/"$structureName".edr -g md0/"$structureName".log || { echo "-> Error: gmx mdrun for md0 failed" ; cleanup md0 ; exit 1; }

