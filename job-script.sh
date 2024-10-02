#!/bin/bash

source /home/tyson/.bashrc

## Test command and ensure job-script was called with the runnb arg we want
if [ $# != 5 ]; then
  echo "Wrong number of arguments"
  echo "  Command as-called: '$0 $*'"
  echo "  Command should be: '$0 <source code repository> <datafilepath> <runnb> <configfile> <treename>'"
  exit 1
fi

SRCDIR=$1
FILEPATH=$2
RUNNB=$3
CONFIGFILE=$4
TREENM=$5

#to test this script, create new directory somewhere other than source code
#this mimics copying working directory to farm node
#then use command:
#/path/to/source/directoy/job-script.sh /path/to/source/directoy/ /cache/clas12/rg-b/production/recon/spring2019/torus-1/pass2/v0/dst/train/jpsi/jpsi_006334.hipo 006334 config.dat eed

#set -x   ## Dump all commands to stderr for extra debugging info
set -e   ## Set bash to die on errors (better this then waste time on a broken job!)

module use /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/modulefiles
#module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles

module purge
module load root/6.30.04
module load sqlite/5.10

#clas12root stuff

module load clas12root/1.8.4
export CLAS12ROOT=/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.8.4/4.1.0/
export ROOT_INCLUDE_PATH=/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.8.4/4.1.0/hipo4:${ROOT_INCLUDE_PATH}
export ROOT_INCLUDE_PATH=/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.8.4/4.1.0/Clas12Banks:${ROOT_INCLUDE_PATH}
export ROOT_INCLUDE_PATH=/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.8.4/4.1.0/Clas12Root:${ROOT_INCLUDE_PATH}

export  CC=/usr/bin/gcc
export  CX=/usr/bin/g++

export PATH="$PATH":"$CLAS12ROOT/bin"

if [[ "$PWD" == *"/home/"* ]]; then
  echo "Job working dir ($PWD) is under /home.  You probably don't want this.  Exitting..."
  exit 1
fi

echo
echo "-- Copy your software to $PWD on the local node -------------------------------"
rsync -avP --exclude "output/" --exclude ".git*" "${SRCDIR}/" "$PWD"  # don't copy output/ dir or git stuff

echo
echo "-- Run main job -----------------------------------"
set +e   # we'll disable 'killing the shell on errors' from here on

# NOTE: Be mindful that any output files *must* have a unique name, or they will
#       step on each other when you run multiple jobs.
#       Here we rename the output based on the run number passed as input

OUTPUTFILE="${TREENM}_${RUNNB}.root"
LOGFILE="${TREENM}_${RUNNB}.log"
clas12root -b 'maketree.cpp("'${FILEPATH}'", "'${CONFIGFILE}'", "'${OUTPUTFILE}'")' > "${LOGFILE}"


echo
echo "-- List Working Directory (END OF JOB) ---------------------------------------------"
echo ${PWD}
ls -laF .

