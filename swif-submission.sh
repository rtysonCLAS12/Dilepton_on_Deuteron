#!/bin/bash

SRCDIR="/w/work/clas12/tyson/Dilepton_on_Deuteron"

CONFIG_FILE_FO="${SRCDIR}/config.dat"

DESTDIR=$(grep '^OutPath' "$CONFIG_FILE_FO" | cut -d' ' -f2-)
RunLists=$(grep '^RunLists' "$CONFIG_FILE_FO" | cut -d' ' -f2-)
TREENM=$(grep '^treename' "$CONFIG_FILE_FO" | cut -d' ' -f2-)
RCDBPath=$(grep '^RCDBPath' "$CONFIG_FILE_FO" | cut -d' ' -f2-)
CCDBPath=$(grep '^CCDBPath' "$CONFIG_FILE_FO" | cut -d' ' -f2-)

echo
echo "SourceCodePath: $SRCDIR"
echo "OutPath: $DESTDIR"
echo "RunLists: $RunLists"
echo "TreeName: $TREENM"
echo "RCDBPath: $RCDBPath"
echo "CCDBPath: $CCDBPath"

# - Use *your* slurm account
#   Find your name here: https://scicomp.jlab.org/scicomp/slurmJob/slurmAccount
SLURM_ACCT='clas12'  # USE YOUR SLURM ACCOUNT/group here!
WORKFLOW="c12root-$TREENM-$USER"

OUTPUT_GLOB="match:${TREENM}_*"   # output file(s) copied to $DESTDIR by swif

JOB_SCRIPT="${SRCDIR}/job-script.sh"

ECHO="echo"   ## Used to echo swif commands for testing
if [ ${1:-unset} == "submit" ]; then
  ECHO=""  ## disables the 'echo' and actually runs the swif commands
fi

echo
echo "Copy rcdb and ccdb databases to source directory"
cp "$RCDBPath" "${SRCDIR}/"
cp "$CCDBPath" "${SRCDIR}/"

echo
echo "Ensure the DESTDIR ($DESTDIR) is available..."
mkdir -p "$DESTDIR"  ## this will exit and die if the dir doesn't exist and can't be created
echo

#Jobs will use local copy of config file in source directory
#just need global path to set up swif jobs
CONFIG_FILE="config.dat"

# NOTE:
# - The -input file can be pulled directly from tape 
#   the file lands in the working directory on the farm node.  For example:
#      -input 'jpsi_006334.hipo' 'mss:/mss/clas12/rg-b/production/recon/spring2019/torus-1/pass2/v0/dst/train/jpsi/jpsi_006334.hipo' \

$ECHO 
# Loop through each file in RunLists
for runlist_file in $RunLists; do
    
    #runlists have path at top of file
    read -r path < "$runlist_file"
    
    # Read each number in the file and iterate over them
    tail -n +2 "$runlist_file" | while IFS= read -r run_number; do

      #get path to hipo file
      hipo_file_path="${path}${run_number}.hipo"
      #create swif job
      $ECHO \
      swif2 add-job \
        -create \
        -workflow "$WORKFLOW" \
        -account "$SLURM_ACCT" \
        -partition 'production' \
        -constraint 'el9'  \
        -cores '1'  \
        -ram '5000M'  \
        -time '24h'  \
        -shell '/bin/bash'  \
        -output "$OUTPUT_GLOB" "$DESTDIR/" \
        "$JOB_SCRIPT" "${SRCDIR}" "$hipo_file_path" "$run_number" "$CONFIG_FILE" "$TREENM"

    done 
done

# Wait for one minute
sleep 60

# run the swif jobs
$ECHO swif2 run -workflow "$WORKFLOW"

### Some convenient swif2 reports and reminders below
if [ -n "$ECHO" ]; then exit 0; fi

echo
echo "-- Active SWIF workflows --------"
echo swif2 list

echo
echo "-- SWIF status for workflow $WORKFLOW --------"
echo swif2 status -workflow "$WORKFLOW"

echo
echo "-- SWIF jobs for workflow $WORKFLOW --------"
echo swif2 status -workflow "$WORKFLOW" -jobs

echo
echo "-- cancel SWIF workflow -------"
echo swif2 cancel -workflow "$WORKFLOW" -delete

