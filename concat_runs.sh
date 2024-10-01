#!/bin/bash

SRCDIR="/w/work/clas12/tyson/Dilepton_on_Deuteron"

CONFIG_FILE_FO="${SRCDIR}/config.dat"

ROOT_DIR=$(grep '^OutPath' "$CONFIG_FILE_FO" | cut -d' ' -f2-)
OUTPUT_FILE=$(grep '^OutputFile' "$CONFIG_FILE_FO" | cut -d' ' -f2-)
TREENM=$(grep '^treename' "$CONFIG_FILE_FO" | cut -d' ' -f2-)

WORKFLOW="c12root-$TREENM-$USER"
echo
echo "cancel swif workflow $WORKFLOW"
swif2 cancel -workflow "$WORKFLOW" -delete

echo
echo "hadd *.root from: $ROOT_DIR"
echo "Output File: $OUTPUT_FILE"

# Initialize an empty file list
file_list=""

# Iterate over each .root file in the specified directory
for root_file in "$ROOT_DIR"/*.root; do
    if [[ -f "$root_file" ]]; then  # Check if it is a regular file
        file_list+="$root_file "     # Add the file to the list
    fi
done

echo
echo "Combining files:"
echo $file_list
echo

# Use hadd to merge the files
if [[ -n $file_list ]]; then
    hadd "$OUTPUT_FILE" $file_list
    echo "Merged files into $OUTPUT_FILE"

    
else
    echo "No .root files found in $ROOT_DIR"
fi

rm $ROOT_DIR/*.root
echo "Removed original .root files."
    

