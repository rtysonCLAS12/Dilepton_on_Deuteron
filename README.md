# Dilepton_on_Deuteron
Code for the analysis of dilepton production on the deuteron at CLAS12

## clas12root script

The script is called *maketree.cpp*. The goal of this script is to process CLAS12 data and produce an output ROOT tree with relevant variables for the events we are interested in: e d -> e+ e- d.

The script uses [clas12root](https://github.com/JeffersonLab/clas12root/tree/master) and its wrappers to the [RCDB](https://github.com/JeffersonLab/rcdb), [CCDB](https://github.com/JeffersonLab/ccdb), and QADB. clas12root is installed with all dependencies on ifarm and can be accessed using modules. An example on how to set up the right environment can be found on lines 27-45 of the ***job-script.sh***. clas12root uses downloaded local copies of the RCDB and CCDB databases, please refer to the [documentation](https://github.com/JeffersonLab/clas12root/tree/master?tab=readme-ov-file#clas12databases) on how to set these up for clas12root using the [PrepareDatabases.C script](https://github.com/JeffersonLab/clas12root/blob/master/RunRoot/PrepareDatabases.C).

The *maketree.cpp* script uses the config.dat config file to set the paths to the local copies of the RCDB and CCDB prepared as discussed in the documentation. Note that there is a default that points to [my](mailto:tyson@jlab.org) copies, but these might disappear or change over time. Best to have your own. The other two config options are the name of the root tree and whether or not to only use golden runs. Golden runs are based on the Quality Assurance which selects runs where the accumulated charge should be well calculated. The QADB is also used to calculate the charge, this will be printed to screen or to a log file (see below).
***N.B.:*** The database files are copied to the source directory when submitting jobs on the farm, this is to avoid reading from /work too many times. When running interactively, the database files will be opened from the path specified in the config.dat file unless the databases already exist in the working directory.

The point of the *maketree.cpp* script is to produce an output ROOT tree with relevant variables for the events we are interested in. The script already has a bunch of useful variables added to the output. More can be added by doing three things. First, scrolling to the bottom of the script, the function [*void initNames(string * varNames)*](https://github.com/rtysonCLAS12/Dilepton_on_Deuteron/blob/797f0f9ea2b9661317a172492bf3d100bb3232b0/maketree.cpp#L500) creates a list of variable names to be added to the tree. You can add more variables to this list. Scrolling back up to the top of the script, the [*int nVars*](https://github.com/rtysonCLAS12/Dilepton_on_Deuteron/blob/797f0f9ea2b9661317a172492bf3d100bb3232b0/maketree.cpp#L85) variable should be changed to be the same length as the list of variables. Scrolling to the middle of the script is where the branches of the tree are filled with the value corresponding to a variable. For example:

      branchValues[getRow("elP", varNames,nVars)] = el.P();

[fills the branch called elP with the electron momentum](https://github.com/rtysonCLAS12/Dilepton_on_Deuteron/blob/797f0f9ea2b9661317a172492bf3d100bb3232b0/maketree.cpp#L206). Variables corresponding to an electron start with el, those for the positron start with po and those for the deuteron start with deut. Information on accessing bank variables in clas12root is available [here](https://github.com/JeffersonLab/clas12root/blob/master/AccesssingBankDataInCpp.txt), with a description of these variables [here](https://clasweb.jlab.org/wiki/index.php/CLAS12_DSTs). Some variables correspond to an event, eg the invariant mass of the e+ e- pair, or the missing mass of e d -> e+ e- d. The run and event number are also recorded.

The script will iterate over runs one by one, see below on how to run the script. It will read the beam energy from the RCDB so this does not need to be modified. The RCDB also contains information on e.g. the torus polarity.

The script then iterates over all events in a run. The only requirements on these events is that they have at least one electron, one positron and one deuteron as IDed by the [CLAS12 event builder](https://www.sciencedirect.com/science/article/pii/S0168900220300784). The script will take all posible combitorials ie combinations of electrons, positrons and deuterons. That is to say, if they are more than one of each particle type, the code will consider all possible combinations. A flag Combis is used to show which combination an entry corresponds to (ie Combis=0 means first combination in an event). The first electron will be the trigger electron, the first positron and deuteron should be the highest momentum particles, so it's probably safe to use Combis=0 as a starting point, although ultimately some selection algorithm should be implemented.  

The code applies corrections to the electron and positron that radiate photons when going through detector material between the target and detector. These photons can be identified by the fact that they have a small polar angular difference with the electron. The momentum of the photon is added back to the electron or positron. The script also creates a branch containing a flag *elTriangCut* or *poTriangCut* which corresponds to electrons or positrons that pass a useful cut to remove pion contamination. See [Section 6.3 of my thesis](https://www.jlab.org/Hall-B/general/thesis/RTyson_thesis.pdf) for more details on the correction and cut. Finally the script also creates a branch containing a flag *elPassDeadPaddlePCAL* or *poPassDeadPaddlePCAL* which corresponds to whether electrons or positrons hit paddles in the ECAL that were identified as having a low efficiency. Its good to remove these for better matching between data and simulation.


***N.B.:*** The clas12root script can be run interactively with e.g.: 

      clas12root -b 'maketree.cpp("/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass2/v0/dst/train/jpsi/jpsi_006334.hipo","config.dat","eedFS_test.root")' > eedFS_test.log

This creates a log file with the output from the script. Altenatively, the output can be printed to screen with e.g.: 

      clas12root -b 'maketree.cpp("/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass2/v0/dst/train/jpsi/jpsi_006334.hipo","config.dat","eedFS_test.root")'

The script can be run on more than one run using wildcards:

      clas12root -b 'maketree.cpp("/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass2/v0/dst/train/jpsi/jpsi_0063*.hipo","config.dat","eedFS_test.root")'

Running the code interactively should be done for debugging purposes only as it is slow (takes ~10h to process the full RG-B spring2019 dataset on ifarm). I would only run the code interactively on a small number of runs. Two good branches to check in the output tree are RunNb and EventNb. The RunNb should correspond to the hipo file and the EventNb should be non zero.

This code should be updated with [Iguana](https://github.com/JeffersonLab/iguana) once Iguana contains e.g. the lepton PID algorithms or lepton corrections. 

## Plotting script

A script to plot some of the key variables has been included as an example and for convenience to check that the code runs correctly. The script produces one pdf with one plot per page, run it with ROOT. You'll need to change two things:

The treeLocRoot on line 27 should be the path to the output tree produced by the clas12root script.

The fileLoc string variable on line 25 should be the path to where you want the output plots to go. The endName variable on line 29 appends a string to the end of the output file name for the plots (i.e. /fileLoc/Vars_endName.pdf). 

You can further change cuts based on the variables in the trees by changing the cut variable on line 37. cutIM on line 49 is the cut used for the e+ e- invariant mass plots.

## SWIF Submission

Two scripts are used to submit jobs on the JLab farm using swif2. Hopefully you won't need to change these much more than what is indicated below. One final script is used to hadd all output trees.

### job-script.sh 
Sets up a single job. This will be run by each individual farm node. It sets up the environment using modules, copies the code from the source repository onto the farm node then runs the maketree code.

***N.B.:*** Please test job-script.sh before launching swif jobs. To emulate a job on a farm node, create a new directory and test job-script.sh from there using:

      /path/to/job-script.sh /path/to/source/code/eed /path/to/data_runnb.hipo runnb config.data treename

The script should end by printing the contents of the working directory. This should contain code, a log file and an output tree, you can check the contents of this to make sure the script worked.

### swif-submission.sh 
Creates and submits the swif jobs. Three things should be modified in the script: 

On line 3 the SRCDIR variable should contain the path to the source code. 

On line 5, CONFIG_FILE_FO should contain the path to a config file used to setup the swif submission.
This config file must have three entries. OutPath which is the path to the output directory. Note that there will be as many output trees as there are jobs.
RunLists which are files containing the path to the hipo data files and run numbers. You can look at eg spring2019.txt as an example. For convenience there already are three run lists corresponding to the three RG-B datasets.
treename which is how the swif workflow and output files will be named, this can be the same as the tree name in the clas12root script.
Note that this config file can be the same as the one read by the *maketree.cpp* clas12root script as in the example config file provided in this repository.

The third thing to change is the farm account on line 19, this should be clas12 or clas or other. See [documentation here](https://jlab.servicenowservices.com/kb_view.do?sys_kb_id=b022cd801b35c110a888ea4ce54bcb18).

The submission script will copy the database files from the path specified in the config file to the source code directory. This is so the job-script copies these files to the farm nodes to avoid reading from /work too many times. 

The submission script will then submit as many jobs as there are entries in the run list files you provide. The script will also output some helpful commands.

***N.B.:*** Without input parameters the script prints to screen the swif commands instead of submitting the jobs. Do this to check everything looks OK. The script will only actually submit the jobs if you pass *submit* as a parameter eg:

      /path/to/swif-submission.sh submit

There are some test run lists already included in the directory. Please test out the swif-submission script with these. They contain one run number that doesn't exist and one run with no golden events to show you what happens in the code due to these possible trip ups. The code doesn't actually crash but the log files contain useful information to debug empty ouput root files.

The [scicomp webpage](https://scicomp.jlab.org/scicomp/swif/active) shows active swif jobs. Useful documentation on swif can be found [here](https://scicomp.jlab.org/docs/swif2) and [here](https://scicomp.jlab.org/cli/swif.html)

### concat_runs.sh
Is used to hadd all output trees together to only have one output tree with all the processed data once ***all*** of your jobs are finished. It will also cancel your swif workflow and delete individual run trees so wait until ***all*** of your jobs are finished. The script will finally delete the rcdb and ccdb databases from the source directory, this is to avoid any confusion if you change the databases at their original location as specified in the config file. You'll need to change the path to the source code directory and config file on lines 3 and 5. The script uses similar information from the config file as the submission script. It also requires to know the path to the final combined output root file.

## ***Good luck!***


