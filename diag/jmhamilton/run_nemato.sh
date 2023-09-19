#!/bin/bash
#SBATCH -t 360
#SBATCH --qos=regular
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH -A m3016
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH -J nemato
#SBATCH -C cpu
#SBATCH -o nemato.out
#SBATCH -e nemato.err
##SBATCH -L scratch
#SBATCH --mail-user=jmhamilton@lanl.gov
#SBATCH --export=ALL

#################
#### READ ME ####
# *** ENSURE No other files named ss* are in the current directory as they will be deleted ***
# This SLURM script will form an array of all poincare bin files in the current directory, then send them to Nemato to process in parallel.
# The poincare bin files are created for many timesteps by running pixplot after a Pixie3D simulation, with output timesteps specified by the -dplot parameter when running pixplot.
# This script will request a standard allocation for the maximum duration (minus 2 minutes for proper exiting).
# The USER should change the --nodes= value to reflect the requested resources that are necessary for the resolution and number of files.
# Use option -h or --help for a read-out of possible options (i.e. ./run_nemato.sh -h).
#################
#### TO DO ######
# Test if pixplot functionality works.
# Change module files to Chicoma GNU environment
#################

date

# Check for user options
if [ $# -eq 0 ]; then
	echo No options provided, defaulting to 8 files per node and pixplot will not run
	echo "(run $0 -h for help)"
	echo ""
fi

# Default option variables
filespernode=8
eqdsk_file=""
pixplot=false
dplot=10.

while [ True ]; do
	if [ "$1" = "-e" -o "$1" = "--eqdsk" ] && [ ! -z "$2" ]; then
		echo "Attempting to use EQDSK file at the relative Path: $2"
		eqdsk_file=$2
		shift 2
	elif [ "$1" = "-f" -o "$1" = "--filespernode" ] && [ ! -z "$2" ]; then
		echo "Attempting $2 file(s) per node"
		filespernode=$2
		shift 2
	elif [ "$1" = "-p" -o "$1" = "--pixplot" ] && [ ! -z "$2" ]; then
		echo "Attempting to first run pixplot, with -dplot = $2"
		pixplot=true
		dplot=$2
		shift 2	
	elif [ "$1" = "-h" -o "$1" = "--help" ]; then
		echo "Usage:"
		echo "sbatch ./run_nemato.sh -h"
		echo "sbatch ./run_nemato.sh -e </relative/path/eqdsk_file> -f <filespernode> -p <dplot>"
		echo ""
		echo "(optional)	-e, --eqdsk	/relative/path/eqdsk_file	: to specify the relative Path to the eqdsk file, default is none"
		echo "(optional)	-f, --filespernode	filespernode	: to specify how many files will be processed per node, default is 8"
		echo "(optional)	-p, --pixplot	dplot	: to first run pixplot with timesteps of dplot (Alfven times) to generate poincare.bin files, default is false and dplot=10."
		echo "(optional)	-h, --help	: help (this message)"
		exit
	elif [ -z "$1" ]; then
		break
	else
		echo ""
		echo User provided a bad option that will be ignored
		echo "(run $0 -h for help)"
		shift 1
	fi
done

# Ensure input arguments have the proper type, trim leading zeros, turn dplot into a float, etc
if [[ $filespernode =~ ^[0-9]+$ ]] && (( filespernode > 0)); then
	filespernode=$(echo $filespernode | sed 's/^0*//')
	echo ""
	echo "Using filespernode = $filespernode"
else	
	echo ""
        echo filespernode must be a positive integer, defaulting to 8 files per node
        filespernode=8
fi
  
if ! [ -r "$eqdsk_file" ]; then
	echo ""
        echo Either no EQDSK file was specified or it did not exist or it was not readable, so none will be used
        eqdsk_file=""
else
	echo ""
	echo Successfully found the specified EQDSK file
fi 

if [ "$pixplot" = true ]; then
	if [[ $dplot =~ ^[0-9]+$ ]] || [[ $dplot =~ ^[+]?[0-9]+\.?[0-9]*$ ]]; then
		dplot=$(echo $dplot | sed 's/^0*//')
		if [[ $dplot =~ ^[0-9]+$ ]]; then
			dplot="${dplot}.0"
		fi
		echo ""
		echo "Using dplot = $dplot"
	else
		echo ""
		echo dplot must be a positive number, defaulting to 10. Alfven times
		dplot=10.
	fi
fi

########### Running Pixplot ###########

if [ "$pixplot" = true ]; then
	module purge
	source ~/bin/load_gcc-12.1.0_chicoma
	pixplotlocation=$(ls ./pixplot*)
	$pixplotlocation -dplot $dplot &
	wait
else
	echo ""
	echo Pixplot option turned off
	echo ""
fi

# Specify line-breaks as array separators, anticipating output format of the 'ls' command
IFS=$'\n'
# Form an array consisting of the locations of all desired files
totalarray=($(ls poinc_t=*.bin))   
# Return IFS to default 
unset IFS
echo "This script will process ${#totalarray[@]} poincare.bin files"

########### Running Nemato ###########

#module purge
# Load appropriate modules based on compiled nemato.mpi.x
#source ~/bin/load_gcc-9.3.0_badger  
#module load git/2.31.1   
#module load mercurial/5.8
#module load gcc/10.3.0   
#module load openmpi/3.1.6
#module load hdf5-parallel/1.10.7
#source ~/bin/load_gcc-12.1.0_chicoma

export OMP_NUM_THREADS=2

# Remove old .bin files used by nemato for plotting, ***ENSURE No other desired files are similarly named in the current directory***
if [ $(ls -l ss* | wc -l) -gt 0 ]; then
        rm ss*
        echo All old files have been deleted
fi

# Form an array to be used to index these files
filesindex=($(seq 0 1 $(($filespernode-1))))
# Initialize the segmented array that will specify which files are to be processed by Nemato in each node
filesarray=()
# Initialize a loop counter
counter=0
# Loop through all files
while [ $counter -lt ${#totalarray[@]} ]; do
        # Re-initialize the segmented array every iteration
        filesarray=()
        # Append the specified number of files to the segmented array to be processed by Nemato, making sure we don't index out of bounds
        for filenumber in ${filesindex[@]}; do
                if [ $(($counter + $filenumber)) -lt ${#totalarray[@]} ]; then
                        filesarray+="${totalarray[$counter + $filenumber]} "
                fi
        done
        # This is the main MPI Nemato command that operates on the specified files
	#sleep $(( $RANDOM % 5 + 1 )) && echo ${filesarray[@]} &
	srun -N 1 --ntasks=8 ./nemato.mpi.x ${filesarray[@]} &
	#mpiexec -np 4 ./nemato.mpi.x ${filesarray[@]} &
        # Update the counter to keep track of how many files have been processed
        counter=$(($counter + $filespernode))
done
# Wait until all background commands are complete before progressing further
wait

########### Generating images and a movie ##########

# Now execute poinc.py to generate .png image files and a poincare.mp4 movie file
#module load python
#module load sandbox ffmpeg
python poinc.py $eqdsk_file

echo ""
echo all done!
