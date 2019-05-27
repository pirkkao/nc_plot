#!/bin/bash
#
#SBATCH -J sv_rmse
#SBATCH -t 10:00:00
#SBATCH -n 4
#SBATCH --mem-per-cpu=8000
#
#
module load bioconda python-env/3.5.3 && export CONDA_ENVS_PATH=$WRKDIR/DONOTREMOVE/taito-conda-envs && source activate plot2

exp=sv_rmse

 date=2016120100
edate=2017112600
dstep=192

mandtg=/homeappl/home/ollinaho/bin/mandtg

# Initialize job list
joblist=""

# Loop over dates
while [ $date -le $edate ]; do

    # Iterate over variables
    # T850 - T, 4
    # U850 - U, 4
    # U200 - U, 1
    # Z500 - z, 3
    for vari in "T,4" "U,4" "U,1" "Z,3"; do
	
	# File names
	if [ $vari == "T,4" ]; then
	    namvari=T850
	elif [ $vari == "U,4" ]; then
	    namvari=U850
	elif [ $vari == "U,1" ]; then
	    namvari=U200
	elif [ $vari == "Z,3" ]; then
	    namvari=Z500
	fi

	# Create jobs for GNU parallel
	#
	job=tmp/data.${exp}_${date}_${namvari}

	# Copy base exp config
	cp configs/data.$exp $job
	
	# sed in some field changes
	sed -i "s|2016120100|$date|g" $job
	sed -i "s|T,4|$vari|"         $job

	# create jobs
	joblist="$joblist ${exp}_${date}_${namvari}"
    done
    
    date=$($mandtg $date + $dstep)
done

# What will be run
parallel --dry-run python3 main.py {1} data tmp ::: $joblist

# Execute
parallel --jobs $SLURM_NTASKS python3 main.py {1} data tmp ::: $joblist

# Clean-up
#rm -f tmp/data.${exp}_*
