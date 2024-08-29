#!/bin/bash
#SBATCH --job-name=RNAseq_align      # Name of job
#SBATCH --mail-type=BEGIN,END,FAIL   # Set when I get emails (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=email@.edu # Email address
#SBATCH --partition=sjmac            # Either sixhour or sjmac
#SBATCH --nodes=1                    # Number of nodes to run on
#SBATCH --ntasks=1                   # Number of tasks to run (1 with >1 CPUS per task is multithreaded)
#SBATCH --cpus-per-task=4            # Number of CPUs per task (>1 = multithreaded)
#SBATCH --mem-per-cpu=4gb            # Memory request (default is 2Gb)
#SBATCH --time=12:00:00              # Time limit in hrs:min:sec (default for sjmac queue is 8 hrs)
#SBATCH --output=./run_log_files/salmon_align_%j.log    # Standard output and error logs
#SBATCH --array=1-48

module load salmon/1.10.1
files="data_raw/readname.mapping.txt"

shortname=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f1`

R1_L1=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f 2`
R1_L2=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f 3`
R1_L3=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f 4`
R1_L4=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f 5`

R2_L1=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f 6`
R2_L2=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f 7`
R2_L3=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f 8`
R2_L4=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f 9`


salmon quant -i salmon_output/dmel_index -l A -1 <(cat $R1_L1 $R1_L2 $R1_L3 $R1_L4) -2 <(cat $R2_L1 $R2_L2 $R2_L3 $R2_L4) -p 8 --validateMappings -o salmon_
