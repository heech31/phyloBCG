#BSUB -J R_QMP3              # sets the job name to R_fix.
#BSUB -L /bin/bash           # uses the bash login shell to initialize the job's execution environment.
#BSUB -W 10:00                # sets to 10 hours the job's runtime wall-clock limit.
#BSUB -n 2                   # assigns 2 core for execution.
#BSUB -R "span[ptile=2]"     # assigns 2 core per node.
#BSUB -R "rusage[mem=4000]"  # reserves ~4.0GB per process/CPU for the job (4.0GB * 2 Cores = 8GB per node)
#BSUB -M 4000                # sets to ~4GB the per process enforceable memory limit.
#BSUB -o qmp3.%J.o           # directs the job's standard output to stdout1.jobid

##OPTIONAL JOB SPECIFICATIONS
#BSUB -P 082757044453        #Set billing account to 082757043477
#BSUB -u hcchung@tamu.edu    #Send all emails to email_address
#BSUB -B -N                  #Send email on job begin (-B) and end (-N)

## Load the necessary modules
module purge
module load R/4.0.0-foss-2020a

## Launch R with proper parameters
Rscript qmp3.R

