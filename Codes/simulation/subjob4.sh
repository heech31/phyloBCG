#BSUB -J R_se                # sets the job name to R_se (SpiecEasi).
#BSUB -L /bin/bash           # uses the bash login shell to initialize the job's execution environment.
#BSUB -W 12:00               # sets to 5 hours the job's runtime wall-clock limit.
#BSUB -n 10                  # assigns 10 core for execution.
#BSUB -R "span[ptile=10]"    # assigns 10 core per node.
#BSUB -R "rusage[mem=1000]"  # reserves ~1.0GB per process/CPU for the job (1.0GB * 10 Cores = 10GB per node)
#BSUB -M 1000                # sets to ~1GB the per process enforceable memory limit.
#BSUB -o spieceasi.%J.o           # directs the job's standard output to stdout1.jobid

##OPTIONAL JOB SPECIFICATIONS
#BSUB -P 082757044453        #Set billing account to 082757043477
#BSUB -u hcchung@tamu.edu    #Send all emails to email_address
#BSUB -B -N                  #Send email on job begin (-B) and end (-N)

## Load the necessary modules
module purge
module load R/4.0.0-foss-2020a

## Launch R with proper parameters
Rscript m4.R

