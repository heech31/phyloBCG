#BSUB -J sens_IGsig2_2               # sets the job name to R_s_IGsig2_2.
#BSUB -L /bin/bash           # uses the bash login shell to initialize the job's execution environment.
#BSUB -W 12:00               # sets to 12 hours the job's runtime wall-clock limit.
#BSUB -n 10                  # assigns 10 core for execution.
#BSUB -R "span[ptile=10]"    # assigns 10 core per node.
#BSUB -R "rusage[mem=2500]"  # reserves ~2.5GB per process/CPU for the job (2.5GB * 10 Cores = 25GB per node)
#BSUB -M 2500                # sets to ~2.5GB the per process enforceable memory limit.
#BSUB -o s_IGsig2_2.%J.o           # directs the job's standard output to stdout1.jobid

##OPTIONAL JOB SPECIFICATIONS
#BSUB -P 082757044453        #Set billing account to 082757043477
#BSUB -u hcchung@tamu.edu    #Send all emails to email_address
#BSUB -B -N                  #Send email on job begin (-B) and end (-N)

## Load the necessary modules
module purge
module load R/4.0.0-foss-2020a

## Launch R with proper parameters
Rscript IGsig2_2.R

