#!/bin/bash
#SBATCH --job-name=reProcessFibroblasts    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lucast3@gene.com     # Where to send mail	
#SBATCH --ntasks=1                # Run on a single CPU
#SBATCH --partition="defq"         #Run on Himem or defq partition
#SBATCH --mem=220gb                     # Job memory request
#SBATCH --cpus-per-task=21            # Number of CPU cores per task
#SBATCH --output=/gstore/scratch/u/lucast3/fibroticmemory/reprocessOut.json   # Standard output and error log

pwd; hostname; date

source activate MACS

ml R/prd

Rscript -e "rmarkdown::render('/gstore/scratch/u/lucast3/fibroticmemory/reProcess.Rmd')"
