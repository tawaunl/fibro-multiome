#!/bin/bash
#SBATCH --job-name=Ulm    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lucast3@gene.com     # Where to send mail
#SBATCH --ntasks=1                # Run on a single CPU
#SBATCH --partition="defq"         #Run on Himem or defq partition
#SBATCH --mem=120gb                     # Job memory request
#SBATCH --cpus-per-task=7            # Number of CPU cores per task
#SBATCH --output=/gstore/scratch/u/lucast3/fibroMultiome/CollectTriOut.json   # Standard output and error log

pwd; hostname; date

Rscript /gstore/scratch/u/lucast3/fibroMultiome/AnalysisScripts/CollecTRI_ULM.R
