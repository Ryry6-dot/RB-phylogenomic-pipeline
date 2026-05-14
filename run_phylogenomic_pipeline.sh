#!/bin/bash
#SBATCH --job-name=phylogenomic_pipe_RB_3
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:0:0
#SBATCH --output=phylogenomic_pipe_W6_Proj_3.out
#SBATCH --error=phylogenomic_pipe_W6_Proj_3.err
#SBATCH --mail-user=boldajir@oregonstate.edu
#SBATCH --mail-type=END

#run phylogenomic pipeline
python phylogenomic_pipeline.py
 