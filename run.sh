#!/bin/bash
#SBATCH --qos=boucher-b
#SBATCH --partition=hpg2-compute
#SBATCH --job-name=dis1en31
##SBATCH --mail-user=kingdgp@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output edges/dis1en31_%A-%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=9gb
#SBATCH --time=4-00:00:00
#SBATCH --array=1-1217

/usr/bin/time -v ./vari_find_paths $SLURM_ARRAY_TASK_ID 90000 1000 31.reads.dbg restriction_nodes1en31pe
