#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=100M
#SBATCH --output=cpp_sanity_check.out
#SBATCH -p batch-csl,batch-skl

cd /scratch/cs/networks/nurmit7/mesu/
echo "here 1"
module load anaconda
echo "here 2"
python cpp_sanity_check.py [2,2] cpp_sanity_check_results_2_2/