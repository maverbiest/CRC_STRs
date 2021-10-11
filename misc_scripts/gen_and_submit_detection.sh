#!/bin/bash

project_dir="${LSFM_CLUSTER_SCRATCH_USER_PATH}/projects/CRC_STRs/"

for i in `seq 1 25`; do
    echo "#!/usr/bin/bash
#SBATCH --job-name=part${i}_tral_detect
#SBATCH --partition=single
#SBATCH --qos=single
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=node009
#SBATCH --chdir=${LSFM_CLUSTER_SCRATCH_USER_PATH}/jobs/tral/detection/run_dir/

module load USS/2020 gcc/7.3.0
module load slurm
module load python

python3 ${project_dir}python/tral_detector_run.py \
-f ${project_dir}data/partitions/part_${i}/part_${i}.fa \
-o ${project_dir}results/repeats/partitions/part_${i}/ \
2> ${project_dir}results/repeats/partitions/part_${i}/debug.log" > launch_scripts/sbatch_detection_part${i}.sh
    sbatch launch_scripts/sbatch_detection_part${i}.sh
    sleep 1
done
