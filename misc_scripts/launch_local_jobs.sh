#!/bin/bash

PROJECT_DIR="${LSFM_CLUSTER_SCRATCH_USER_PATH}/projects/CRC_STRs/"
MAX_SEQS=10


for i in {2,3,4,7,13,14,15,16,17,20,24,25}; do
    echo "#!/bin/bash
#SBATCH --job-name=part${i}_tral
#SBATCH --partition=single
#SBATCH --qos=single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --constraint=skylake-sp
#SBATCH --export=ALL,LSFM_CLUSTER_KEEP_LOCAL_SCRATCH=false

module load USS/2020 gcc/7.3.0
module load slurm
module load python

# make directories in local scratch
mkdir \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python
mkdir -p \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/data/
mkdir -p \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/results/

# output directory variables
LOCAL_OUTPUT=\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/results/
GLOBAL_OUTPUT=$PROJECT_DIR/results/repeats/partitions/part_${i}/

# copy script + input data to local scratch on node
cp $PROJECT_DIR/python/tral_detector_run_local.py \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python/
cp $PROJECT_DIR/data/partitions/part_${i}/part_${i}.fa \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/data/

# change dir to local scratch on node
cd \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH

# run detection on local scratch
python3 python/tral_detector_run_local.py \
-f \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/data/part_${i}.fa \
-m $MAX_SEQS \
-g \$GLOBAL_OUTPUT \
-l \$LOCAL_OUTPUT 2> \$LOCAL_OUTPUT/debug.log

# sync local output to global output directory
rsync -av \$LOCAL_OUTPUT \$GLOBAL_OUTPUT" > launch_scripts/sbatch_detection_local_part${i}.sh
    sbatch launch_scripts/sbatch_detection_local_part${i}.sh
    sleep 1
done
