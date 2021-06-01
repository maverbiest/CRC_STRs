#!/bin/bash

PROJECT_DIR="${LSFM_CLUSTER_SCRATCH_USER_PATH}/projects/CRC_STRs/"
MAX_SEQS=20
# {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}
for i in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}; do
    echo "#!/bin/bash
#SBATCH --job-name=part${i}_scoring
#SBATCH --partition=single
#SBATCH --qos=single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL,LSFM_CLUSTER_KEEP_LOCAL_SCRATCH=true
#SBATCH --chdir=/cfs/earth/scratch/verb/jobs/tral/scoring_filtering/run_dir/

module load USS/2020 gcc/7.3.0
module load slurm
module load python
export PYTHONUSERBASE=$LSFM_CLUSTER_SCRATCH_USER_PATH/python
set -x
# make directories in local scratch
mkdir \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python
mkdir -p \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/data/
mkdir -p \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/results/

rsync -av \$PYTHONUSERBASE/ \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python/
export PYTHONUSERBASE=\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python

rsync -av /net/home/verb/.tral \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/
for CONFIG in \$(find \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/.tral -name config.ini); do
    sed -i -e 's|/cfs/earth/scratch/verb/localapps|'\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/localapps'|' -e 's|cfs/earth/scratch/verb/python/|'\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH'/python/|' -e 's|~/.tral|'\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH'/.tral|' \$CONFIG
done
sed -i '/^CONFIG_DIR/s|=.*|= \"'\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/.tral'\"|' \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python/lib/python3.6/site-packages/tral/paths.py

# change dir to local scratch on node
cd \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH

# run detection on local scratch
python3 $LSFM_CLUSTER_SCRATCH_USER_PATH/projects/CRC_STRs/python/scoring_filtering/tral_pvalues.py \
-d $LSFM_CLUSTER_SCRATCH_USER_PATH/projects/CRC_STRs/results/repeats/partitions/part_${i}
    "> launch_scripts/sbatch_repeat_scoring_part${i}.sh
    # sbatch launch_scripts/sbatch_detection_local_part${i}.sh
    # sleep 60
done