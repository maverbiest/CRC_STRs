#!/bin/bash

PROJECT_DIR="${LSFM_CLUSTER_SCRATCH_USER_PATH}/projects/CRC_STRs/"
MAX_SEQS=20

for i in {2,3,4,13,14,16,24}; do
    echo "#!/bin/bash
#SBATCH --job-name=part${i}_tral
#SBATCH --partition=single
#SBATCH --qos=single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --constraint=skylake-sp
#SBATCH --export=ALL,LSFM_CLUSTER_KEEP_LOCAL_SCRATCH=true
#SBATCH --chdir=/cfs/earth/scratch/verb/jobs/tral/local_scratch/run_dir/

module load USS/2020 gcc/7.3.0
module load slurm
module load python
export PYTHONUSERBASE=$LSFM_CLUSTER_SCRATCH_USER_PATH/python
set -x
# make directories in local scratch
mkdir \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python
mkdir -p \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/data/
mkdir -p \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/results/

# output directory variables
LOCAL_OUTPUT=\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/results/
GLOBAL_OUTPUT=$PROJECT_DIR/results/repeats/partitions/part_${i}/

rsync -av \$PYTHONUSERBASE/ \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python/
export PYTHONUSERBASE=\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python

rsync -av /cfs/earth/scratch/verb/localapps \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/
rsync -av /net/home/verb/.tral \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/
for CONFIG in \$(find \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/.tral -name config.ini); do
	sed -i 's|/cfs/earth/scratch/verb/localapps|'\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/localapps'|' \$CONFIG
done
sed -i '/^CONFIG_DIR/s|=.*|= \"'\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/.tral'\"|' \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python/lib/python3.6/site-packages/tral/paths.py

# copy script + input data to local scratch on node
cp /cfs/earth/scratch/verb/projects/CRC_STRs//python/tral_detector_run_local.py \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python/
cp $PROJECT_DIR/data/partitions/part_${i}/part_${i}.fa \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/data/

# change dir to local scratch on node
cd \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH

# run detection on local scratch
python3 python/tral_detector_run_local.py -f \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/data/part_${i}.fa -m ${MAX_SEQS} -g \$GLOBAL_OUTPUT -l \$LOCAL_OUTPUT 2> \$LOCAL_OUTPUT/debug.log

# sync local output to global output directory
rsync -av \$LOCAL_OUTPUT \$GLOBAL_OUTPUT
    "> launch_scripts/sbatch_detection_local_part${i}.sh
    sbatch launch_scripts/sbatch_detection_local_part${i}.sh
    sleep 10
done