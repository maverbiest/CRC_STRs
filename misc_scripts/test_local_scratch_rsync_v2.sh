#!/bin/bash
#SBATCH --job-name=test_tral-GLOBAL_python_localapps_LOCAL_OUTPUT
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
mkdir $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python
mkdir -p $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/data/
mkdir -p $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/results/
# mkdir -p $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH//CRC_STRs/results/repeats/partitions

# output directory variables
LOCAL_OUTPUT=$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/results/
GLOBAL_OUTPUT=/cfs/earth/scratch/verb/projects/CRC_STRs/results/repeats/partitions/test
# GLOBAL_LOCAL_OUTPUT=$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/results/repeats/partitions/test/
# rsync -av $GLOBAL_OUTPUT $GLOBAL_LOCAL_OUTPUT
rsync -av $PYTHONUSERBASE/ $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python/
export PYTHONUSERBASE=$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python

rsync -av /cfs/earth/scratch/verb/localapps $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/
rsync -av /net/home/verb/.tral $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/
for CONFIG in $(find $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/.tral -name config.ini); do
	sed -i 's|/cfs/earth/scratch/verb/localapps|'$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/localapps'|' $CONFIG
done
sed -i '/^CONFIG_DIR/s|=.*|= "'"$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH"'/.tral"|' $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python/lib/python3.6/site-packages/tral/paths.py

# copy script + input data to local scratch on node
cp /cfs/earth/scratch/verb/projects/CRC_STRs//python/tral_detector_run_local.py $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python/
cp /cfs/earth/scratch/verb/projects/CRC_STRs/data/test/fasta/tiny.fa $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/data/

#rsync -av /cfs/earth/scratch/verb/CRC_STRs//results/repeats/partitions/test /cfs/earth/scratch/verb/CRC_STRs/results/repeats/partitions/

# change dir to local scratch on node
cd $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH

# run detection on local scratch
# python3 python/tral_detector_run_local.py -f $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/data/tiny.fa -m 20 -g $GLOBAL_LOCAL_OUTPUT -l $LOCAL_OUTPUT 2> $LOCAL_OUTPUT/debug.log
python3 python/tral_detector_run_local.py -f $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/CRC_STRs/data/tiny.fa -m 20 -g $GLOBAL_OUTPUT -l $LOCAL_OUTPUT 2> $LOCAL_OUTPUT/debug.log

# sync local output to global output directory
# rsync -av $GLOBAL_LOCAL_OUTPUT $GLOBAL_OUTPUT
# rsync -av --exclude='repeats' $LOCAL_OUTPUT $GLOBAL_OUTPUT
rsync -av $LOCAL_OUTPUT $GLOBAL_OUTPUT
