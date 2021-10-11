#!/bin/bash

PROJECT_DIR="${LSFM_CLUSTER_SCRATCH_USER_PATH}/projects/CRC_STRs/"
SCORE_TYPE="phylo_gap01"
PVALUE=0.05
DIVERGENCE=0.1
MAX_SEQS=100

for i in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}; do
    echo "#!/bin/bash
#SBATCH --job-name=part${i}_refine
#SBATCH --partition=single
#SBATCH --qos=single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node[001-014]
#SBATCH --constraint=skylake-sp
#SBATCH --export=ALL,LSFM_CLUSTER_KEEP_LOCAL_SCRATCH=true
#SBATCH --chdir=/cfs/earth/scratch/verb/jobs/tral/refine/run_dir

module load USS/2020 gcc/7.3.0
module load slurm
module load python
export PYTHONUSERBASE=$LSFM_CLUSTER_SCRATCH_USER_PATH/python
set -x
# make directories in local scratch
mkdir \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python
mkdir -p \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/localapps/bin
mkdir -p \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/localapps/libexec
mkdir -p \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/run/data/
mkdir -p \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/run/results/

rsync -av \$PYTHONUSERBASE/ \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python/
export PYTHONUSERBASE=\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python

rsync -av /cfs/earth/scratch/verb/localapps/bin/ \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/localapps/bin
rsync -av /cfs/earth/scratch/verb/localapps/libexec/mafft \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/localapps/libexec/
export MAFFT_BINARIES=\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/localapps/libexec/mafft

rsync -av /net/home/verb/.tral \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/

for CONFIG in \$(find \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/.tral -name config.ini); do
	sed -i 's|/cfs/earth/scratch/verb/localapps|'\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/localapps'|' \$CONFIG
done
sed -i '/^CONFIG_DIR/s|=.*|= \"'\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/.tral'\"|' \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/python/lib/python3.6/site-packages/tral/paths.py

# copy script + input data to local scratch on node
rsync -av $PROJECT_DIR/python/scoring_filtering/ \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/run/
cp $PROJECT_DIR/data/partitions/part_${i}/part_${i}.fa \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/run/data/

# set local and global output directory variables
LOCAL_OUTPUT=\$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/run/results/
GLOBAL_OUTPUT=$PROJECT_DIR/results/repeats/refined/part_${i}/

# change dir to local scratch on node
cd \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH

# run detection on local scratch
python3 run/tral_refine.py \
-r $PROJECT_DIR/results/repeats/filtered_p0.05_d0.075/part_${i} \
-f \$LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/run/data/part_${i}.fa \
-s ${SCORE_TYPE} \
-p ${PVALUE} \
-d ${DIVERGENCE} \
-m ${MAX_SEQS} \
-g \$GLOBAL_OUTPUT \
-l \$LOCAL_OUTPUT 2> \$LOCAL_OUTPUT/debug.log

# sync local output to global output directory
rsync -av \$LOCAL_OUTPUT \$GLOBAL_OUTPUT
    "> launch_scripts/sbatch_detection_local_part${i}.sh
    # sbatch launch_scripts/sbatch_detection_local_part${i}.sh
    # sleep 10
done