#!/usr/bin/env bash

DATABASE_PATH="/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/db/test_brca2_M.db"
GTF_PATH="/cfs/earth/scratch/verb/projects/CRC_STRs/data/test/genome_annot/gencode_small_brca2_M.gtf"
REPEAT_DIR="/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/db/test_data"

echo "Generating DataBase here: ${DATABASE_PATH}";
python3 /cfs/earth/scratch/verb/projects/CRC_STRs/python/db_utils/setup_db.py \
--database ${DATABASE_PATH} && 
echo "DataBase setup done" &&

python3 /cfs/earth/scratch/verb/projects/CRC_STRs/python/db_utils/gtf_to_sqlite.py \
--database ${DATABASE_PATH} \
--gtf ${GTF_PATH} &&
echo "GTF succesfully inserted into DataBase" &&

python3 /cfs/earth/scratch/verb/projects/CRC_STRs/python/db_utils/insert_repeats.py \
--database ${DATABASE_PATH} \
--repeat_dir ${REPEAT_DIR} \
--score_type "phylo_gap01" &&
echo "Repeats succesfully inserted into DataBase" &&

echo "DataBase setup succesfully and populated!"
