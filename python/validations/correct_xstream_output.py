#!/usr/bin/env python3

import argparse
from collections import defaultdict
import os
import sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import pandas as pd

from scoring_filtering.tral_pvalues import load_repeatlists

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mismatches", "-m", type=str, required=True, help="tsv file containing information on which STRs do not match the reference sequence"
    )
    parser.add_argument(
        "--tr_directory", "-t", type=str, required=True, help="Directory containing pickled RepeatLists to validate against reference sequence"
    )

    return parser.parse_args()

def correct_xstream_repeat(repeat, correction_seq):
    new_msa = []
    character_counter = -1
    for unit in repeat.msa:
        new_unit = ""
        for nucleotide in unit:
            if nucleotide != "-":
                character_counter += 1
            if nucleotide == "X":                
                nucleotide = correction_seq[character_counter]
            new_unit += nucleotide
        new_msa.append(new_unit)

    repeat.msa = new_msa


def main():
    args = cla_parser()

    df = pd.read_csv(args.mismatches, sep="\t")
    target_genes = set(df.seq_of_origin)

    # remove the two mitochondrial genes that have different index offset due to start at 0
    exceptions = {"chrM_0-4262_fw", "chrM_0-5511_fw"}
    target_genes.difference_update(exceptions)

    partition_dirs = [os.path.join(args.tr_directory, i) for i in os.listdir(args.tr_directory)\
         if i.startswith("part_")\
             and os.path.isdir(os.path.join(args.tr_directory, i))]
    
    for partition in partition_dirs:
        for file_name, repeat_list in load_repeatlists(partition, targets=target_genes):
            for index, row in df.loc[(df["seq_of_origin"] == file_name.replace(".pickle", ""))].iterrows():     
                if not row["msa"] == "[]" and row["TRD"] == "XSTREAM":
                    for repeat in repeat_list.repeats:
                        if repeat.TRD == "XSTREAM" and repeat.begin == row["begin"]:
                            correct_xstream_repeat(repeat, row["ref_seq"])
                            break                    
                else:
                    for repeat in repeat_list.repeats:
                        if repeat.TRD == row["TRD"] and repeat.begin == row["begin"]:
                            repeat.msa = repeat.text.split(" ")
                            break                                       
            repeat_list.write(output_format="pickle", file=os.path.join(partition, file_name))

if __name__ == "__main__":
    main()
