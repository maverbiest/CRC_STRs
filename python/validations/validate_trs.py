#!/usr/bin/env python3

import argparse
import os
import sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

from scoring_filtering.tral_pvalues import load_repeatlists

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--reference", "-r", type=str, required=True, help="Path to the reference genome that the TRs will be validated against"
    )
    parser.add_argument(
        "--tr_directory", "-t", type=str, required=True, help="Directory containing pickled RepeatLists to validate against reference sequence"
    )

    return parser.parse_args()


def main():
    args = cla_parser()

    chromosomes = {f"chr{i}": "" for i in range(1, 23)}
    chromosomes = {**chromosomes, **{"chrX": "", "chrY": "", "chrM": ""}}

    mismatch_count = 0
    fasta_parser = SeqIO.parse(args.reference, "fasta")
    for record in fasta_parser:
        try:
            chromosomes[record.id] = record
        except KeyError:
            continue
    
    mismatch_dict = {
        "seq_of_origin": [],
        "TRD": [],
        "begin": [],
        "repeat_region_length": [],
        "msa": [],
        "ref_seq": []
    }

    partition_dirs = [os.path.join(args.tr_directory, i) for i in os.listdir(args.tr_directory)\
         if i.startswith("part_")\
             and os.path.isdir(os.path.join(args.tr_directory, i))]

    for partition in partition_dirs:
        for file_name, repeat_list in load_repeatlists(partition):
            chrom_of_origin = file_name.split("_")[0]

            if "_fw" in file_name:
                strand = "fw"
            elif "_rv" in file_name:
                strand = "rv"                
            else:
                raise Exception(f"Could not determine strand for repeatlist from file {file_name}")
                        
            for repeat in repeat_list.repeats:   
                if strand == "fw":
                    # only valid for fw strand genes
                    gene_begin = int(file_name.split("_")[1].split("-")[0])
                    chrom_begin = gene_begin + repeat.begin - 2 
                else:
                    # only valid for rv strand genes
                    gene_begin = int(file_name.split("_")[1].split("-")[0])
                    gene_end = int(file_name.split("_")[1].split("-")[1])
                    repeat_end = repeat.begin + repeat.repeat_region_length - 1
                    chrom_begin = gene_end - repeat_end 

                repeat_string = "".join(repeat.msa).replace("-", "")    
                chrom_string = str(chromosomes[chrom_of_origin].seq[chrom_begin : chrom_begin + repeat.repeat_region_length])
                if strand == "rv":
                    chrom_string = str(Seq(chrom_string).reverse_complement())

                if not repeat_string == chrom_string:
                    mismatch_dict["seq_of_origin"].append(file_name.split(".")[0])
                    mismatch_dict["TRD"].append(repeat.TRD)
                    mismatch_dict["begin"].append(repeat.begin)
                    mismatch_dict["repeat_region_length"].append(repeat.repeat_region_length)
                    mismatch_dict["msa"].append(repeat.msa)
                    mismatch_dict["ref_seq"].append(chrom_string)
                    mismatch_count += 1

    df = pd.DataFrame(mismatch_dict)
    df.to_csv(path_or_buf="mismatches.tsv", sep="\t", index=False)

    print(f"Validation done, encountered -- {mismatch_count} -- mismatches between repeats and reference sequence")

if __name__ == '__main__':
    main()
