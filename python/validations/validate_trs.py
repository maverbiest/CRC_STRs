#!/usr/bin/env python3

import argparse
import os
import sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Bio import SeqIO
from Bio.Seq import Seq

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

    chromosomes = ["chr10"]

    for chromosome in chromosomes:
        fasta_parser = SeqIO.parse(args.reference, "fasta")
        for record in fasta_parser:
            if record.id == chromosome:
                chr_record = record
                break
        
        for file_name, repeat_list in load_repeatlists(args.tr_directory):
            if not file_name.split("_")[0] == chromosome:
                continue
            # chromosome = file_name.split("_")[0]
            # if not chromosome.startswith("chr"):
            #     raise Exception(f"Could not determine chromosome for repeatlist from {file_name}")

            if "_fw" in file_name:
                strand = "fw"
            elif "_rv" in file_name:
                strand = "rv"                
            else:
                raise Exception(f"Could not determine strand for repeatlist from file {file_name}")
            
            for repeat in repeat_list.repeats:   
                if not repeat.TRD == "PHOBOS":
                    continue
                if strand == "fw":
                    gene_begin = int(file_name.split("_")[1].split("-")[0])
                    chrom_begin = gene_begin + repeat.begin - 2 # only valid for fw strand genes
                else:
                    ### UNTESTED CODE ###
                    #TODO test
                    gene_begin = int(file_name.split("_")[1].split("-")[0])
                    gene_end = int(file_name.split("_")[1].split("-")[1])
                    repeat_end = repeat.begin + repeat.repeat_region_length - 1
                    chrom_begin = gene_end - repeat_end # only valid for rv strand genes
                repeat_string = "".join(repeat.msa).replace("-", "")
                if strand == "rv":
                    repeat_string = str(Seq(repeat_string).reverse_complement())
                chrom_string = chr_record.seq[chrom_begin : chrom_begin + repeat.repeat_region_length]
                print(f"{repeat_string}\n{chrom_string}")
                print(repeat_string == chrom_string)
                print()

if __name__ == '__main__':
    main()
