#!/usr/bin/env python3

import argparse
import os

from Bio import SeqIO


def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta", "-f", type=str, required=True, help="Input fasta file to be split into smaller files"
    )
    parser.add_argument(
        "--n_parts", "-n", type=int, required=True, help="Number of parts to split input file into"
    )
    parser.add_argument(
        "--n_sequences", "-s", type=int, required=True, help="Number of sequences in input fasta. IMPORTANT: MAKE SURE THIS VALUE IS RIGHT FOR YOU INPUT"
    )
    parser.add_argument(
        "--output_base", "-o", type=str, required=True, help="Base dir of partition directory structure"
    )

    return parser.parse_args()


def main():
    args = cla_parser()

    file_handle = args.fasta
    n_parts = args.n_parts
    seqs_per_file = int(args.n_sequences / n_parts) + 1
    output_base = args.output_base

    
    fasta_parser = SeqIO.parse(file_handle, "fasta")

    for i in range(1, n_parts + 1):
        with open(os.path.join(output_base, "part_{}/part_{}.fa".format(i, i)), "w") as f:
            for j in range(1, seqs_per_file + 1):
                try:
                    SeqIO.write(next(fasta_parser), f, "fasta")
                except StopIteration:
                    print("finished all")
                    break                    
            else:
                continue
            break

if __name__ == "__main__":
    main()
