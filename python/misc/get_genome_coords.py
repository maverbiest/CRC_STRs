#!/usr/bin/env python3
"""
Parse a gencode.gene file (e.g. gencode.gene.info.v22.tsv at https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)
Using this, create two files (one for forward strand, one for reverse), that contain genomic coordinates of 
interest to be extracted from a reference genome using SAM tools faidx ('samtools faidx --help' for details)
Also allow for the inclusion of a specified length of 5'UTR and/or 3'UTR
"""

import argparse
import os

import constants

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(        
        "--gene_info",  "-g", type=str, required=True, help="genocode.gene file containing locations of genes in reference genome"
    )
    parser.add_argument(
        "--output_dir",  "-o", type=str, default="./", help="directory to deposit output files (dir must exist)"
    )
    parser.add_argument(
        "--output_base", "-b", type=str, default="region_file", help="base for output file name, specific suffixes will be added based on arguments"
    )
    parser.add_argument(
        "--utr5", "-u5", type=int, default=0, help="int: How many bp of the 5'UTR should be included? (default=0)"
    )
    parser.add_argument(
        "--utr3", "-u3", type=int, default=0, help="int: How many bp of the 3'UTR should be included? (default=0)"
    )    

    return parser.parse_args()


def create_gene_coord_files(gene_info, output_dir, output_file_base, utr5, utr3):
    if not os.path.exists(output_dir):
        raise ValueError("Specified output directory does not exist")

    forward = []
    reverse = []    

    with open(gene_info, "r") as f:
        next(f) # skip header
        for line in f:
            line_split = line.strip().split("\t")
            chromosome, strand = line_split[2], line_split[5]
            begin, end = int(line_split[3]), int(line_split[4])
            if strand == "+":
                begin -= utr5
                end += utr3
                if begin < 0:
                    begin = 0
                if end > constants.CHROMOSOME_LENGTHS[chromosome]:
                    end = constants.CHROMOSOME_LENGTHS[chromosome]
                forward.append((chromosome, begin, end))
            elif strand == "-":
                begin -= utr3
                end += utr5
                if begin < 0:
                    begin = 0
                if end > constants.CHROMOSOME_LENGTHS[chromosome]:
                    end = constants.CHROMOSOME_LENGTHS[chromosome]
                reverse.append((chromosome, begin, end))
            else:
                raise ValueError("Unexpected strand identifier encountered: '{}' instead of +/-".format(strand))

    if forward:
        forward_file_handle = os.path.join(output_dir, get_outfile_name(output_file_base, "forward", utr5, utr3))
        write_list_to_file(forward_file_handle, forward)
    if reverse:
        reverse_file_handle = os.path.join(output_dir, get_outfile_name(output_file_base, "reverse", utr5, utr3)) 
        write_list_to_file(reverse_file_handle, reverse)

def get_outfile_name(base, strand, utr5, utr3):
    return "{}_{}_{}-{}.txt".format(base, strand, utr5, utr3)

def write_list_to_file(output_file, entry_list, format="bed"):
    if not format in {"bed", "samtools"}:
        raise ValueError(f"format identifier must be one of 'bed, samtools', not '{format}'")
    with open(output_file, "w") as f:
        for entry in entry_list:
            if format == "bed":
                f.write("{}\t{}\t{}\n".format(entry[0], entry[1], entry[2])) 
            elif format == "samtools":   
                f.write("{}:{}-{}\n".format(entry[0], entry[1], entry[2]))


def main():
    args = cla_parser()
    create_gene_coord_files(
        gene_info=args.gene_info, 
        output_dir=args.output_dir, 
        output_file_base=args.output_base, 
        utr5=args.utr5,
        utr3=args.utr3
    )   

if __name__ == "__main__":
    main()
