#!/usr/bin/env python3
"""
Collect repeatlists from a directory and filter the associated repeats based on command line arguments.
The fitlered repeatlists will then be serialized to a new output location.
"""

import argparse
import os
import pickle
import time

from tral_pvalues import load_repeatlists

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input", "-i", type=str, required=True, help="Directory containing pickled RepeatLists to score"
    )
    parser.add_argument(
        "--output", "-o", type=str, required=True, help="Directory where scored RepeatLists will be stored"
    )
    parser.add_argument(
        "--pvalue", "-p", type=float, required=True, help="P value threshold (upper bound)"
    )
    parser.add_argument(
        "--divergence", "-d", type=float, required=True, help="Divergence threshold (upper bound)"
    )
    parser.add_argument(
        "--units", "-u", type=float, required=False, help="(Optional) Number of repeat units required (lower bound)"
    )

    return parser.parse_args()

def main():
    args = cla_parser()   

    input_dir = args.input 
    output_dir = args.output

    for file_name, repeat_list in load_repeatlists(input_dir):
        # filter out repeats that do not pass thresholds for pvalue and divergence
        repeat_list_filt = repeat_list.filter(
            "pvalue",
            "phylo",
            args.pvalue)

        repeat_list_filt = repeat_list_filt.filter(
            "divergence",
            "phylo",
            args.divergence)
    
        # optional: filtering for number of repeat units
        if args.units:            
            repeat_list_filt = repeat_list_filt.filter(
                "attribute",
                "n_effective",
                "min",
                args.units)

        # Clustering
        # De novo repeats are clustered for overlap (common ancestry). In case of overlap, best repeat
        # (i.e. lowest p-value and lowest divergence) is retained.
        criterion_list = [("pvalue", "phylo"), ("divergence", "phylo")]
        repeat_list_clust = repeat_list_filt.filter("none_overlapping", ["common_ancestry"], criterion_list)

        new_file_name = file_name.split(".")[0] + "_filt.pickle"
        output_path = os.path.join(output_dir, new_file_name)
        repeat_list_clust.write(output_format="pickle", file=output_path)        

if __name__ == "__main__":
    main()
