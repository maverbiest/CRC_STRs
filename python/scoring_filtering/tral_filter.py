#!/usr/bin/env python3
"""
Collect repeatlists from a directory and filter the associated repeats based on command line arguments.
The fitlered repeatlists will then be serialized to a new output location.
"""

import argparse
import os

from tral.repeat.repeat import Repeat
from tral.repeat_list.repeat_list import RepeatList

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
        "--model", "-m", type=str, required=True, help="Model to use for filtering repeats. Options: phylo, phylo_gap01, phylo_gap001"
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

def filter_and_correct_tails(repeat_list, model, pval, div):
    filtered_list = []
    for repeat in repeat_list.repeats:
        if repeat.d_pvalue[model] >= pval or repeat.d_divergence[model] >= div:
            trailing_unit = repeat.msa[-1].replace("-", "")
            if len(trailing_unit) <= 0.5 * repeat.l_effective and len(repeat.msa) >= 3:
                # print(repeat, repeat.repeat_region_length)
                trunc_repeat = Repeat(
                    repeat.msa[0:-1],
                    begin=repeat.begin,
                    sequence_type="DNA",
                    scoreslist=[model], 
                    calc_score=True, 
                    calc_pvalue=True
                    )
                # print(trunc_repeat, trunc_repeat.repeat_region_length)
                if trunc_repeat.d_pvalue[model] < pval and trunc_repeat.d_divergence[model] < div:
                    filtered_list.append(trunc_repeat)
            continue
        filtered_list.append(repeat)
    return filtered_list

def main():
    args = cla_parser()

    input_dir = args.input 
    output_dir = args.output
    # input_dir = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/repeats/localscratch"
    # output_dir = "/tmp"    
    for file_name, repeat_list in load_repeatlists(input_dir):
        repeat_list_filt = RepeatList(filter_and_correct_tails(repeat_list, model="phylo_gap01", pval=0.05, div=0.01))
        
        # filter out repeats that do not pass thresholds for pvalue and divergence
        # repeat_list_filt = repeat_list.filter(
        #     "pvalue", 
        #     args.model,
        #     args.pvalue)

        # repeat_list_filt = repeat_list_filt.filter(
        #     "divergence",
        #     args.model,
        #     args.divergence)
    
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
        criterion_list = [("pvalue", args.model), ("divergence", args.model)]
        repeat_list_clust = repeat_list_filt.filter("none_overlapping", ["common_ancestry"], criterion_list)

        new_file_name = file_name.split(".")[0] + "_filt.pickle"
        output_path = os.path.join(output_dir, new_file_name)
        repeat_list_clust.write(output_format="pickle", file=output_path)        

if __name__ == "__main__":
    main()
