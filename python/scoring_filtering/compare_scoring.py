#!/usr/bin/env python3
import argparse
import pickle

import pandas as pd

from tral.repeat_list.repeat_list import RepeatList
from tral_pvalues import load_repeatlists


def compare_divergence_cutoffs(source_seq, repeat_list, score_type, cutoff_one, cutoff_two):
    repeat_dict = {
        "sequence": [],
        "TRD": [],
        "begin": [],
        "l_effective": [],
        "msa": [],
        "pvalue": [],
        "divergence": []
    }
    for repeat in repeat_list.repeats:
        if repeat.repeat_region_length > 15:
            continue
        if repeat.d_pvalue[score_type] < 0.05 and repeat.d_divergence[score_type] < cutoff_two:
            if repeat.d_divergence[score_type] >= cutoff_one:
                repeat_dict["sequence"].append(source_seq)
                repeat_dict["TRD"].append(repeat.TRD)
                repeat_dict["begin"].append(repeat.begin)
                repeat_dict["l_effective"].append(repeat.l_effective)
                repeat_dict["msa"].append(" ".join(repeat.msa_original))
                repeat_dict["pvalue"].append(repeat.d_pvalue[score_type])
                repeat_dict["divergence"].append(repeat.d_divergence[score_type])
    
    return pd.DataFrame.from_dict(repeat_dict)
    
        
def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--repeat_dir", "-r", type=str, required=True, help="Directory where pickled RepeatLists are stored"
    )
    parser.add_argument(
        "--score_type", "-s", type=str, required=True, help="Score type to compare divergence cutoffs for"
    )
    parser.add_argument(
        "--score_one", "-o", type=float, required=True, help="The first (lowest) cutoff value"
    )    
    parser.add_argument(
        "--score_two", "-t", type=float, required=True, help="The second (higher) cutoff value to compare to the first"
    )
    parser.add_argument(
        "--output", "-u", type=str, default=None, help="(optional) path where difference between cutoffs will be written"
    )

    return parser.parse_args()


def main():
    args = cla_parser()

    # Double check to see if collected objects are RepeatList
    difference_df = pd.DataFrame()
    for file_name, repeatlist in load_repeatlists(args.repeat_dir):
        source_sequence = file_name.replace(".pickle", "")
        current_diff_df = compare_divergence_cutoffs(
                source_sequence, 
                repeatlist, 
                args.score_type, 
                args.score_one, 
                args.score_two
                )
        if difference_df.empty:
            difference_df = current_diff_df
        else:
            difference_df = difference_df.append(current_diff_df, ignore_index=True)
    
    difference_df.sort_values(by=["l_effective", "begin"], ascending=[True, True], inplace=True)

    if args.output:
        difference_df.to_csv(args.output, sep="\t", index=False)
    else:
        print(difference_df.head())
        print(difference_df.shape)
    

if __name__ == "__main__":
    main()
