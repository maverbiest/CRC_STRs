#!/usr/bin/env python3

import argparse
import pickle
import time

# from tral_detector_run import check_output_dir

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--directory", "-d", type=str, required=True, help="Directory containing pickled RepeatLists to score"
    )
    parser.add_argument(
        "--output", "-o", type=str, required=True, help="Directory where scored RepeatLists will be stored"
    )

    return parser.parse_args()


def main():
    # args = cla_parser()    
    test_file = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/repeats/chr1_1275436-1292029_fw.pickle"
    test_file_longer = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/repeats/chr1_31291982-31364953_fw.pickle"

    # to_skip = check_output_dir(args.output)
    with open(test_file_longer, "rb") as f:
        test_rlist = pickle.load(f)
    
    num_repeats = len(test_rlist.repeats)    
    filt_count_1 = 0
    filt_count_2 = 0
    filt_count_3 = 0

    pval_1 = 0.05
    pval_2 = 0.01
    pval_3 = 0.001
    for repeat in test_rlist.repeats:
        repeat.calculate_scores(scoreslist=["phylo"])
        repeat.calculate_pvalues(scoreslist=["phylo"])
        if repeat.d_pvalue["phylo"] <= pval_1:
            filt_count_1 += 1
        if repeat.d_pvalue["phylo"] <= pval_2:
            filt_count_2 += 1
        if repeat.d_pvalue["phylo"] <= pval_3:
            filt_count_3 += 1
    
    print("Repeats before filter: {}\nRepeats where p <= {}: {}".format(
        num_repeats, 
        [pval_1, pval_2, pval_3], 
        [filt_count_1, filt_count_2, filt_count_3]))

if __name__ == "__main__":
    main()
