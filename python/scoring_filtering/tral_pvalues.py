#!/usr/bin/env python3
"""
Calculate divergence and p-values for repeats from all pickled repeatlists in a directory.
Serialize scored repeats again.
"""

import argparse
import os
import pickle

from tral.repeat_list.repeat_list import RepeatList


def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--directory", "-d", type=str, required=True, help="Directory containing pickled RepeatLists to score"
    )
    parser.add_argument(
        "--model", "-m", type=str, required=True, help="Model to use for calculating pvalues. Options: phylo, phylo_gap01, phylo_gap001"
    )
    parser.add_argument(
        "--output", "-o", type=str, required=False, help="Directory where scored RepeatLists will be stored (if the same as input, files will be overwritten)"
    )    

    return parser.parse_args()

def load_repeatlists(directory):
    # collect all pickle files from input directory
    file_names = [i for i in os.listdir(directory) if i.endswith(".pickle")]
    
    for file_name in file_names:
        full_path = os.path.join(directory, file_name)

        with open(full_path, "rb") as f:
            repeat_list = pickle.load(f)
        
        # Double check to see if collected object is RepeatList
        if not isinstance(repeat_list, RepeatList):
            continue
        
        yield(file_name, repeat_list)


def main():
    args = cla_parser()   

    input_dir = args.directory 
    # if output dir was not specified, overwrite input file
    if not args.output:
        output_dir = args.directory
    else:
        output_dir = args.output

    if not args.model in {"phylo", "phylo_gap01", "phylo_gap001"}:
        raise ValueError("Model must be one of: {phylo, phylo_gap01, phylo_gap001} ")
    model = args.model

    for file_name, repeat_list in load_repeatlists(input_dir):
        for repeat in repeat_list.repeats:
            repeat.calculate_scores(scoreslist=[model])
            repeat.calculate_pvalues(scoreslist=[model])

        output_handle = os.path.join(output_dir, file_name)
        repeat_list.write(output_format="pickle", file=output_handle)


if __name__ == "__main__":
    main()
