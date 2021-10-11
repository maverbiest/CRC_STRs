#!/usr/bin/env python3
""" Script to check whether repeats where all units have the same length can have gaps in them or not
If these types of repeats never have gaps, it may be a new criterion to filter on when deciding
whether we need to realign a repeat or not.
"""
import argparse
import os
import sys
# ugly sys hack for parent dir imports
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from scoring_filtering.tral_pvalues import load_repeatlists

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--repeat_dir", "-r", type=str, required=True, help="Directory where pickled RepeatLists are stored"
    )

    return parser.parse_args()


def main():
    args = cla_parser()

    repeat_dir = args.repeat_dir

    equal_lengths_gap = 0
    equal_lengths_no_gap = 0
    equal_lengths_longer_than_5 = 0
    total = 0
    unit_nums = dict()

    for file_name, repeatlist in load_repeatlists(repeat_dir):
        for repeat in repeatlist.repeats:
            gaps = False
            unequal_lengths = False
            first_unit_len = len(repeat.msa[0].replace("-", ""))
            for unit in repeat.msa:
                if "-" in unit:
                    gaps = True
                if len(unit.replace("-", "")) != first_unit_len:
                    unequal_lengths = True
            
            if not unequal_lengths:
                if len(repeat.msa) >= 5:                    
                    equal_lengths_longer_than_5 += 1                    
                if gaps:
                    if len(repeat.msa) >= 4:
                        print(repeat)
                    # print(repeat)
                    # print(repeat.msa)
                    try:
                        unit_nums[len(repeat.msa)] += 1
                    except KeyError:
                        unit_nums[len(repeat.msa)] = 1
                    # if not len(repeat.msa) == 2:
                    #     print(repeat)
                    #     print(repeat.msa_original)
                    equal_lengths_gap += 1
                    # if equal_lengths_gap == 100:
                    #     exit()
                else:
                    equal_lengths_no_gap += 1

            total += 1

    # test = [["AAACAA", "AAACAA", "AAACAA"], ["CGC-GA", "CGC-GA", "CGCT-A"], ["TCG", "TCG" "TCG", "TC-"]]
    # for repeat in test:
    #     gaps = False
    #     unequal_lengths = False
    #     first_unit_len = len(repeat[0].replace("-", ""))
    #     for unit in repeat:
    #         if "-" in unit:
    #             gaps = True
    #         if len(unit.replace("-", "")) != first_unit_len:
    #             unequal_lengths = True
        
    #     if not unequal_lengths:
    #         if gaps:
    #             equal_lengths_gap += 1
    #         else:
    #             equal_lengths_no_gap += 1

    #     total += 1
    
    print(f"Found a total of {total} repeats")
    print(f"{equal_lengths_gap + equal_lengths_no_gap} had equal unit lengths, {equal_lengths_gap} of which contained at least one gap in their msa")
    print(f"Equal length and len(msa) >= 5: {equal_lengths_longer_than_5}")
    print(unit_nums)


if __name__ == "__main__":
    main()
