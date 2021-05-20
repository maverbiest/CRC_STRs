#!/usr/bin/env python3

import argparse
import pickle

from tral.repeat_list.repeat_list import RepeatList

def compare_repeat_lists(list_one, one_type, list_two, two_type):
    list_one_sigs = dict()
    for repeat in list_one.repeats:
        if repeat.d_divergence[one_type] < 0.01 and repeat.d_pvalue[one_type] < 0.05:
            try:
                list_one_sigs[repeat.begin].append(repeat)
            except KeyError:
                list_one_sigs[repeat.begin] = [repeat]
    print(f"Significant TRs under score model {one_type}: {len(list_one_sigs)}")

    diff = []
    for repeat in list_two.repeats:
        if repeat.d_divergence[two_type] >= 0.01 or repeat.d_pvalue[two_type] >= 0.05:
            try:
                if any([repeat.repeat_region_length == i.repeat_region_length and repeat.TRD == i.TRD \
                    for i in list_one_sigs[repeat.begin]]):
                    diff.append(repeat)
            except KeyError:
                continue
    print(f"Significant TRs under score model {two_type}: {len(list_one_sigs) - len(diff)}")
    return diff
        
def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--score_one", "-s", type=str, required=True, help="The first pickled repeatlist containing scored repeats"
    )
    parser.add_argument(
        "--score_one_type", "-y", type=str, required=True, help="The second pickled repeatlist containing scored repeats, will be compared to first"
    )
    parser.add_argument(
        "--score_two", "-t", type=str, required=True, help="The second pickled repeatlist containing scored repeats, will be compared to first"
    )
    parser.add_argument(
        "--score_two_type", "-p", type=str, required=True, help="The second pickled repeatlist containing scored repeats, will be compared to first"
    )

    return parser.parse_args()


def main():
    args = cla_parser()

    with open(args.score_one, "rb") as f:
        list_one = pickle.load(f)
    
    with open(args.score_two, "rb") as f:
        list_two = pickle.load(f)
    
    # list_two = sorted(list_two.repeats, key = lambda x: x.begin)
    # for repeat in list_two:
    #     print(repeat, repeat.TRD)
        # if repeat.l_effective == 2:
        #     print(repeat, repeat.TRD)
    # for repeat in list_one.repeats:
    #     if repeat.begin == 147167:
    #         print(repeat)
    # exit()

    # Double check to see if collected objects are RepeatList
    if not isinstance(list_one, RepeatList) or not isinstance(list_two, RepeatList) :
        raise Exception("One or both paths did not point to pickled RepeatList!")
    print(f"Total number of TRs: {len(list_one.repeats)}")
    diff = compare_repeat_lists(list_one, args.score_one_type, list_two, args.score_two_type)
    
    for repeat in diff:
        print(repeat, repeat.TRD)
        # if repeat.TRD == "PHOBOS" or repeat.TRD == "TRF":
        #     print(repeat, repeat.TRD)
    

if __name__ == "__main__":
    main()
