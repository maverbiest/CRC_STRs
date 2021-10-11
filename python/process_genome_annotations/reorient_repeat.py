#!/usr/bin/env python3
import os
import pickle

from Bio.Seq import Seq
from numpy.core.arrayprint import repr_format
from tral.repeat.repeat import Repeat


def main():
    to_reverse = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/reverse_repeats/part_25/chrM_14149-16569_rv_refined.pickle"
    basename = os.path.basename(to_reverse)

    with open(to_reverse, 'rb') as f:
        repeat_list = pickle.load(f)
        for repeat in repeat_list.repeats:
            seq_begin = int(basename.split("_")[1].split("-")[0])
            seq_end = int(basename.split("_")[1].split("-")[1])
            seq_len = seq_end - seq_begin + 1
            repeat_end = repeat.begin + repeat.repeat_region_length - 1
            new_repeat_begin = seq_len - repeat_end

            rc_msa = str(Seq(",".join(repeat.msa)).reverse_complement())
            rc_msa = rc_msa.split(",")

            new_repeat = Repeat(
                msa=rc_msa,
                begin=new_repeat_begin,
                sequence_type="DNA",
                # scoreslist=repeat.scoreslist,
                calc_score=False,
                calc_pvalue=False
            )
            new_repeat.d_score = repeat.d_score
            new_repeat.d_pvalue = repeat.d_pvalue
            new_repeat.d_divergence = repeat.d_divergence
            
            print(repeat.__dict__)
            print()
            print(new_repeat.__dict__)
            exit()


if __name__ == "__main__":
    main()
