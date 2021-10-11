#!/usr/bin/env python3
# import sys
# import os
# try:
#     slurm_jobid = os.environ["SLURM_JOBID"]
#     sys.path.insert(0, os.path.join("/data/scratch", slurm_jobid, "python/lib/python3.6/site-packages/"))
# except:
#     pass
import argparse
import pickle

from tral.repeat import repeat_align
from tral.hmm.hmm import HMM

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--rlist", "-r", type=str, required=True, help="Path to repeatlist file"
    )
    parser.add_argument(
        "--realign", "-e", action='store_true'
    )
    parser.add_argument(
        "--hmm", "-m", action='store_true'
    )

    return parser.parse_args()

def main():
    args = cla_parser()

    with open(args.rlist, "rb") as f:
        repeat_list = pickle.load(f)

    print("starting")
    for i, repeat in enumerate(repeat_list.repeats):
        if args.hmm:
            for i in range(0, 9):
                hmm = HMM.create(input_format='repeat', repeat=repeat)
        if args.realign:
            new_msa = repeat_align.realign_repeat(
                repeat.msa,
                "mafft",
                "DNA",
                rate_distribution="constant"
            )
        print(f"Finished repeat number {i}")
    print("\nfinished")

if __name__ == "__main__":
    main()
