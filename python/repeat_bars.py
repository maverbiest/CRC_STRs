#!/usr/bin/env python3

import argparse

import pandas as pd
import matplotlib.pyplot as plt

from scoring_filtering.tral_pvalues import load_repeatlists

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--directory", "-d", type=str, required=True, help="Directory where pickled RepeatLists to visualize are stored"
    )
    parser.add_argument(
        "--img_name", "-i", type=str, required=True, help="Dersired path and file name to where (png) file will be generated"
    )

    return parser.parse_args()

def main():
    args = cla_parser()
    directory = args.directory
    img_name = args.img_name
    if not img_name.endswith(".png"):
        raise ValueError("Please specify output image as png format")
    # directory = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/repeats/selection/filtered"

    data_dict = {
        "gene": []
    }
    for i in range(1, 16):
        data_dict[i] = []

    index = 0
    for file_name, repeat_list in load_repeatlists(directory):
        # get genomic region name, initialize new row in dictionary
        gen_region = file_name.replace(".pickle", "")        
        data_dict['gene'].append(gen_region)
        for i in range(1, 16):
            data_dict[i].append(0)

        # add counts for each unit length to dictionary
        for repeat in repeat_list.repeats:
            try:
                data_dict[repeat.l_effective][index] += 1
            except KeyError:
                print(f"Found l_effective {repeat.l_effective} for genomic region {gen_region}. TRD: {repeat.TRD}")
        index += 1

    df = pd.DataFrame(data_dict)

    sums = df.drop(columns=['gene']).sum(axis = 0)
    total = sums.to_numpy().sum()

    ax = sums.plot.bar(grid=True)
    ax.text(0.9, 0.9, f'total={total}', ha="right", va="center", transform=ax.transAxes)

    plt.xlabel('Unit length')
    plt.ylabel('Counts')
    plt.grid(axis='x', alpha=0.5)
    plt.savefig(img_name, bbox_inches="tight")

if __name__ == "__main__":
    main()
