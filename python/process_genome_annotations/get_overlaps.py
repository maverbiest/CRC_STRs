#!/usr/bin/env python3

import argparse
import copy
import json

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f", "--filenames", type=str, required=True, help="File with filenames of genomic regions to detect overlap in. One file name per line"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="File where the info on merged regions will be written (json format)"
    )

    return parser.parse_args()

def main():
    args = cla_parser()
    regions = args.filenames

    # initialize required variables
    overlap_dict = {f"chr{i}": [] for i in [*range(1, 23), "M", "X", "Y"]}
    current_chrom = ""
    current_overlap = {"starts": [], "ends": [], "files": [], "range": None}

    with open(regions, 'r') as f:        
        for line in f:
            # filenames will look like './part_13/chr1_64091-70008_fw_refined.pickle'
            file_name = line.strip()
            genome_region = line.strip().split("/")[-1].split("_")[:2]
            genome_region = [genome_region[0], *genome_region[1].split("-")]
            if not current_chrom == genome_region[0]:
                # we arrived at a new chromosome
                # check if we are in a merged region, if so: add to overlap_dict
                if len(current_overlap["files"]) > 1:
                    current_overlap["range"] = (min(current_overlap["starts"]), max(current_overlap["ends"]))
                    overlap_dict[current_chrom].append(copy.deepcopy(current_overlap))
                # reset current chromosome and overlap values
                current_chrom = genome_region[0]
                current_overlap["starts"] = [int(genome_region[1])]
                current_overlap["ends"] = [int(genome_region[2])]
                current_overlap["files"] = [file_name]
                current_overlap["range"] = None
                continue
            # check if the current entry (line) overlaps with the previous region
            if int(genome_region[1]) <= max(current_overlap["ends"]):
                # overlap: append to currently growing overlap region
                current_overlap["starts"].append(int(genome_region[1]))
                current_overlap["ends"].append(int(genome_region[2]))
                current_overlap["files"].append(file_name)
            else:
                # no overlap
                if len(current_overlap["files"]) > 1:
                    # if we are currently in a merged region (>1 file), set range and add to overlap_dict
                    current_overlap["range"] = (min(current_overlap["starts"]), max(current_overlap["ends"]))
                    overlap_dict[current_chrom].append(copy.deepcopy(current_overlap))
                # reset current overlap values
                current_overlap["starts"] = [int(genome_region[1])]
                current_overlap["ends"] = [int(genome_region[2])]
                current_overlap["files"] = [file_name]
                current_overlap["range"] = None

    # write to json
    with open(args.output, "w") as o:
        json.dump(overlap_dict, o)

if __name__ == "__main__":
    main()
