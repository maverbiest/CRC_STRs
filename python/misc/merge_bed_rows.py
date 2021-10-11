#!/usr/bin/env python3
import argparse
import os

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-b", "--bed", type=str, required=True, help="Sorted bed file where rows describing the same TR locus will be merged"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Name of deduplicated output file that will be generated"
    )

    return parser.parse_args()

def deduplicate(bed_file, output_file):
    cur_row, prev_row = [], []
    
    with open(bed_file, 'r') as f:
        for line in f:
            cur_row = line.split("\t")
            # print(line)
            try:
                if cur_row[1] == prev_row[1]:
                    cur_row = merge_rows(prev_row, cur_row)
                else:
                    # write prev_row to file here
                    print(prev_row)
                prev_row = cur_row
            except IndexError:
                # Should only happen on the first line in the file
                prev_row = cur_row


def merge_rows(prev_row, cur_row):
    prev_info = prev_row[7].split(":")
    prev_info[0] += f",{cur_row[7].split(':')[0]}"    
    
    prev_row[7] = ":".join(prev_info)
    if prev_row[4] != cur_row[4]:
        print("\t".join(cur_row), end="")
        print("\t".join(prev_row), end="")
        # print()
        # exit()    
    # if len(prev_row[7].split(':')[0].split(',')) > 2:
    #     print(cur_row)
    #     print(prev_row)
    #     exit()
    return prev_row


def main():
    args = cla_parser()

    try:
        deduplicate(args.bed, args.output)
    except FileNotFoundError:
        print("Specified bed file was not found")
        exit(1)

if __name__ == "__main__":
    main()
