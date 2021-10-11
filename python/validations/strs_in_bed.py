#!/usr/bin/env python3

def main():
    regions_bed = "/cfs/earth/scratch/verb/projects/CRC_STRs/data/TCGA_genes/gene_utr_mate_pair_merged_6000.txt"
    str_coords = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/db/regions_consensus_lai_sun_thresh.bed"

    regions_dict = dict()
    with open(regions_bed, "r") as f:
        for line in f:
            line = line.split("\t")
            try:
                regions_dict[line[0]].append((int(line[1]), int(line[2])))
            except KeyError:
                regions_dict[line[0]] = [(int(line[1]), int(line[2]))]

    with open(str_coords, "r") as f:
        for line in f:
            line = line.split("\t")
            in_region = False
            for region in regions_dict[line[0]]:                
                if region[0] <= int(line[1]) and region[1] >= int(line[2]):
                    in_region = True
                    break
            if not in_region:
                raise Exception("STR not in gene range!!!\n{}".format(line))
        print("All STRs are within the gene regions")
    

if __name__ == "__main__":
    main()
