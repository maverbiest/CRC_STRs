#!/usr/bin/env python3

import os
from Bio import SeqIO


def main():
    file_list = [
        "/cfs/earth/scratch/verb/projects/CRC_STRs/data/TCGA_genes/protein_coding/tcga_pc_genes_fw_5000.fa",
        "/cfs/earth/scratch/verb/projects/CRC_STRs/data/TCGA_genes/protein_coding/tcga_pc_genes_rv_5000.fa"
    ]

    for fasta in file_list:        
        addition = "fw"
        if "tcga_pc_genes_rv_5000.fa" in fasta:
            addition = "rv"
        
        fasta_parser = SeqIO.parse(fasta, "fasta")        

        output_dir  = fasta.split("/")[0:-1]
        output_file_name = "reformat_" + fasta.split("/")[-1]
        output_handle = os.path.join("/", *output_dir, output_file_name)

        with open(output_handle, "w") as f:
            for record in fasta_parser:
                header = record.id
                new_header = header.replace(":", "_")
                if addition == "rv":
                    new_header = new_header.replace("/rc", "")
                new_header = "{}_{}".format(new_header, addition)
                record.id = new_header
                record.description = new_header
                SeqIO.write(record, f, "fasta")


if __name__ == "__main__":
    main()
