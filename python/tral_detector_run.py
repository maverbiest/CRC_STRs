#!/usr/bin/env python3

import os
import argparse
import logging
import logging.config

from Bio import SeqIO

from tral.paths import config_file, PACKAGE_DIRECTORY
from tral import configuration
from tral.sequence import sequence


def detect_trs(fasta_handle, output_dir, detectors=None, write=True):
    logging.config.fileConfig(config_file("logging.ini"))
    log = logging.getLogger('root')

    CONFIG_GENERAL = configuration.Configuration.instance().config
    CONFIG = CONFIG_GENERAL["repeat_list"]
    seq_type = CONFIG_GENERAL["sequence_type"]
    if not detectors:
        detectors = CONFIG_GENERAL["sequence"]["repeat_detection"][seq_type]

    print("""
    -----------------------------------------------
    TRAL initialized with the following parameters:
    Sequence type: {}
    Tandem repeat detectors: 
        {}
    -----------------------------------------------
    Parameters can be set in ~/.tral/config.ini
    Commencing run
    -----------------------------------------------
    """.format(seq_type, detectors))

    finished_sequences = check_output_dir(output_dir)
    sequences = SeqIO.parse(fasta_handle, "fasta")

    total_tr_counter = 0
    seq_counter = 0
    for record in sequences:
        seq_name = record.id

        # Check whether output for this genomic region already exists. If so, skip
        if seq_name in finished_sequences:
            continue
        
        seq_counter += 1
        print("Started work on sequence number {} ({})".format(seq_counter, seq_name))
        

        # Use sequence reformatted header from fasta files as name
        seq = sequence.Sequence(seq=str(record.seq), name=seq_name)       

        denovo_list = seq.detect(
            denovo=True, 
            sequence_type=seq_type,
            detection = {"detectors": detectors}     
            )

        tr_count = len(denovo_list.repeats)
        print("------ Detected {} repeat(s) in sequence\n".format(tr_count))

        total_tr_counter += tr_count

        seq.set_repeatlist(denovo_list, tag="denovo")

        if write:
            output_file_name = os.path.join(output_dir, seq_name + ".pickle")
            seq.get_repeatlist(tag="denovo").write(output_format="pickle", file=output_file_name)
        else:
            sorted_trs = sorted(seq.get_repeatlist(tag="denovo").repeats, key=lambda x: x.begin)
            for tr in sorted_trs:
                print(tr)
    print("""
    ------------------------------------------------
    Run finished!
    Detected a total of {} TRs across {} sequence(s)
    ------------------------------------------------
    """.format(total_tr_counter, seq_counter))


def check_output_dir(output_dir):
    """
    Checks the output directory for files generated in previous runs, these can be skipped later by detect_trs()
    Checking is done quite naively, only looking for files ending in '.pickle' (so no support for .pcl, .pkl ...)
    Parameters:
    output_dir (str):   Directory to check for output from previous runs

    Retruns:
    finished_sequences (set):   Set of genomic regions that can be skipped by detect_trs()
    """
    finished_sequences = {i.replace(".pickle", "") for i in os.listdir(output_dir) if i.endswith(".pickle")}
    return finished_sequences


def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta_file", "-f", type=str, required=True, help="Path to input fasta file"
    )
    parser.add_argument(
        "--output_dir", "-o", type=str, required=True, help="Path to directory where output TRs will be deposited"
    )

    return parser.parse_args()


def main():
    try:
        args = cla_parser()
        detect_trs(args.fasta_file, args.output_dir)
    except:
        test_path = "/cfs/earth/scratch/verb/projects/CRC_STRs/data/test/refining/seqs/tiny.fa"
        # test_path = "/cfs/earth/scratch/verb/projects/CRC_STRs/data/test/avg_protein_coding_gene.fa"
        test_output = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/"
        detectors = ["TRF", "XSTREAM", "Phobos"]
        # detectors = ["Phobos"]
        detect_trs(test_path, test_output, write=False, detectors=detectors)


if __name__ == "__main__":
    main()
