#!/usr/bin/env python3

import sys
import os
import pickle
import argparse
import logging
import logging.config

from Bio import SeqIO

from tral.paths import config_file, PACKAGE_DIRECTORY
from tral import configuration
from tral.sequence import sequence
from tral.repeat_list import repeat_list
from tral.hmm import hmm

from run_tral_DNA import filter_repeatlist


def main():
    logging.config.fileConfig(config_file("logging.ini"))
    log = logging.getLogger('root')

    CONFIG_GENERAL = configuration.Configuration.instance().config
    CONFIG = CONFIG_GENERAL["repeat_list"]
    score = "phylo"
    seq_type = "AA"
    detectors = ["HHrepID", "T-REKS", "TRUST", "XSTREAM"]

    sequences = SeqIO.parse("/cfs/earth/scratch/verb/projects/CRC_STRs/data/test/O43318.fasta", "fasta") 

    all_denovo_repeats = 0
    all_filtered_repeats = 0

    counter = 1
    for record in sequences:
        counter += 1
        seq_name = record.id.replace(">", "")

        # name is sequence header
        seq = sequence.Sequence(seq=str(record.seq), name=seq_name)

        denovo_list = seq.detect(
            denovo=True, 
            sequence_type=seq_type,
            detection = {"detectors": detectors},
            repeat={"calc_pvalue": True, "calc_score": True}     
            )

        denovo_list = repeat_list.RepeatList([denovo_list.repeats[2]])

        if not denovo_list or len(denovo_list.repeats) == 0:
            continue

        for TR in denovo_list.repeats:
            TR.model = None

        seq.set_repeatlist(denovo_list, "denovo_all")

        # add number of denovo found repeats
        all_denovo_repeats += len(seq.get_repeatlist("denovo_all").repeats)

        ##########################################################################
        # Filtering TRs
        denovo_list_filtered = filter_repeatlist(seq.get_repeatlist("denovo_all"), score=score)

        if not denovo_list_filtered or len(denovo_list_filtered.repeats) == 0:
            continue


if __name__ == "__main__":
    main()
