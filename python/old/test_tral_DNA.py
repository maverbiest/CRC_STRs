#!/usr/bin/env python3

import sys
import os
import pickle
import argparse
import logging
import logging.config

from Bio import SeqIO
from collections import defaultdict

from tral.paths import config_file, PACKAGE_DIRECTORY
from tral import configuration
from tral.sequence import sequence
from tral.repeat_list import repeat_list
from tral.repeat import repeat_pvalue
from tral.hmm import hmm

from run_tral_DNA import filter_repeatlist, manual_empirical_list


def main():
    logging.config.fileConfig(config_file("logging.ini"))
    log = logging.getLogger('root')
    log.setLevel("ERROR")
    log.propagate = False
    log.parent.setLevel("ERROR")

    CONFIG_GENERAL = configuration.Configuration.instance().config
    CONFIG = CONFIG_GENERAL["repeat_list"]
    score = CONFIG["model"]    
    seq_type = CONFIG_GENERAL["sequence_type"]
    detectors = CONFIG_GENERAL["sequence"]["repeat_detection"][seq_type]

    sequences = SeqIO.parse("/cfs/earth/scratch/verb/projects/CRC_STRs/data/test/tiny.fa", "fasta") 

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
            repeat={"calc_score": True}     
            )

        denovo_list = repeat_list.RepeatList([denovo_list.repeats[0]])

        if not denovo_list or len(denovo_list.repeats) == 0:
            continue

        for TR in denovo_list.repeats:
            TR.model = None
            if score in {"parsimony", "pSim"}:
                TR.calculate_pvalues()
            else:
                TR.d_pvalue = defaultdict(int)
                empirical_list = manual_empirical_list(TR, score, seq_type)
                empirical_list = sorted(empirical_list)
                TR.d_pvalue[score] = repeat_pvalue.pvalue_from_empirical_list(
                    tandemrepeat = TR,
                    score_type = score,
                    score_value = TR.d_score[score],    
                    empirical = empirical_list
                )  
                # if score_type in ['entropy']:
                #     if score_value == -1:
                #         TR.d_pvalue[score] = 1
                #     else:
                #         TR.d_pvalue[score] = bisect.bisect_right(empirical_list, score_value) / len(empirical_list)
                # else:
                #     TR.d_pvalue[score] = 1 - bisect.bisect_left(empirical_list, score_value) / len(empirical_list)

        seq.set_repeatlist(denovo_list, "denovo_all")

        # add number of denovo found repeats
        all_denovo_repeats += len(seq.get_repeatlist("denovo_all").repeats)
        sorted_trs = sorted(seq.get_repeatlist("denovo_all").repeats, key=lambda x: x.begin)
        for i in sorted_trs:
            print(i)

        ##########################################################################
        # Filtering TRs
        denovo_list_filtered = filter_repeatlist(seq.get_repeatlist("denovo_all"), score=score)

        if not denovo_list_filtered or len(denovo_list_filtered.repeats) == 0:
            continue


if __name__ == "__main__":
    main()
