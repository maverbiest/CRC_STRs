#!/usr/bin/env python3
"""
Script to detect Tandem Repeats in a specified fasta file containing protein sequences using TRAL.
Output will be generated in a specified output directory (which will be made if it does not exist).
Output will consist of one .tsv file and one .pkl file for each protein in which a TR is detected. These files can be
merged into one file containing all Tandem Repeats by running the separate 'merge_tral_results.py' script.

NOTE: this script is heavily based on 'TR_in_multiple_Protein.py' (https://github.com/matteodelucchi/CRC_TRs)
and 'tandem_repeat_annotation_scripts.py'
    (https://github.com/acg-team/swissrepeats/tree/master/src/tandem_repeat_annotation_workflows)

Author: Max Verbiest
Contact: max.verbiest@zhaw.ch
"""

import sys
import os
import pickle
import argparse
import numpy as np
import logging
import logging.config

from Bio import SeqIO
from collections import defaultdict

from tral.paths import config_file, PACKAGE_DIRECTORY
from tral import configuration
from tral.sequence import sequence
from tral.repeat import repeat_pvalue
from tral.repeat_list import repeat_list
from tral.hmm import hmm


def find_protein_repeats(sequences_file, result_dir, refine=True):
    """Detect TRs in all protein entries in a specified .fasta file

    Parameters:
    sequences_file (str):   path to .fasta file containing sequences for TR detection
    result_dir (str):       path to directory where results will be deposited (one file will be generated
                            per protein entry in the input file)
    """
    # TODO implement checking of potential results dir for sequences to skip
    logging.config.fileConfig(config_file("logging.ini"))
    log = logging.getLogger('root')

    CONFIG_GENERAL = configuration.Configuration.instance().config
    CONFIG = CONFIG_GENERAL["repeat_list"]
    score = CONFIG["model"]    
    seq_type = CONFIG_GENERAL["sequence_type"]
    detectors = CONFIG_GENERAL["sequence"]["repeat_detection"][seq_type]

    print("""
    -----------------------------------------------
    TRAL initialized with the following parameters:
    Sequence type: {}
    Tandem repeat detectors: {}
    Tandem repeat scoring method: {}
    Refine de novo repeats using cpHMMs: {}
    -----------------------------------------------
    Parameters can be set in ~/.tral/config.ini
    Commencing run
    -----------------------------------------------
    """.format(seq_type, detectors, score, refine))

    sequences = SeqIO.parse(sequences_file, "fasta")    

    all_denovo_repeats = 0
    all_filtered_repeats = 0

    counter = 1
    for record in sequences:
        print("Started work on sequence number {}".format(counter))
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

        if not denovo_list or len(denovo_list.repeats) == 0:
            continue

        for TR in denovo_list.repeats:
            TR.model = None
            if score in {"parsimony", "pSim"}:
                TR.calculate_pvalues()
            else:
                manual_pvalue(TR, score, seq_type)
            # if score in {"parsimony", "pSim"}:
            #     TR.calculate_pvalues()
            # else:      
            #     TR.d_pvalue = defaultdict(int)
            #     empirical_list = manual_empirical_list(TR, score, seq_type)
            #     empirical_list = sorted(empirical_list)
            #     TR.d_pvalue[score] = repeat_pvalue.pvalue_from_empirical_list(
            #         tandemrepeat = TR,
            #         score_type = score,
            #         score_value = TR.d_score[score],    
            #         empirical = empirical_list
            #     )  

        seq.set_repeatlist(denovo_list, "denovo_all")

        # add number of denovo found repeats
        all_denovo_repeats += len(seq.get_repeatlist("denovo_all").repeats)

        ##########################################################################
        # Filtering TRs
        denovo_list_filtered = filter_repeatlist(seq.get_repeatlist("denovo_all"), score=score)

        if not denovo_list_filtered or len(denovo_list_filtered.repeats) == 0:
            continue

        ##########################################################################
        # Clustering
        # De novo TRs are clustered for overlap (common ancestry). Only best =
        # lowest p-Value and lowest divergence are retained.
        criterion_list = [("pvalue", score)]
        if not score in {"entropy", "parsimony", "pSim"}:
            criterion_list.append(("divergence", score))
        denovo_list_filtered = denovo_list_filtered.filter(
            "none_overlapping", ["common_ancestry"], criterion_list)
        seq.set_repeatlist(denovo_list_filtered, "denovo_filtered")

        ##########################################################################
        # Optionally, de novo TRs are refined with HMMs
        if refine:
            final_list = refine_repeatlist(seq, "denovo_filtered", score, seq_type)
            seq.set_repeatlist(final_list, "denovo_final")
        else:
            seq.set_repeatlist(denovo_list_filtered, "denovo_final")

        if len(seq.get_repeatlist("denovo_final").repeats) == 0:
            continue

        ##########################################################################
        # Save Tandem Repeats TODO: implement pickle binary dump
        # Create output directory if not already exists.
        try:
            if not os.path.isdir(result_dir):
                os.makedirs(result_dir)
        except:
            raise Exception(
                "Could not create path to result directory: {}".format(
                    os.path.dirname(result_dir)))

        # create filename
        output_tsv_file = os.path.join(result_dir, seq_name + ".tsv")

        # save TR-file in tsv format 
        write_file(seq.get_repeatlist("denovo_final").repeats, output_tsv_file, score=score)

        all_filtered_repeats += len(seq.get_repeatlist("denovo_final").repeats)
        print("\n***", seq_name, "***")
        print("denovo repeats:", len(denovo_list.repeats))
        print("repeats after filtering and clustering:",
              len(seq.get_repeatlist("denovo_final").repeats))

        for i in range(len(seq.get_repeatlist("denovo_final").repeats)):
            print(seq.get_repeatlist("denovo_final").repeats[i])

    return print("\nThere where {} repeats found de novo.".format(all_denovo_repeats),
                 "After filtering and clustering there where only {} repeats left.\n".format(
                     all_filtered_repeats))


def filter_repeatlist(tr_list, score, pvalue_threshold=0.05, divergence_threshold=0.3, n_threshold=2.5, score_threshold=0.3):
    """Filter a list of tandem repeats based on specified criteria

    Parameters
    tr_list (tral.repeat_list.RepeatList):
                            List of repeats to filter
    score (str):            Score model used to evaluate TRs
    pvalue_threshold (float):
                            P value cutoff for filter
    divergence_threshold (float):
                            Divergence cutoff for filter
    n_threshold (float):    Minimum number of repeat unit repetitions for repeat to be kept

    Returns
    tr_list_filtered (RepeatList):
                            Filtered list of repeats
    """
    # filtering for pvalue
    tr_list_filtered = tr_list.filter(
        "pvalue",
        score,
        pvalue_threshold)

    if not score in {"entropy", "parsimony", "pSim"}:
        # filtering for divergence, only works with model-based TR scoring (?)
        tr_list_filtered = tr_list_filtered.filter(
            "divergence",
            score,
            divergence_threshold)
    elif not score_threshold is None:     
        tr_list_filtered = repeat_list.RepeatList([tr for tr in tr_list_filtered.repeats if tr.d_score[score] <= score_threshold])   
    
    # filtering for number of repeat units
    tr_list_filtered = tr_list_filtered.filter(
        "attribute",
        "n_effective",
        "min",
        n_threshold)

    return tr_list_filtered


def refine_repeatlist(seq, tag, score, seq_type):
    """Take a sequence with one or more repeat lists associated to it, and a tag specifying repeat list of interest
    (one sequence can have multiple repeat lists associated to it, each identifiable with a tag). Then, create a cpHMM
    for each tandem repeat, detect TRs in sequence again and see if this improves (refines) the de novo TR

    seq (tral.sequence.Sequence):
                    A sequence with one or more RepeatLists associated to it    
    tag (str):      Which RepeatList associated to seq should be used?
    score (str):    Score used to evaluate TRs

    Returns
    refined_list (list):
                    A list containing HMM refined tandem repeats
    """

    refined_list = []
    for TR in seq.get_repeatlist(tag).repeats:
        use_refined = False
        denovo_hmm = hmm.HMM.create(input_format='repeat', repeat=TR)
        # Run HMM on sequence
        denovo_refined_list = seq.detect(lHMM=[denovo_hmm])
        if denovo_refined_list and denovo_refined_list.repeats:
            TR_refined = denovo_refined_list.repeats[0]
            TR_refined.TRD = TR.TRD
            TR_refined.model = "cpHMM"
            TR_refined.calculate_scores(scoreslist=[score])
            manual_pvalue(TR_refined, score, seq_type)
            # Check whether new and old TR overlap. Check whether new TR is
            # significant. If not both, put unrefined TR into final.
            if repeat_list.two_repeats_overlap(
                    "shared_char",
                    TR,
                    TR_refined):
                tmp_repeatlist = repeat_list.RepeatList([TR_refined])
                tmp_repeatlist_filtered = filter_repeatlist(tmp_repeatlist, score=score)
                if tmp_repeatlist_filtered.repeats:
                    use_refined = True
        if use_refined:
            refined_list.append(TR_refined)
        else:
            refined_list.append(TR)
    return repeat_list.RepeatList(refined_list)


def manual_empirical_list(repeat, score_type, sequence_type):
    """
    Copied and slightly altered from tral.repeat.repeat_pvalue to handle the fact that model based 
    null-distribution files for DNA sequences look different than those for AA sequences in the downloadable 
    TRAL pvalue data files.
    The package function tries to access the column of scores in the null-distribution file through index ["arr_0"],
    which is not present for DNA files so there will be a KeyError. This current custom method instead looks for
    the column 'value' in these files, however I am uncertain these values are valid.
    """
    l_max = 99
    n_max = 49
    if score_type == 'entropy':
        l_max = 20
        n_max = 5

    total_repeat_length_max = 1000
    l = min(repeat.l_effective, l_max)
    n = min(repeat.n, n_max)
    if l * n > total_repeat_length_max:
        # LOG.debug(
        #     "l: %d and n: %d are bigger than the total_repeat_length_max: %d" %
        #     (l, n, total_repeat_length_max))
        # Dirty little hack: as we do not have data for this pair of l and n,
        # apply nearest data file.
        while l * n > total_repeat_length_max:
            if l > n:
                l -= 1
            else:
                n -= 1

    try:
        file = config_file("data",
                           "pvalue",
                           sequence_type,
                           score_type,
                           str(l) + '_' + str(n) + '.npz')
    except ValueError:
        filename = os.path.join(CONFIG_DIR,
                                "data",
                                "pvalue",
                                sequence_type,
                                score_type,
                                str(l) + '_' + str(n) + '.npz')
        raise ValueError("Complete pdf file %s to calculate the pvalue does not exist!" % filename)

    empirical_list_all = np.load(file)
    # It is necessary to close the numpy filehandle.
    # Otherwise, there will be too many open operating system filehandles
    # if you run this function many times (more than ulimit -n allows, that is)
    if sequence_type == "AA":
        empirical_list_part = empirical_list_all['arr_0']
    elif sequence_type == "DNA":
        empirical_list_part = empirical_list_all['value']
    empirical_list_all.close()

    return empirical_list_part

def manual_pvalue(repeat, score, seq_type):
    """
    Custom pvalue calculation to avoid calling 'empirical_list' from base TRAL, which will cause errors.
    """
    if score in {"parsimony", "pSim"}:
        repeat.calculate_pvalues()
    else:      
        repeat.d_pvalue = defaultdict(int)
        empirical_list = manual_empirical_list(repeat, score, seq_type)
        empirical_list = sorted(empirical_list)
        repeat.d_pvalue[score] = repeat_pvalue.pvalue_from_empirical_list(
            tandemrepeat = repeat,
            score_type = score,
            score_value = repeat.d_score[score],    
            empirical = empirical_list
        )


def write_file(tr_list, destination, score):
    """Write all tandem repeats from a list of Repeats to a specified output file (.tsv format)

    Parameters
    tr_list (list):
                    A basic Python list (NOT RepeatList) of Repeat instances to be written to file
    destination (str):
                    Name of file to be made
    """
    # Sort TRs based on order of appearance in the input sequence
    sorted_trs = sorted(tr_list, key=lambda x: x.begin)

    header = "\t".join(["begin",
                        "msa_original",
                        "l_effective",
                        "n_effective",
                        "repeat_region_length",
                        "score",
                        "pvalue"])
    if not score in {"entropy", "parsimony", "pSim"}:
        header += "\tdivergence"
    with open(destination, "w") as f:
        f.write(header)
        for tr in sorted_trs:
            line = [
                str(i) for i in [
                    tr.begin,
                    ",".join(tr.msa),
                    tr.l_effective,
                    tr.n_effective,
                    tr.repeat_region_length,
                    tr.d_score[score],
                    tr.d_pvalue[score]
                ]
            ]
            if not score in {"entropy", "parsimony", "pSim"}:
                line.append(str(tr.d_divergence[score]))
            f.write("\n" + "\t".join(line))



def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta", "-f", type=str, required=True, help="Path to file containing protein sequence(s)"
    )
    parser.add_argument(
        "--outdir", "-o", type=str, required=True, help="Path to output directory"
    )
    parser.add_argument(
        "--skip_refine", "-s", action="store_false", help="Flag to control refining of de novo TRs using cpHMMs"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parser()

    find_protein_repeats(args.fasta, args.outdir, args.skip_refine)
