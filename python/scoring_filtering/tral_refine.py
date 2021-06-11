#!/usr/bin/env python3
import argparse
import os

from tral.sequence.sequence import Sequence
from tral.repeat.repeat import Repeat
import tral.repeat.repeat_align as repeat_realign
from tral.hmm.hmm import HMM
from tral.hmm import hmm_viterbi
from Bio import SeqIO
import numpy as np

from tral_pvalues import load_repeatlists

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--repeat_dir", "-r", type=str, required=True, help="Directory where pickled RepeatLists are stored"
    )
    parser.add_argument(
        "--fasta_file", "-f", type=str, required=True, help="Fasta file containing the sequence that the TRs originate from. Fasta file is expected to contain multiple sequences, out of which the one matching the repeatlist will be selected"
    )
    parser.add_argument(
        "--score_type", "-s", type=str, required=True, help="Type of score to use"
    )
    

    return parser.parse_args()

def refine_repeatlist(seq, tag, score, seq_type):
    """ Take a sequence with one or more repeat lists associated to it, and a tag specifying repeat list of interest
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
    skipped_count, realigned_count = 0, 0

    print(f"Length of sequence {len(seq.seq)}")
    print(f"Total number of repeats in RepeatList: {len(seq.get_repeatlist(tag).repeats)}")
    for TR in seq.get_repeatlist(tag).repeats:
        use_refined = False
        denovo_hmm = HMM.create(input_format='repeat', repeat=TR)
        most_likely_path = denovo_hmm.viterbi(seq.seq)
        if not most_likely_path:
            continue

        unaligned_msa = hmm_viterbi.hmm_path_to_non_aligned_tandem_repeat_units(
                    seq.seq,
                    most_likely_path,
                    denovo_hmm.l_effective)

        if new_units_same_as_old(TR, unaligned_msa):
            skipped_count += 1            
            print(TR)
            print("Use old!")
            print(unaligned_msa)
            continue

        # skip_realign, checked_msa = skip_realign_current_check(unaligned_msa)
        skip_realign, checked_msa = skip_realign_new_check(unaligned_msa)        
        if not skip_realign:
            realigned_count += 1
            # checked_msa = repeat_realign.realign_repeat(checked_msa,
            #                                             "mafft",
            #                                             "DNA")
        else:
            print(TR)
            print(checked_msa)
            skipped_count += 1        
        
    print(f"skipped: {skipped_count}, realigned: {realigned_count}")

def skip_realign_current_check(msa):
    length_tuple = (len(unit) for unit in msa)
    msa = fill_out_msa_units(msa, max(length_tuple))
    
    if len(msa) <= 2 and len(msa[0]) <= 10:
        return True, msa
    return False, msa

def skip_realign_new_check(msa):
    length_tuple = tuple((len(unit) for unit in msa))

    if check_if_homogeneous(msa):
        return True, fill_out_msa_units(msa, max(length_tuple))
    
    msa = fill_out_msa_units(msa, max(length_tuple))

    # existing check in TRAL
    if len(msa) <= 2 and len(msa[0]) <= 10:
        return True, msa
    
    # are all units of length 1?
    if len(msa[0]) == 1:
        return True, msa

    # are all units of the same length?
    if len(set(length_tuple)) == 1:
        return True, msa

    return False, msa

def fill_out_msa_units(repeat_units, max_unit_length):
    pre_processed_units = []
    for unit in repeat_units:
        filled_unit = f"{unit}{'-' * (max_unit_length - len(unit))}"
        pre_processed_units.append(filled_unit)
    return pre_processed_units

def new_units_same_as_old(unrefined_repeat, hmm_msa):
    # in order for the new and old msa to match, they have to contain an equal number of units
    if not len(unrefined_repeat.msa) == len(hmm_msa):
        return False
    
    # do the nucleotides/amino acids in each unit (excluding gaps) of the old msa match the new msa? 
    for i, unit in enumerate(unrefined_repeat.msa):
        if not unit.replace("-", "") == hmm_msa[i]:
            return False

    # all units match, therefore it is safe to use the old repeat and skip further refinement/realignment
    return True

def check_if_homogeneous(msa):
    ### WORK IN PROGRESS ###
    array = np.array([list(unit) for unit in msa], dtype=object).transpose()
    array = array.transpose()
    unique_chars_per_col = tuple((len(set("".join(col))) for col in array))
    # unique_chars_per_col = tuple((len(set("".join(char))) for char in zip(*msa)))
    return all(i == 1 for i in unique_chars_per_col)


def main():
    # args = cla_parser()

    # repeat_dir = args.repeat_dir
    # fasta_file = args.fasta_file
    # score_type = args.score_type

    repeat_dir = "/cfs/earth/scratch/verb/projects/CRC_STRs/data/test/refining/repeats/"
    fasta_file = "/cfs/earth/scratch/verb/projects/CRC_STRs/data/test/refining/seqs/tiny.fa"
    score_type = "phylo_gap01"

    for file_name, repeatlist in load_repeatlists(repeat_dir):
        seq_of_origin = file_name.replace(".pickle", "")
        if "_filt" in seq_of_origin:
            seq_of_origin = seq_of_origin.replace("_filt", "")

        seq_parser = SeqIO.parse(fasta_file, "fasta")
        for record in seq_parser:
            if record.id == seq_of_origin:
                seq = Sequence(seq=str(record.seq), name=seq_of_origin, sequence_type="DNA")
                break
        
        seq.set_repeatlist(repeatlist, "unrefined")
        refined_list = refine_repeatlist(
            seq, 
            "unrefined",
            score_type, 
            "DNA"
            )

if __name__ == "__main__":
    main()
