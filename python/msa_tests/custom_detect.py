#!/usr/bin/env python3
import logging
import pickle
import re

from tral import configuration
from tral.repeat import repeat, repeat_align
from tral.repeat_list import repeat_list
from tral.hmm import hmm, hmm_viterbi
from tral.sequence import repeat_detection_run, sequence_io

CONFIG = configuration.Configuration.instance().config
REPEAT_CONFIG = CONFIG["repeat"]
LOG = logging.getLogger(__name__)

def detect(sequence, lHMM=None, denovo=None, realignment='mafft', sequence_type='AA', rate_distribution=REPEAT_CONFIG['castor_parameter']['rate_distribution'], user_path=None, **kwargs):
    """ Detects tandem repeats on ``sequence.seq`` from 2 possible sources.

    A list of ``Repeat`` instances is created for tandem repeat detections
    on the sequence from two possible sources:

    * Sequence profile hidden Markov models ``HMM``
    * de novo detection algorithms.

    Args:
        hmm (hmm.HMM): A list of ``HMM`` instances.
        denovo (bool): boolean
        realignment (str): either "mafft", "proPIP" or None
        *kwargs: Parameters fed to denovo TR prediction and/or Repeat
            instantiation. E.g. ``repeat = {"calc_score": True}``

    Returns:
        A ``RepeatList`` instance
    """
    realignment_types = ['mafft', 'proPIP', None]
    if realignment not in realignment_types:
        raise ValueError("Invalid realignment type. Expected one of: {}".format(realignment_types))

    if denovo:
        if 'detection' in kwargs:
            predicted_repeats = repeat_detection_run.run_detector(
                [sequence],
                **kwargs['detection'])[0]
        else:
            raise ValueError()
            # predicted_repeats = repeat_detection_run.run_detector([sequence])[0]

        LOG.debug("predicted_repeats: {}".format(predicted_repeats))
        repeats = []

        for jTRD, jlTR in predicted_repeats.items():
            for iTR in jlTR:
                if 'repeat' in kwargs:
                    iTR = repeat.Repeat(iTR.msa, begin=iTR.begin,
                                        **kwargs['repeat'])
                else:
                    iTR = repeat.Repeat(iTR.msa, begin=iTR.begin)

                # Consider only tandem repeats that have a repeat unit
                # predicted to be at least one character long.
                if iTR.l_effective > 0:

                    # Save l, n, MSA, TRD, scores, sequence_type, position
                    # in sequence of given type
                    iTR.TRD = jTRD

                    # Sanity check repeat and set begin coordinate for
                    # all repeats
                    if not sequence.repeat_in_sequence(iTR):
                        LOG.debug("The tandem repeat is not part of"
                                    "the sequence. Detector: %s", iTR.TRD)
                        continue

                    repeats.append(iTR)

        # Realignment
        for iTR in repeats:
            if len(iTR.msa) > 2 or len(iTR.msa[0]) > 10:  # TODO: change from 2 to 3
                # only TRs with at least 3 units with a length > 10 characters will be realigned
                if len(iTR.msa[0]) <= 2:
                    # Repeats with unit length 1 and 2 will not be realigned             
                    continue    
                iTR.msa = repeat_align.realign_repeat(iTR.msa,
                                                        realignment,
                                                        sequence_type,
                                                        rate_distribution=rate_distribution,
                                                        user_path=user_path)
        return repeat_list.RepeatList(repeats)

    else:
        raise Exception("Either require denovo detection, or provide an",
                        "HMM")
