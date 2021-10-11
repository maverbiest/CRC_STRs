#!/usr/bin/env python3
import os
import argparse
from typing import Iterable, Optional
from collections import Counter

import numpy as np
from Bio.Seq import Seq

from setup_db import Gene, Repeat
from gtf_to_sqlite import connection_setup

class BedMaker(object):        

    def __init__(self, 
                db_session, 
                chromosomes: Iterable[str]={"auto"}, 
                consensus_only: bool=False, 
                thresholds: Optional[dict]=None) -> None:
        self.db_session = db_session
        self.set_chroms(chromosomes)    
        self.consensus_only = consensus_only
        if not thresholds:
            ## Default threshold values taken from Lai & Sun, 2009
            self.thresholds = {
                    1: 9, # unit length: minimal number of units
                    2: 4,
                    3: 4,
                    4: 4,
                    5: 4,
                    6: 4
                }
        else:
            self.set_thresholds(thresholds)        
        self.gene_selection: Optional[dict] = None  # {gene.id: gene.strand}

    def set_thresholds(self, thresh_dict: dict) -> None:
        # Check if all keys and values in supplied dictionary are integers
        if not all(isinstance(i, int) for i in [*thresh_dict.keys(), *thresh_dict.values()]):
            raise ValueError("Thresholds must be specified as dict(int(unit_len): int(min_units))")
        self.thresholds = thresh_dict

    def set_chroms(self, chromosomes: Iterable[str]) -> None:
        """ Generate the set of chromsome identifiers from which repeats will be included

        Parameters
        chromosomes (Iterable[str]):    Iterable specifying which chromosomes to include in the 
                                        generated bed file. Allowed values: 'auto', 'sex', 'mito'
        """

        if len(chromosomes) == 0 or not all(chrom in {"auto", "sex", "mito"} for chrom in chromosomes):
            raise ValueError("Specify desired target chromosomes as (subset of): {'auto', 'sex', 'mito'}")

        targets = []
        if "auto" in chromosomes:
            targets += [f"chr{i}" for i in range(1, 23)]
        if "sex" in chromosomes:
            targets += ["chrX", "chrY"]
        if "mito" in chromosomes:
            targets += ["chrM"]
        
        self.chromosomes = targets

    def set_gene_selection(self) -> None:
        """ For each gene in the DB originating from a chromosome in 
        self.chromosomes, extract the gene.id and gene.strand to save in 
        self.gene_selection as key-value pairs {gene.id: gene.strand ... } 
        """

        self.gene_selection = dict()
        for gene in self.db_session.query(Gene).all():
            if not gene.chromosome in self.chromosomes:
                continue
            self.gene_selection[gene.id] = gene

    def threshold_filter(self):
        """ Generator function that yields repeats containing a sequence of consecutive 
        units longer than the predefined thresholds in self.thresholds
        (see self.longest_consensus_stretch() for details on consensus unit etc.)

        Returns
        repeat.gene_id (str)    ID that the gene the repeat originates from has in the DB
        longest_stretch (dict)  A dictionary containing information and genomic coordinates 
                                of the longest stretch of consensus units found in the repeat region                     
        """

        if not self.gene_selection:
            raise AttributeError("self.gene_selection must be set before generating bed file")
        
        # for current_gene_id in self.gene_selection.keys():
        for gene in self.gene_selection.values():
            # for repeat in self.db_session.query(Repeat).filter_by(gene_id=current_gene_id):
            for repeat in self.db_session.query(Repeat).filter_by(gene_id=gene.id):
                try: 
                    bed_tr = self.get_bed_tr(repeat)
                    if len(bed_tr.units) >= self.thresholds[len(bed_tr.consensus_unit)]:
                        yield bed_tr
                    # longest_stretch = self.get_bed_tr(repeat)
                    # if longest_stretch["n"] >= self.thresholds[repeat.l_effective]:
                        # yield repeat.gene_id, longest_stretch
                except KeyError:
                    # repeat.l_effective is not part of the self.threshold range currently active
                    continue

    def send_to_bed(self, output_file: str) -> None:
        """ Generate an output file in bed format that can be used by
        GangSTR to genotype STRs in alignment data
        """

        with open(output_file, "w") as o:
            for bed_tr in self.threshold_filter():
                    bed_string = bed_tr.get_bed_line()
                    o.write(bed_string)

    def get_consensus_unit(self, msa: str) -> str:    
        """ Given a comma-separated set of units representing a multiple sequence
        alignment (from DB), determine the consensus unit at each column of the alignment and 
        return the consensus unit. Insertion columns where the number of 
        insertions >= number of nucleotides are skipped. For each column, most common nt
        is selected as consensus, in case of tie: pick one at random.

        Parameters
        msa (str):  String representation of a multiple sequence alignment, with units
                    delimited by commas
        
        Returns  
        consensus_unit (str): 
                    A consensus unit derived from the provided msa
        """    
        # Convert msa into list of lists, convert into np.Array() and transpose
        msa_matrix_t = np.array([list(unit) for unit in msa.split(",")]).transpose()
        consensus_unit = ""

        for msa_col in msa_matrix_t:
            # if half or more of the msa column are gap entries, skip column
            if np.count_nonzero(msa_col == '-') >= 0.5 * len(msa_col):
                continue

            # discard gap entries, get most common nt
            msa_col = msa_col[msa_col != '-']
            consensus_unit += Counter(msa_col).most_common(1)[0][0] # most_common() will return e.g. [('A': 5)]

        return consensus_unit

    def get_bed_tr(self, db_repeat):
        """ Determine the longest consecutive sequence of consensus units in the
        repeat region. Also checks if the repeat comes from the reverse strand of the genome.
        If so, the reverse complement of the db_repeat.msa is generated.

        Parameters
        db_repeat:       A repeat region as represented in the database.

        Returns
        max_stretch (BedTR):
                        BedTR instance describing the longest sequence of units with the same length as
                        the consensus unit (though not necessarily perfect matches) found in the repeat region. 
        """

        msa = db_repeat.msa
        if self.gene_selection[db_repeat.gene_id].strand == "rv":
            # if repeat is from reverse strand, get the reverse complement of the msa
            msa = str(Seq(msa).reverse_complement())
        consensus_unit = self.get_consensus_unit(msa)

        # initialize BedTR instances to keep track of longest stretch of repeat units
        current_bed_tr = BedTR(self.gene_selection[db_repeat.gene_id].chromosome, consensus_unit, db_repeat.id)
        max_bed_tr = BedTR(self.gene_selection[db_repeat.gene_id].chromosome, consensus_unit, db_repeat.id)

        current_pos = db_repeat.begin
        
        for unit in msa.split(","):
            current_unit = unit.replace("-", "")
            if len(current_unit) == len(consensus_unit):    # current unit must be same length as consensus to be considered
                # if we only allow perfect units, check if current unit is a match. If not: reset
                if self.consensus_only and not current_unit == consensus_unit:
                    current_bed_tr = BedTR(self.gene_selection[db_repeat.gene_id].chromosome, consensus_unit, db_repeat.id)
                    current_pos += len(current_unit)
                    continue
                # code below gets a bit cumbersome to prevent shallow copies of BedTR.units lists TODO: rework using copy.deepcopy()
                # if starting a new current_bed_tr, set the begin position to current position, else just append new unit
                if not current_bed_tr.begin:
                    current_bed_tr.begin = current_pos
                    current_bed_tr.units.append(current_unit)
                else:
                    current_bed_tr.units.append(current_unit)

                # if current stretch is longer than the longest observed so far, update the max_bed_tr
                # TODO: add __lt__, __gt__ etc. for BedTR class?
                if len(current_bed_tr.units) > len(max_bed_tr.units):
                    max_bed_tr.units.append(current_bed_tr.units[-1])
                    max_bed_tr.begin = current_bed_tr.begin
            else:
                # current unit was rejected, reset the current_bed_tr
                current_bed_tr = BedTR(self.gene_selection[db_repeat.gene_id].chromosome, consensus_unit, db_repeat.id)
            current_pos += len(current_unit)        

        if len(max_bed_tr.units) == 0:
            # No unit in the repeat matched the consensus unit, set end position to begin
            max_bed_tr.end = max_bed_tr.begin
        else:
            # calculate end position of the longest stretch, set purity metric
            max_bed_tr.end = max_bed_tr.begin + len(max_bed_tr.consensus_unit) * len(max_bed_tr.units)
            max_bed_tr.set_purity()

        return max_bed_tr


class BedTR(object):   
    supported_formats =  {"GangSTR"}

    def __init__(self, chromosome: str, consensus_unit: str, db_tr: int, out_format=None) -> None:
        self.chromosome = chromosome
        self.consensus_unit = consensus_unit
        self.db_tr = db_tr
        if not out_format:
            self.set_out_format("GangSTR")
        else:
            self.set_out_format(out_format)
        self.units = []
        self.begin = None
        self.end = None
        self.purity = None

    def set_out_format(self, new_format: str) -> None:
        if not new_format in self.supported_formats:
            raise ValueError(f"Specified output format '{new_format}' is not supported")
        self.out_format = new_format

    def get_bed_line(self) -> str:
        if self.out_format == "GangSTR":
            return "\t".join([
                        str(self.chromosome),
                        str(self.begin),
                        str(self.end),                    
                        str(len(self.consensus_unit)),
                        str(self.consensus_unit),
                        "\t", # GangSTR allows for optional 6th col with off-target regions. Just placeholder \t here for now
                        f"{self.db_tr}:{','.join(self.units)}:{self.purity}", # custom info field linking bed entry to DB
                        "\n"
                    ])

    def set_purity(self):
        tr_seq_len = 0
        n_mismatches = 0
        for unit in self.units:
            tr_seq_len += len(unit)
            for idx, nt in enumerate(unit):
                if not nt == self.consensus_unit[idx]:
                    n_mismatches += 1
        self.purity = round(1 - (n_mismatches / tr_seq_len), 2)
        

    def __str__(self):
        return f"BedTR from repeat '{self.db_tr}':\n{self.units}"

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-d", "--database", type=str, required=True, help="Path to database for which .bed file will be generated"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Output file will be generated here. Format will be 'bed-like' as described for the GangSTR\
            tandem repeat genotyping tool (https://github.com/gymreklab/gangstr)" 
    )

    return parser.parse_args()



def main():
    args = cla_parser()
    
    engine, session = connection_setup(args.database, echo=False)
    bedmaker = BedMaker(db_session=session, chromosomes=["auto", "sex", "mito"])
    
    thresholds = {i: 2 for i in range(1, 16)} # accept all repeats with at least 2 consensus-length motifs in tandem
    bedmaker.set_thresholds(thresholds)

    bedmaker.set_gene_selection()

    bedmaker.send_to_bed(args.output)

if __name__ == "__main__":
    main()
