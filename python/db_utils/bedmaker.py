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
    # Default threshold values for STRs, taken from Lai & Sun, 2003
    # format: {unit length: minimal number of units}
    default_thresholds = {
        1: 9, 
        2: 4,
        3: 4,
        4: 4,
        5: 4,
        6: 4
    }
    # All allowed values for chromosomes
    allowed_chromosomes = {f"chr{i}" for i in range(1, 23)}.union({"chrX", "chrY", "chrM"})
   
    def __init__(self, 
                db_session, 
                # chromosomes: Iterable[str]={"auto"}, 
                consensus_only: bool=False, 
                thresholds: dict=default_thresholds) -> None:
        self.db_session = db_session
        # self.set_chroms(chromosomes)    
        self.consensus_only = consensus_only
        self.set_thresholds(thresholds)
        # if not thresholds:            
        #     self.thresholds = self.default_thresholds
        # else:
        #     self.set_thresholds(thresholds)        
        self.gene_selection: Optional[dict] = None  # {gene.id: gene.strand}

    def set_thresholds(self, thresh_dict: dict) -> None:
        # Check if all keys and values in supplied dictionary are integers
        if not all(isinstance(i, int) for i in [*thresh_dict.keys(), *thresh_dict.values()]):
            raise ValueError("Thresholds must be specified as dict(int(unit_len): int(min_units))")
        self.thresholds = thresh_dict

    def set_gene_selection(self, chromosomes=None, genes=None) -> None:
        """ For each gene in the DB originating from a chromosome in 
        self.chromosomes, extract the gene.id and gene.strand to save in 
        self.gene_selection as key-value pairs {gene.id: gene, ... } 
        """
        if chromosomes and genes:
            raise ValueError("Specify either set of chromosomes or set of genes to set gene selection, not both.")
        self.gene_selection = dict()

        if chromosomes:
            if not len(chromosomes.difference(self.allowed_chromosomes)) == 0:
                raise ValueError(f"Unknown chromosome identifier(s) specified. Allowed chromosomes: {self.allowed_chromosomes}")
            for chromosome in chromosomes:
                # Collect all genes originating from the current chromosome, add them to gene_selection
                db_genes = self.db_session.query(Gene).filter_by(chromosome=chromosome).all()
                for db_gene in db_genes:
                    self.gene_selection[db_gene.id] = db_gene
        elif genes:
            for gene in genes:
                db_gene = self.db_session.query(Gene).filter_by(ensembl_gene=gene).one()
                self.gene_selection[db_gene.id] = db_gene
        else:
            raise ValueError("Specify either a set of target chromosomes or a set of target genes")


    def threshold_filter(self, consensus_thresh=False):
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
            for repeat in self.db_session.query(Repeat).filter_by(gene_id=gene.id):                
                bed_tr_list = self.get_bed_trs(repeat)
                for bed_tr in bed_tr_list:
                    if consensus_thresh:
                        bed_tr_len = bed_tr.longest_cs
                    else:
                        bed_tr_len = len(bed_tr.units)
                    try:                                                         
                        if bed_tr_len >= self.thresholds[len(bed_tr.consensus_unit)]:
                            yield bed_tr
                    except KeyError:
                        # repeat.l_effective is not part of the self.threshold range currently active
                        continue

    def send_to_bed(self, output_file: str, consensus_thresh=False) -> None:
        """ Generate an output file in bed format that can be used by
        GangSTR to genotype STRs in alignment data
        """

        with open(output_file, "w") as o:
            for bed_tr in self.threshold_filter(consensus_thresh=consensus_thresh):
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
        current_bed_tr = BedTR(self.gene_selection[db_repeat.gene_id].chromosome, consensus_unit, db_repeat.id, out_format="GangSTR")
        max_bed_tr = BedTR(self.gene_selection[db_repeat.gene_id].chromosome, consensus_unit, db_repeat.id, out_format="GangSTR")

        current_pos = db_repeat.begin
        
        for unit in msa.split(","):
            current_unit = unit.replace("-", "")
            if len(current_unit) == len(consensus_unit):    # current unit must be same length as consensus to be considered
                # if we only allow perfect units, check if current unit is a match. If not: reset
                if self.consensus_only and not current_unit == consensus_unit:
                    current_bed_tr = BedTR(self.gene_selection[db_repeat.gene_id].chromosome, consensus_unit, db_repeat.id, out_format="GangSTR")
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
                current_bed_tr = BedTR(self.gene_selection[db_repeat.gene_id].chromosome, consensus_unit, db_repeat.id, out_format="GangSTR")
            current_pos += len(current_unit)        

        if len(max_bed_tr.units) == 0:
            # No unit in the repeat matched the consensus unit, set end position to begin
            max_bed_tr.end = max_bed_tr.begin
        else:
            # calculate end position of the longest stretch, set purity metric
            max_bed_tr.end = max_bed_tr.begin + len(max_bed_tr.consensus_unit) * len(max_bed_tr.units) - 1
            max_bed_tr.set_purity()

        return max_bed_tr

    def get_bed_trs(self, db_repeat):
        """ Get all stretches of units with len == len(consensus_unit) from the input db_repeat. Will even
        return BedTRs of single units and leave all filtering up to the calling function.

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
        bed_tr_list = []
        current_bed_tr = BedTR(self.gene_selection[db_repeat.gene_id].chromosome, consensus_unit, db_repeat.id, out_format="GangSTR")
        current_pos = db_repeat.begin

        for unit in msa.split(","):            
            current_unit = unit.replace("-", "")
            if len(current_unit) == len(consensus_unit):    # current unit must be same length as consensus to be considered
                # if we only allow perfect units, check if current unit is a match. If not: reset
                if self.consensus_only and not current_unit == consensus_unit:
                    if len(current_bed_tr.units) > 0:
                        bed_tr_list.append(current_bed_tr)
                    current_bed_tr = BedTR(self.gene_selection[db_repeat.gene_id].chromosome, consensus_unit, db_repeat.id, out_format="GangSTR")
                    current_pos += len(current_unit)
                    continue
                if len(current_bed_tr.units) == 0:
                    current_bed_tr.begin = current_pos
                current_bed_tr.units.append(current_unit)
            else:
                if len(current_bed_tr.units) > 0:
                    bed_tr_list.append(current_bed_tr)
                current_bed_tr = BedTR(self.gene_selection[db_repeat.gene_id].chromosome, consensus_unit, db_repeat.id, out_format="GangSTR")
            current_pos += len(current_unit)

        # Check if last entry in bed_tr_list equals the current_bed_tr
        if len(current_bed_tr.units) > 0:
            try:
                if current_bed_tr != bed_tr_list[-1]:
                    bed_tr_list.append(current_bed_tr)
            except IndexError:
                bed_tr_list.append(current_bed_tr)


        for bed_tr in bed_tr_list:
            if len(bed_tr.units) == 0:
                # No unit in the repeat matched the consensus unit, set end position to begin
                bed_tr.end = bed_tr.begin
                continue
            # calculate end position of the longest stretch, set purity metric
            bed_tr.end = bed_tr.begin + len(consensus_unit) * len(bed_tr.units) - 1
            bed_tr.set_purity()
            bed_tr.set_longest_cs_stretch()

        return bed_tr_list


class BedTR(object):   
    supported_formats =  {"GangSTR"}

    def __init__(self, chromosome: str, consensus_unit: str, db_tr: int, out_format: str) -> None:
        self.chromosome = chromosome
        self.consensus_unit = consensus_unit
        self.db_tr = db_tr
        self.set_out_format(out_format)
        self.units = []
        self.begin = None
        self.end = None
        self.purity = None
        self.longest_cs = None

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
                        f"{self.db_tr}:{','.join(self.units)}:{self.purity}:{self.longest_cs}", # custom info field linking bed entry to DB
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

    def set_longest_cs_stretch(self):
        if self.purity == 1.0:
            self.longest_cs = len(self.units)
            return
        cur_stretch = 0
        longest_stretch = 0
        for unit in self.units:            
            if unit == self.consensus_unit:                
                cur_stretch += 1
                if cur_stretch > longest_stretch:
                    longest_stretch = cur_stretch
                continue
            cur_stretch = 0
        self.longest_cs = longest_stretch
        
    def __str__(self):
        return f"BedTR from repeat '{self.db_tr}': {self.units}"

    def __eq__(self, other):
        if not isinstance(other, BedTR):
            return False
        if not self.db_tr == other.db_tr:
            return False
        if not self.begin == other.begin:
            return False
        if not self.end == other.end:
            return False
        if not self.consensus_unit == other.consensus_unit:
            return False
        return True
    
    def __ne__(self, other):
        return not self.__eq__(other)

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-d", "--database", type=str, required=True, help="Path to database for which .bed file will be generated"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Output file will be generated here. Format will be 'bed-like' as described for the GangSTR\
            tandem repeat genotyping tool (https://github.com/gymreklab/gangstr)" 
    )
    parser.add_argument(
        "-g", "--genes", type=str, required=False, help="File with ensembl gene names for which STR loci should be retrieved."
    )
    parser.add_argument(
        "-p", "--perfect", action='store_true', help="True/False flag to control whether only perfect STRs will be considered (default: False)"
    )

    return parser.parse_args()



def main():
    args = cla_parser()
    
    engine, session = connection_setup(args.database, echo=False)
    bedmaker = BedMaker(db_session=session, consensus_only=args.perfect)
    
    # thresholds = {i: 2 for i in range(1, 16)} # accept all repeats with at least 2 consensus-length motifs in tandem
    thresholds = {
        1: 9,
        2: 4,
        3: 4,
        4: 3,
        5: 3,
        6: 3
    }
    bedmaker.set_thresholds(thresholds)

    if args.genes:
        if not os.path.isfile(args.genes):
            raise ValueError(f"No file of gene identifiers found at '{args.genes}'")
        with open(args.genes, 'r') as f:
            target_genes = {line.strip() for line in f}
            bedmaker.set_gene_selection(genes=target_genes)
    else:
        bedmaker.set_gene_selection(chromosomes=bedmaker.allowed_chromosomes)


    bedmaker.send_to_bed(args.output, consensus_thresh=True)

if __name__ == "__main__":
    main()
