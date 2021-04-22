#!/usr/bin/env python3
import argparse
import os
import sys

from sqlalchemy.engine import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import NoResultFound

from setup_db import Gene, Transcript, Exon, Repeat
from tral_pvalues import load_repeatlists
from misc.constants import UPSTREAM

def get_gene_from_repeatlist(session, file_name):
    chromosome = file_name.split("_")[0]
    if not chromosome.startswith("chr"):
        raise Exception(f"Could not determine chromosome for repeatlist from {file_name}")

    if "_fw" in file_name:
        strand = "fw"
    elif "_rv" in file_name:
        strand = "rv"
    else:
        raise Exception(f"Could not determine strand for repeatlist from file {file_name}")

    if strand =="fw":
        gene_begin = int(file_name.split("_")[1].split("-")[0]) + UPSTREAM
        gene_end = int(file_name.split("_")[1].split("-")[1])
    else:
        ### UNTESTED CODE ###
        gene_begin = int(file_name.split("_")[1].split("-")[0])
        gene_end = int(file_name.split("_")[1].split("-")[1]) - UPSTREAM

    gene = session.query(Gene).filter_by(chromosome=chromosome, strand=strand, begin=gene_begin, end=gene_end).one()
    return gene

def add_repeat(gene, repeat):
    # repeat contains coordinates relative to the extracted genomic region that it was detected in
    ## thus, need to remap to chromosomal coordinates before adding to DB
    if gene.strand == "fw":
        chrom_begin = (gene.begin - UPSTREAM) + repeat.begin - 1 # only valid for fw strand genes
    else:
        ### UNTESTED CODE ###
        repeat_end = repeat.begin + repeat.repeat_region_length - 1
        chrom_begin = (gene.end + UPSTREAM) - repeat_end + 1 # only valid for rv strand genes
    
    # initialize instance of database Repeat
    db_repeat = Repeat(
        source = repeat.TRD,
        msa = ",".join(repeat.msa), # convert msa from list() to ',' separated str()
        begin = chrom_begin,
        end = chrom_begin + repeat.repeat_region_length - 1, # calculate end position
        l_effective = repeat.l_effective,
        n_effective = repeat.n_effective,
        region_length = repeat.repeat_region_length,
        p_value = repeat.d_pvalue['phylo'],
        divergence = repeat.d_divergence['phylo']
    )
    
    # relate Repeat to gene
    gene.repeats.append(db_repeat)

    # relate Repeat to all Transcript(s) of the Gene where (partial) overlap exists
    for transcript in gene.transcripts:
        if db_repeat.begin >= transcript.begin and db_repeat.begin <= transcript.end:
            transcript.repeats.append(db_repeat)
        elif db_repeat.end <= transcript.end and db_repeat.end >= transcript.begin:
            transcript.repeats.append(db_repeat)

def main():
    # args = cla_parser()
    # db_handle = args.database
    # gtf_handle = args.gtf

    db_path = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/db/test_brca2.db"
    input_path = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/repeats/selection"
    # check if database exists
    if not os.path.exists(db_path):
        raise FileNotFoundError("No DB was found at the specified handle")

    # connect to DB, initialize and configure session
    engine = create_engine("sqlite:///{}".format(db_path), echo=False)   
    Session = sessionmaker(bind=engine)
    session = Session()

    for file_name, repeat_list in load_repeatlists(input_path):
        if not file_name.endswith("_filt.pickle"):
            continue
        # retrieve the Gene that the repeat_list belongs to
        gene = get_gene_from_repeatlist(session, file_name)

        # make db entry for each repeat and add relationships to Gene and Transcripts
        for repeat in repeat_list.repeats:
            add_repeat(gene, repeat)

    # for t in gene.transcripts:
    #     print(t, t.end - t.begin, len(t.repeats))
    # print(len([r for r in gene.repeats if r.end < gene.begin]))
    
    # commit to DB
    session.commit()

if __name__ == "__main__":
    main()
