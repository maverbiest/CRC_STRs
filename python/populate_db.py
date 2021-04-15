#!/usr/bin/env python3
import argparse
import os

import gtfparse
from sqlalchemy.engine import create_engine
from sqlalchemy.orm import sessionmaker

from setup_db import Gene, Transcript

# Session = sessionmaker()

def get_genome_annotations(gtf_handle, protein_coding=True):
    """ Parsing and optional filtering of gtf genome annotation file into pd.DataFrame
    """
    if not os.path.exists(gtf_handle):
        raise FileNotFoundError("No gtf genome annotation file was found at specified handle")
    gtf_df = gtfparse.read_gtf(gtf_handle)

    if protein_coding:
        # Select only protein coding genes from the gtf
        gtf_df = gtf_df.loc[(gtf_df["gene_type"] == "protein_coding")]

    return gtf_df


def add_genes(session, gtf_df):
    """ Get desired field values for all genes from a gtf data frame. For each gene in the gtf file,
    create a setup_db.Gene instance. Finally, add all encountered Genes to the session

    Parameters 
    session:        A SQLAlchemy session that is connected to a database where information from the genome
                    annotation file will be stored
    gtf_df (pd.DataFrame):  
                    Pandas data frame of gtf genome annotation file. Specifically one produced
                    using gtfparse.read_gtf()
    """
    gene_list = []
    features = ["gene_id", "seqname", "strand", "start", "end"]
    for index, row in gtf_df.loc[(gtf_df["feature"] == "gene")].iterrows():
        # Create dictionary mapping features to values, may be excessive however it will prevent 
        ## accessing the strand value by index later as it needs to be modified
        values = {feature: row[feature] for feature in features}
        # Strands are listed as '+' and '-' in gtf files
        if values["strand"] == "+":
            values["strand"] = "fw"
        elif values["strand"] == "-":
            values["strand"] = "rv"
        else:
            raise ValueError("Unknown strand identifier encountered, gtf indexes may be off")
        gene_list.append(
            Gene(
                ensembl_gene = values["gene_id"],
                chromosome = values["seqname"],
                strand = values["strand"],
                begin = values["start"],
                end = values["end"]
            )
        )
    session.add_all(gene_list)


def add_transcripts(session, gtf_df):
    """ Get desired field values for all transcripts from a gtf data frame, add them to their respective
    genes in the session

    Parameters
    session:        A SQLAlchemy session that is connected to a database where information from the genome
                    annotation file will be stored
    gtf_df (pd.DataFrame):  
                    Pandas data frame of gtf genome annotation file. Specifically one produced
                    using gtfparse.read_gtf()
    """
    for index, row in gtf_df.loc[(gtf_df["feature"] == "transcript")].iterrows():     
        ensembl_gene = row["gene_id"]
        gene = session.query(Gene).filter_by(ensembl_gene = ensembl_gene).one()
        gene.transcripts.append(
            Transcript(
                ensembl_transcript=row["transcript_id"],
                begin=row["start"],
                end=row["end"]
            )
        )


def main():
    db_handle = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/db/test.db"

    # check if database exists
    if not os.path.exists(db_handle):
        raise FileNotFoundError("No DB was found at the specified handle")

    # connect to DB, configure and initialize session
    engine = create_engine("sqlite:///{}".format(db_handle), echo=False)
    
    Session = sessionmaker(bind=engine)
    session = Session()

    gtf_handle = "/cfs/earth/scratch/verb/projects/CRC_STRs/data/test/genome_annot/gencode_small.gtf"
    gtf_df = get_genome_annotations(gtf_handle, protein_coding=True)

    add_genes(session, gtf_df)
    add_transcripts(session, gtf_df)
    
    for i in session.query(Gene).all():
        print(i)
        print(i.transcripts)
        print()

    session.rollback()


if __name__ == "__main__":
    main()
