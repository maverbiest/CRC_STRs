#!/usr/bin/env python3
import argparse
import os

import gtfparse
from sqlalchemy.engine import create_engine
from sqlalchemy.orm import sessionmaker

from setup_db import Gene, Transcript

Session = sessionmaker()

def main():
    db_handle = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/db/test.db"

    if not os.path.exists(db_handle):
        raise FileNotFoundError("No DB was found at the specified handle")

    # connect to DB, configure Session and initialize
    engine = create_engine("sqlite:///{}".format(db_handle), echo=False)
    Session.configure(bind=engine)
    session = Session()

    test_gene = Gene(ensembl_gene="ENSG00000279928.1", chromosome="chr1", strand="fw", begin=182393, end=184158)
    test_transcript = Transcript(ensembl_transcript="ENST00000624431.1", begin=182393, end=184158)
    test_gene.transcripts = [test_transcript]
    session.add(test_gene)

    res = session.query(Gene).filter_by(ensembl_gene = "ENSG00000279928.1").one()

    print(res)
    print(res.transcripts)

    session.rollback()


if __name__ == "__main__":
    main()
