#!/usr/bin/env python3
import argparse
import os

import gtfparse
from sqlalchemy.engine import create_engine
from sqlalchemy.orm import declarative_base, relationship, sessionmaker
from sqlalchemy import Column, Integer, String, ForeignKey

Base = declarative_base()

class Gene(Base):
    __tablename__ = "genes"

    id = Column(Integer, primary_key=True)
    ensembl_gene = Column(String, nullable=False)
    chromosome = Column(String, nullable=False)
    strand = Column(String, nullable=False)
    begin = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)    

    def __repr__(self):
        return "Gene(ensembl_gene={}, chromosome={}, strand={}, begin={}, end={})".format(
            self.ensembl_gene,
            self.chromosome,
            self.strand,
            self.begin,
            self.end
        )


class Transcript(Base):
    __tablename__ = "transcripts"

    id = Column(Integer, primary_key=True)
    ensembl_transcript = Column(String, nullable=False)
    begin = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)

    gene_id = Column(Integer, ForeignKey("genes.id"))
    gene = relationship("Gene", back_populates="transcripts")

    def __repr__(self):
        return "Transcript(ensembl_transcript={}, begin={}, end={})".format(
            self.ensembl_transcript,
            self.begin,
            self.end
        )

# Add relationship directive to Gene class
Gene.transcripts = relationship("Transcript", order_by=Transcript.id, back_populates="gene")

def main():
    db_handle = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/db/test.db"

    if os.path.exists(db_handle):
        raise ValueError("A database already exists at specified handle, exiting!")

    engine = create_engine("sqlite:///{}".format(db_handle), echo = True)
    Base.metadata.create_all(engine)
    
    # Session.configure(bind=engine)
    # gtf_handle = "/cfs/earth/scratch/verb/projects/CRC_STRs/data/test/genome_annot/gencode_small.gtf"


if __name__ == "__main__":
    main()
