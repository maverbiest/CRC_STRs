#!/usr/bin/env python3
"""
Classes for the creation, updating and querying of databases for the CRC TRs project

Contains classes: 
    DatabaseUtil, DatabaseConstructor, DatabaseUpdater, DatabaseReader

Author: Max Verbiest
Contact: max.verbiest@zhaw.ch
"""

import os
import sqlite3
import contextlib

import gtfparse


class DatabaseUtil(object):
    """ Abstract class that will be subclassed by more specific database utility classes.
    Implements execution of SQL commands, in a way that automatically commits and allows for fetching. 
    No checking for constraints on columns or checking of inputs is performed. This will be left to subclasses or user. 
    """

    def __init__(self, db_handle):
        """ All classes must be initialized using an existing SQLite database. The only exception is the DatabaseConstructor
        subclass, which will be used to create a database at the specified handle.
        Upon initialization, a connection to the database will be made (self.connection) that will remain open until the
        object is garbage collected.
        """
        self.db_handle = db_handle
        if not os.path.isfile(self.db_handle):
            if not isinstance(self, DatabaseConstructor):
                raise FileNotFoundError("Database '{}' was not found".format(self.db_handle))
        self.connection = sqlite3.connect(self.db_handle)    

    def execute_sql(self, sql, parameters=tuple(), fetch=False):
        """ Execute SQL command through self.connection
        Parameters
        sql (str):      SQL command to be executed
        parameters (tuple):     
                        parameter tuple that will be inserted into sql. Must match the number of '?' in sql
                        or will fail. Even single parameters must be specified as tuple e.g. (x, )
        fetch (bool)    If True, fetchall() of cursor will be returned after sql execution
        """
        with self.connection as conn: # auto-commits
            with contextlib.closing(conn.cursor()) as cursor: # auto-closes
                cursor.execute(sql, parameters)
                if fetch:
                    return cursor.fetchall()

    def execute_many_sql(self, sql, parameter_list=[], fetch=False):
        """ Execute many SQL command through self.connection
        Parameters
        sql (str):      SQL command to be executed more than once
        parameter_list (list(tuple)):     
                        List of parameter tuples that will be inserted into sql. sql will be executed once for each
                        tuple in the list. Each tuple must match the number of '?' in sql or will fail.
                        Even single parameters must be specified as tuple e.g. [(x, ), (y, ) ...]
        fetch (bool)    If True, fetchall() of cursor will be returned after sql execution
        """
        with self.connection as conn: # auto-commits
            with contextlib.closing(conn.cursor()) as cursor: # auto-closes
                cursor.executemany(sql, parameter_list)
                if fetch:
                    return cursor.fetchall()

    def __str__(self):
        return "{} instance connected to '{}'".format(type(self).__name__, self.db_handle)

    def __repr__(self):
        return "{}(db_handle='{}')".format(type(self).__name__, self.db_handle)


class DatabaseConstructor(DatabaseUtil):
    """ Class used to construct a database, but not populate it.
    Provides SQLite commands for table creation. Can create new database, or add missing tables
    to existing database. Methods ending in '_sql()' just return formatted SQL strings that
    can then be executed.
    """
    def __init__(self, db_handle):
        """Connection to new or existing database is made, check which tables need to be added"""
        super().__init__(db_handle=db_handle)
        self.to_create = self.get_missing_tables()

    def get_missing_tables(self):
        """Check which tables out of the required tables are missing in the connected database
        
        Returns
        required_tables.difference(detected_tables) (set):
                        The names of tables that are required, but are not currently in the connected database
        """
        required_tables = {"genes", "transcripts", "trs", "trs_transcripts", "exons", "exons_transcripts"}

        detected_tables = self.execute_sql("SELECT name FROM sqlite_master WHERE type='table';", fetch=True)
        detected_tables = {table[0] for table in detected_tables}

        return required_tables.difference(detected_tables)

    def construct_db(self):
        """Create all missing tables in the connected database"""
        if not self.to_create:
            print("DB at '{}' appears up-to-date, no further tables created".format(self.db_handle))
            return

        if "genes" in self.to_create:
            self.execute_sql(self.create_genes_table_sql())
            print("Created genes table")
        if "transcripts" in self.to_create:
            self.execute_sql(self.create_transcripts_table_sql())
            print("Created transcripts table")
        if "trs" in self.to_create:
            self.execute_sql(self.create_trs_table_sql())
            print("Created trs table")
        if "trs_transcripts" in self.to_create:
            self.execute_sql(self.create_trs_transcripts_table_sql())
            print("Created trs_transcripts table")
        if "exons" in self.to_create:
            self.execute_sql(self.create_exons_table_sql())
            print("Created exons table")
        if "exons_transcripts" in self.to_create:
            self.execute_sql(self.create_exons_transcripts_table_sql())
            print("Created exons_transcripts table")          
    
    def create_genes_table_sql(self):
        """ Create genes table
        ensembl_gene (str): Ensembl identifier of genomic region
        chromosome (str):   Chromosome identifier of type 'chr_[0-22]' or 'chr_{X,Y,M}'
        strand (str):       Which strand the gene is located on
        begin (int):        Position of the genes first nucleotide on the chromosome
        end (int):          Posistion of the genes last nucleotide on the chromosome
        upstream (int):     How many nucleotides upstream of the gene (i.e. promoter) were included (if any)        
        """
        return """
        CREATE TABLE genes
            (
                id          INTEGER PRIMARY KEY,
                ensembl_gene 
                            TEXT    NOT NULL    UNIQUE,
                chromosome  TEXT    NOT NULL,
                strand      TEXT    NOT NULL    
                                CHECK (strand IN ('fw', 'rv')),
                begin       INTEGER NOT NULL,
                end         INTEGER NOT NULL                    
            );
        """
        

    def create_transcripts_table_sql(self):
        return """
        CREATE TABLE transcripts
            (
                id          INTEGER PRIMARY KEY,
                gene_id     INTEGER NOT NULL,
                ensembl_transcript
                            TEXT    NOT NULL    UNIQUE,
                begin       INTEGER NOT NULL,
                end         INTEGER NOT NULL,
                FOREIGN KEY (gene_id)
                    REFERENCES genes(id)
            );
        """
    
    def create_trs_table_sql(self):
        """ Create TRs table
        The combination of gene_id and TR begin are enforced to be unique, i.e. no two TRs are allowed to
        begin at the same position in a gene. This means that TRs should be clustered (using TRAL) before 
        inserting them into the database
        """
        return """
        CREATE TABLE trs
            (
                id          INTEGER PRIMARY KEY,
                gene_id     INTEGER NOT NULL,
                begin       INTEGER NOT NULL,
                end         INTEGER NOT NULL,
                l_effective REAL    NOT NULL,
                n_effective REAL    NOT NULL,
                region_length
                            INTEGER NOT NULL,
                pvalue      REAL,
                divergence  REAL,
                FOREIGN KEY (gene_id)
                    REFERENCES genes(id),
                UNIQUE (gene_id, begin)
            );
        """

    def create_trs_transcripts_table_sql(self):
        return """
        CREATE TABLE trs_transcripts
            (
                tr_id       INTEGER NOT NULL,
                transcript_id     
                            INTEGER NOT NULL,
                FOREIGN KEY (tr_id)
                    REFERENCES trs(id),                                        
                FOREIGN KEY (transcript_id)
                    REFERENCES transcripts(id)
            );
        """

    def create_exons_table_sql(self):
        """ Create exons table
        start_cd (int):    first position of the start codon (if any)
        stop_cd (int):     first position of the stop codon (if any)
        """
        return """
        CREATE TABLE exons
            (
                id          INTEGER PRIMARY KEY,
                ensembl_exon
                            TEXT    NOT NULL    UNIQUE,
                begin       INTEGER NOT NULL,
                end         INTEGER NOT NULL,
                cds         BOOLEAN NOT NULL    
                                CHECK (cds IN (0, 1)),
                start_cd    INTEGER,
                stop_cd     INTEGER       
            );
        """

    def create_exons_transcripts_table_sql(self):
        return """
        CREATE TABLE exons_transcripts
            (
                exon_id     INTEGER NOT NULL,
                transcript_id     
                            INTEGER NOT NULL,
                FOREIGN KEY (exon_id)
                    REFERENCES exons(id),                                        
                FOREIGN KEY (transcript_id)
                    REFERENCES transcripts(id)
            );
        """

class DatabaseUpdater(DatabaseUtil):
    """ Class used to insert into (update) an existing database.
    Provides SQLite commands for updating tables. Methods ending in '*_sql()' just return formatted SQL 
    strings that can then be executed. Also contains 'first time setup' method (gtf_to_sqlite()) to convert 
    a gtf genome annotation file into a SQLite database.
    """
    def __init__(self, db_handle, gtf=None):
        """ Initialize DatabaseUpdater
        Parameters:
        db_handle (str):    Handle of existing DB with all required tables (see DatabaseConstructor)
        gtf (str):          Handle of gtf genome annotation file
        """
        super().__init__(db_handle=db_handle)
        self.gtf=gtf

    def set_gtf(self, gtf_handle):
        """Set gtf file to use"""
        self.gtf = gtf_handle

    def gtf_to_sqlite(self, gtf_handle=None, protein_coding=True):
        """ Parse gtf file and populate sqlite database
        This method is only meant as a 'first time setup' to populate an empty database, and therefore uses executemany without error handling
        The executemany opertations will fail if constraints (e.g UNIQUE) set on the tables are violated, 
        causing all insertions into that specific table to be rolled back. 

        Parameters
        gtf_handle (str):       Handle of gtf genome annotation file that will be used to construct DB
        protein_coding (bool):  If True, only protein coding elements from the gtf file will be used
        """
        if not gtf_handle:
            if self.gtf:
                gtf_handle = self.gtf
            else:
                raise AttributeError("No GTF file was specified or linked to instance of DatabaseUpdater")
        gtf_df = gtfparse.read_gtf(gtf_handle)

        if protein_coding:
            # Select only protein coding genes from the gtf
            gtf_df = gtf_df.loc[(gtf_df["gene_type"] == "protein_coding")]

        gene_list = self.get_gene_list(gtf_df)        
        self.execute_many_sql(self.insert_gene_sql(), gene_list)

        transcript_list = self.get_transcript_list(gtf_df)        
        self.execute_many_sql(self.insert_transcript_sql(), transcript_list)

        # Get the list of exon values, as well as a dictionary mapping ensembl gene names to a list
        ## of database transcript ids (one exon can belong to >1 transcript)
        exon_list, exon_transcript_id_dict = self.get_exon_information(gtf_df)
        self.execute_many_sql(self.insert_exon_sql(), exon_list)

        for ensembl_exon, db_transcript_ids in exon_transcript_id_dict.items():    
            # Now, for each Ensembl exon name, the database exon id needs to be retrieved  
            exon_table_id = self.execute_sql(
                    "SELECT id FROM exons WHERE ensembl_exon == ?;", 
                    (ensembl_exon,),
                    fetch = True
                    )[0][0]
            # Insert the exon and transcript ids into exon_transcripts join table
            ## There can be multiple ids in transcript_ids
            for db_transcript_id in db_transcript_ids:
                self.execute_sql(
                    self.insert_exons_transcripts_sql(),
                    (exon_table_id, db_transcript_id)
                )

    def get_gene_list(self, df):
        """ Get desired field values for all genes from a gtf data frame, return as list fo tuples to be 
        used for SQLite executemany

        Parameters
        df (pd.DataFrame):      Pandas data frame of gtf genome annotation file. Specifically one produced
                                using gtfparse.read_gtf()
        
        Returns
        gene_list (list(tuple)):
                                A list of tuples where each tuple contains the relevant value for a gene from
                                the input data frame
        """
        gene_list = []
        features = ["gene_id", "seqname", "strand", "start", "end"]
        for index, row in df.loc[(df["feature"] == "gene")].iterrows():
            values = [row[feature] for feature in features]
            # Strands are listed as '+' and '-' in gtf files
            if values[2] == "+":
                values[2] = "fw"
            elif values[2] == "-":
                values[2] = "rv"
            else:
                raise ValueError("Unknown strand identifier encountered, gtf indexes may be off")
            gene_list.append(tuple(values))
        return gene_list
    
    def get_transcript_list(self, df):
        """ Get desired field values for all transcripts from a gtf data frame, return as list fo tuples to be 
        used for SQLite executemany

        Parameters
        df (pd.DataFrame):      Pandas data frame of gtf genome annotation file. Specifically one produced
                                using gtfparse.read_gtf()
        
        Returns
        transcript_list (list(tuple)):
                                A list of tuples where each tuple contains the relevant value for a transcript from
                                the input data frame
        """
        transcript_list = []
        features = ["transcript_id", "start", "end"]
        for index, row in df.loc[(df["feature"] == "transcript")].iterrows():            
            values = [row[feature] for feature in features]
            gene_ensembl = (row["gene_id"], )
            result = self.execute_sql(                    
                    sql = "SELECT id FROM genes WHERE ensembl_gene == ?;",
                    fetch = True,
                    parameters = gene_ensembl                    
                )
            gene_id = result[0][0]

            values.insert(0, gene_id)
            transcript_list.append(tuple(values))
        return transcript_list

    def get_exon_information(self, df):
        """ Get desired field values for all exons from a gtf data frame, return as list fo tuples to be 
        used for SQLite executemany. Exons are trickier than genes or transcripts as their information is 
        spread across multiple rows. Thus, a dictionary is constructed containing relevant info for each exon,
        which is later converted to a list of tuples. Also return a dictionary mapping exon ensembl names (exon_id)
        to database transcript.id values, as this will be needed to populate the transcripts_exons join table later

        Parameters
        df (pd.DataFrame):      Pandas data frame of gtf genome annotation file. Specifically one produced
                                using gtfparse.read_gtf()
        
        Returns
        transcript_list (list(tuple)):
                                A list of tuples where each tuple contains the relevant value for a transcript from
                                the input data frame
        """        
        exon_dict = dict()
        exon_transcript_id_dict = dict()
        for index, row in df[df["feature"].isin({"exon", "CDS", "start_codon", "stop_codon"})].iterrows():            
            exon_id = row["exon_id"]
            if row["feature"] == "exon":                
                exon_dict[exon_id] = {
                    "begin": row["start"], 
                    "end": row["end"], 
                    "cds": False, 
                    "start_cd": None, 
                    "stop_cd": None
                }   

                transcript_id = self.execute_sql(
                    "SELECT id FROM transcripts WHERE ensembl_transcript == ?;", 
                    (row["transcript_id"],),
                    fetch = True
                    )[0][0]
                try:
                    exon_transcript_id_dict[exon_id].append(transcript_id)
                except KeyError:
                    exon_transcript_id_dict[exon_id] = [transcript_id]

            elif row["feature"] == "CDS":
                exon_dict[exon_id]["cds"] = True

            elif row["feature"] == "start_codon":
                exon_dict[exon_id]["start_cd"] = row["start"]

            elif row["feature"] == "stop_codon":
                exon_dict[exon_id]["stop_cd"] = row["start"]
        
        exon_list = []
        features = ["begin", "end", "cds", "start_cd", "stop_cd"]
        for exon_id in exon_dict.keys():
            values =  [exon_dict[exon_id][feature] for feature in features] 
            # values =  [v for k, v in exon_dict[exon_id].items()] # seems to preserve proper order but feels unsafe
            exon_list.append((exon_id, *values))
        return exon_list, exon_transcript_id_dict

    def insert_gene_sql(self):
        return """
        INSERT INTO genes (ensembl_gene, chromosome, strand, begin, end)  
        VALUES (?, ?, ?, ?, ?);
        """

    def insert_transcript_sql(self):
        return """
        INSERT INTO transcripts (gene_id, ensembl_transcript, begin, end)  
        VALUES (?, ?, ?, ?);
        """

    def insert_exon_sql(self, start_cd=None, stop_cd=None):
        """ An exon can belong to multiple transcripts, therefore 'INSERT OR IGNORE' because
        one exon can appear in a gtf file more than once
        """
        return """
        INSERT OR IGNORE INTO exons (ensembl_exon, begin, end, cds, start_cd, stop_cd)
        VALUES (?, ?, ?, ?, ?, ?);
        """

    def insert_exons_transcripts_sql(self):
        return """
        INSERT INTO exons_transcripts (exon_id, transcript_id) 
        VALUES (?, ?);
        """
        
    def insert_tr(self):
        raise NotImplementedError()     

    def __repr__(self):
        return "{}(db_handle='{}', gtf_handle={})".format(type(self).__name__, self.db_handle, self.gtf)


class DatabaseReader(DatabaseUtil):
    """ Class used to read information from an existing database.
    Mosly used to collect more complex methods and SQLite commands that are often needed when interacting
    with the databases in this project.
    """

    def get_gene_transcripts(self, gene):
        """ Get all transcripts belonging to a gene
        Paramters
        gene (str):     Ensemble ID of gene for which all transcripts will be returned         
        
        Returns
        transcripts (list(tuple)):
                        List of tuples where each tuple contains all information in the database
                        associated with a transcript associated with input gene
        """
        if isinstance(gene, str):
            gene = (gene, )
        transcripts = self.execute_sql(self.get_gene_transcripts_sql(), parameters=gene, fetch=True)
        return transcripts

    def get_gene_transcripts_sql(self):
        return """
        SELECT * FROM transcripts WHERE gene_id = (SELECT id FROM genes WHERE ensembl_gene = ?);
        """
    
    def get_transcript_exons(self, transcript, protein_coding=False):
        """ Get all exons belonging to a transcript, sorted on exon start position
        Parameters
        transcript (str):       Ensemble ID of transcript for which all exons will be returned
        protein_coding (bool):  If True, only exons that contain coding sequence will be returned: exons
                                that are purely 5' UTR are filtered out

        Returns
        exons (list(tuple)):    List of tuples where each tuple contains all information in the database
                                associated with an exon associated with input transcript (sorted on exon start site)          
        """
        if isinstance(transcript, str):
            transcript = (transcript, )   
        exons = self.execute_sql(self.get_transcript_exons_sql(), parameters=transcript, fetch=True)
        if protein_coding:
            # the [4] column for each exon will contain a 1 if the exon contains a CDS, otherwise 0
            exons = [exon for exon in exons if bool(exon[4])]
        strand = self.execute_sql(
            "SELECT strand FROM genes WHERE id = (SELECT gene_id FROM transcripts WHERE ensembl_transcript = ?);", 
            transcript, 
            fetch=True)[0][0]
        
        if strand == "fw":
            return list(sorted(exons, key=lambda x : x[2]))
        elif strand == "rv":
            return list(sorted(exons, key=lambda x : x[2], reverse=True))

    def get_transcript_exons_sql(self):
        return """
        SELECT exons.* FROM transcripts
        JOIN exons_transcripts ON transcripts.id = exons_transcripts.transcript_id
        JOIN exons ON exons_transcripts.exon_id = exons.id
        WHERE transcripts.ensembl_transcript = ?;
        """
        

def main():
    db_handle = "CRC_STRs/results/test/db/test.db"

    # constructor = DatabaseConstructor(db_handle=db_handle)
    # constructor.construct_db()

    # gtf = "/cfs/earth/scratch/verb/projects/CRC_STRs/data/test/genome_annot/gencode_small.gtf"
    # updater = DatabaseUpdater(db_handle=db_handle, gtf=gtf)
    # updater.gtf_to_sqlite()

    reader = DatabaseReader(db_handle=db_handle)
    transcript = reader.execute_sql("SELECT ensembl_transcript FROM transcripts;", fetch=True)[0]
    exons = reader.get_transcript_exons(transcript=transcript, protein_coding=True)
    for i in exons:
        print("{} -> {}".format(i[2], i[3]))


if __name__ == "__main__":
    main()
