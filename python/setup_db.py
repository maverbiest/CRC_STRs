#!/usr/bin/env python3

import argparse

from database_utils import database_utils

def cla_parser():
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "--gtf", "-g", type=str, required=True, help="GTF file to convert into SQLite database"
    )
    parser.add_argument(
        "--output", "-o", type=str, required=True, help="Output handle where SQLite database will be generated"
    )

    return parser.parse_args()


def main():
    args = cla_parser()

    constructor = database_utils.DatabaseConstructor(db_handle=args.output)
    constructor.construct_db()

    updater = database_utils.DatabaseUpdater(db_handle=args.output, gtf=args.gtf)
    updater.gtf_to_sqlite()


if __name__ == "__main__":
    main()
