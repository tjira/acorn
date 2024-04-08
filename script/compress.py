#!/usr/bin/env python

import argparse as ap, gzip as gz, shutil as sh

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn File Compressor", description="File compresor script for the Quantum Acorn package.", add_help=False)

    # add optional arguments
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-l", "--level", type=int, default=9, help="Level of compression.")

    # add file arguments
    parser.add_argument("file", help="File to compress.")

    # parse arguments
    args = parser.parse_args()

    # compress the file
    sh.copyfileobj(open(args.file, "rb"), gz.open(args.file + ".gz", "wb", compresslevel=args.level)) # type: ignore

    # overwrite the original file
    sh.move(args.file + ".gz", args.file)
