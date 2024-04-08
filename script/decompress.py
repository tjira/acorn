#!/usr/bin/env python

import argparse as ap, gzip as gz, shutil as sh

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn File Decompressor", description="File decompresor script for the Quantum Acorn package.", add_help=False)

    # add optional arguments
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")

    # add file arguments
    parser.add_argument("file", help="File to decompress.")

    # parse arguments
    args = parser.parse_args()

    # decompress the file
    sh.copyfileobj(gz.open(args.file, "rb"), open(args.file + ".gz", "wb")) # type: ignore

    # overwrite the original file
    sh.move(args.file + ".gz", args.file)
