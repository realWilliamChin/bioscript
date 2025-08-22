#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import os, sys
import argparse
import pysam
from loguru import logger


def parse_input():
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', dest="infile", help="input fasta file")
    parser.add_argument(
        "--outputs", nargs="+", help="list of output files", required=True
    )
    args = parser.parse_args()
    return args


def split_fasta(infile, outputs):
    NIDS = len(outputs)
    fasta = pysam.FastaFile(infile)
    outs = [open(f, "w") for f in outputs]
    outidx = 0
    for name in fasta.references:
        seq = fasta.fetch(name)
        outs[outidx].write(">{}\n{}\n".format(name, seq))
        outidx += 1
        if outidx == NIDS:
            outidx = 0
    for out in outs:
        out.close()


def main():
    args = parse_input()
    split_fasta(args.infile, args.outputs)


if __name__ == "__main__":
    main()