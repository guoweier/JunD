# Weier Guo Python
# Create: 01/02/2025
# INTRODUCTION: search bins with chimeric reads. 

## PACKAGES ##
import sys, math, os, time
import argparse 

## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="search chireads")
    parser.add_argument("--sampath", type=str, default="sam/", help="Path for sam files (default=sam/).")
    parser.add_argument("--controlfile", type=str, default="NA", help="Input sam file (default=NA).")
    parser.add_argument("--binbysamfile", type=str, default="binbysam.txt", help="Output bin file.")
    parser.add_argument("--binsize", type=int, default=100000, help="Bin size.")
    parser.add_argument("--breaks", type=str, default = False, help="Insert breaks.")
    parser.add_argument("--ploidy", type = int, default=2, help="Ploidy multiplier (default=2).")
    return parser.parse_args()