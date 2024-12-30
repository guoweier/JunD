# Weier Guo Python 
# Create: 12/2/2024
# Introduction: check reference fasta has the appropriate chrom names (Chr01, Chr02, etc)

## PACKAGES ##
import sys, math, os, time
import argparse 





## CALL ARGUMENTS ##
def parse_arguments():
    parser.add_argument("--reference", type=str, help="Input reference fasta file.")
    return parser.parse_args()

