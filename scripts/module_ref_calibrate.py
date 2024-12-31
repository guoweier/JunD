# Weier Guo Python 
# Create: 12/2/2024
# Introduction: check reference fasta has the appropriate chrom names (Chr01, Chr02, etc)

## PACKAGES ##
import sys, math, os, time
import argparse 
import subprocess 

## CHECK REF NAME ##
def calibrate_chrname(reference, output):
    with open(reference, "r") as ref, open(output, "w") as outfile:
        seq_id = 1
        for line in ref:
            if line.startswith(">"):
                if "chromosome" in line:
                    outfile.write(f">Chr{seq_id}\n")
                    seq_id += 1
                elif line.split("\n")[0] == ">Chr"+str(seq_id):
                    outfile.write(line)
            else:
                outfile.write(line)


## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="Calibrate reference chromosomes.")
    parser.add_argument("--reference", type=str, help="Input reference fasta file.")
    parser.add_argument("--output", type=str, help="Output file.")
    return parser.parse_args()

if __name__ in "__main__":
    args = parse_arguments()
    calibrate_chrname(args.reference, args.output)
