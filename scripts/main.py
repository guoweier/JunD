# Weier Guo Python 
# Create: 12/2/2024
## usage: currently this scripts needs to be run in the folder with all fq files. 
#### UPDATE: to specify the location of fq files. #####

## PACKAGES ##
import sys, math, os, time
import argparse 
import subprocess 

def run_modules(scripts_and_args):
    for script, args in scripts_and_args:
        try:
            command = ["python3.12", script]+args
            print(f"Running: {' '.join(command)}")
            subprocess.run(command)
        except Exception as e:
            print(f"Failed to run {script}: {e}")

## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="Main")
    # calibrate reference
    #parser.add_argument("--reference", type=str, help="Input reference fasta file.")
    #parser.add_argument("--output", type=str, help="Output calibrated reference file.")
    # mapping 
    parser.add_argument("--database", type=str, help="Input database file for mapping (default=output from previous step).")
    ##### To Update: remove this parameter for automatic ref indexing #######
    #parser.add_argument("--index", type=str, default ='n', help="whether do reference indexing step. n=NO (default), y=YES.")
    parser.add_argument("--mode", type = str, default='u', help="When mapping as uninterleaved, Read type, u = UNINTERLEAVED (default), i = INTERLEAVED")
    parser.add_argument("--thread", type=int, default=8, help="How many threads to use during alignment.")
    # binbysam
    #parser.add_argument("--binbysam", type=str, default="y", help="Whether to perform binbysam. y=YES (default), n=NO.")
    parser.add_argument("--sampath", type=str, default="sam/", help="Path for sam files (default=sam/).")
    parser.add_argument("--controlfile", type=str, default="NA", help="Input sam file (default=NA).")
    parser.add_argument("--binbysamfile", type=str, required=True, help="Output bin file.")
    parser.add_argument("--binsize", type=int, default=100000, help="Bin size.")
    parser.add_argument("--breaks", type=str, default = False, help="Insert breaks.")
    parser.add_argument("--ploidy", type = int, default=2, help="Ploidy multiplier (default=2).")
    return parser.parse_args()

if __name__ in "__main__":
    args = parse_arguments()
    scripts_and_args = [
        #("../scripts/module_ref_calibrate.py", ["--reference", args.reference, "--output", args.output]),
        ("../scripts/module_fqtobam.py", ["--database", args.database, "--mode", args.mode, "--thread", str(args.thread)]),
        ("../scripts/module_binbysam.py", ["--sampath", args.sampath, "--controlfile", args.controlfile, "--binbysamfile", args.binbysamfile, "--binsize", str(args.binsize), "--breaks", args.breaks, "--ploidy", str(args.ploidy)])
    ]
    run_modules(scripts_and_args)