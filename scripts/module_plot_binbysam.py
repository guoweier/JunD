# Weier Guo Python
# Create: 01/13/2025
# INTRODUCTION: plotting binbysam relative read coverage for all samples. 
#### UPDATE: to include all matplotlib parameters. 

## PACKAGES ##
import sys, math, os, time
import argparse 
from matplotlib import pyplot as plt 
import pandas as pd

def create_dir(outputdir):
    if os.path.isdir(outputdir) is False:
         os.system(f"mkdir {outputdir}")

def get_dataframe(input):
    df = pd.read_csv(input, delimiter="\t")
    return df 

def draw_figures(input, figwidth, figheight, figalpha, outputdir):
    df = get_dataframe(input)
    chroms = df["Chrom"].unique()
    colnames = list(df.columns)
    covcols = [value for value in colnames if "/" in value]
    for col in covcols:
        samplename = col.split("/")[0]
        for chrom in chroms:
            subset = df[df["Chrom"]==chrom]
            plt.figure(figsize=(figwidth, figheight))
            plt.scatter(subset["Strt"], subset[col], alpha=figalpha)
            plt.xlabel(f"{chrom}")
            plt.ylabel(f"{col}")
            plt.tight_layout()

            filename = f"{outputdir}/{samplename}_{chrom}.png"
            plt.savefig(filename)
            plt.close()


## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="plot binbysam")
    parser.add_argument("--input", type=str, default="binbysam.txt", help="Input binbysam file.")
    parser.add_argument("--outputdir", type=str, default="binbysam_plot", help="Output directory.")
    parser.add_argument("--figwidth", type=float, default=6, help="Each figure width (default=6).")
    parser.add_argument("--figheight", type=float, default=3, help="Each figure height (default=3).")
    parser.add_argument("--figalpha", type=float, default=0.8, help="Each dot alpha (default=0.7).")
    return parser.parse_args()

if __name__ in "__main__":
    args = parse_arguments()
    create_dir(args.outputdir)
    draw_figures(args.input, args.figwidth, args.figheight, args.figalpha, args.outputdir)