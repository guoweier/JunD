# Weier Guo Python
# Create: 01/03/2025
# INTRODUCTION: set all control files chireads to a chosen threshold.

## PACKAGES ##
import sys, math, os, time
import argparse 

def get_samples(input):
    with open(input, "r") as chiread:
        head = chiread.readline()
        head = head.split("\n")[0]
        samples = head.split("\t")
    return samples 

def filter_control(input, control, number, t1file):
    samples = get_samples(input)
    for ct in control:
        if ct not in control:
            raise ValueError(f"Unsupported control: {ct}")
    ctindex = [i for i, value in enumerate(samples) if value in control]
    with open(input, "r") as chiread, open(t1file, "w") as o:
        header = chiread.readline()
        o.write(header)
        for line in chiread:
            line = line.split("\n")[0]
            ln = line.split("\t")
            values = [float(value) for i, value in enumerate(ln) if i in ctindex]
            if all(value <= float(number) for value in values):
                o.write(line+"\n")

## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="Threshold1: set controls threshold.")
    parser.add_argument("--input", type=str, help="Input chiread file.")
    parser.add_argument("--control", type=str, nargs="+", help="Input control names.")
    parser.add_argument("--number", type=int, default=0, help="Choose control threshold (default=0).")
    parser.add_argument("--t1file", type=str, help="Output file with controls=chosen threshold.")
    return parser.parse_args()

if __name__ in "__main__":
    args = args = parse_arguments()
    filter_control(args.input, args.control, args.number, args.t1file)