# Weier Guo Python
# Create: 01/03/2025
# INTRODUCTION: run pseudo-junction coverage for getting a threshold2 for chimeric reads bin selection. 

## PACKAGES ##
import sys, math, os, time
import argparse 
from collections import defaultdict

## GLOBAL ARGUMENTS ##
chromlist = defaultdict()
bins = defaultdict(lambda: defaultdict(int))

## READ ALL SAM FILES ##
def read_sam(sampath):
    li = os.listdir(sampath)
    return li

def find_pseudo(sampath, binsize, minqual):
    fileset = read_sam(sampath)
    for fname in fileset:
        print(fname)
        f = open(sampath+fname)
        while 1:
            l = f.readline()
            if l == '':
                break 
            if l[0] == '@' and l[1] != "P":
                h = l.split()
                name = h[1].replace("SN:","")
                nlen = int(h[2].replace("LN:", ""))
                chromlist[name] = nlen 
                continue 
            if l[0] == "@" and l[1] == "P":
                continue 
            x = l.split("\t")
            tup = []
            tup.append([x[2], int(x[3]), int(x[4]), x[9]])
            if x[2][0] != "C":
                continue 
            else:
                if int(x[4]) < minqual:
                    continue
                else:
                    if (int(x[3]) // binsize) != (int(x[3]) + len(x[9])) // binsize:
                        name = x[2] + "_" + str((int(x[3]) // binsize + 1) * binsize)
                        bins[name][fname] += 1
                    else:
                        continue
        f.close()
    return bins 

def write_output(sampath, binsize, minqual, mincov):
    fileset = read_sam(sampath)
    bins = find_pseudo(sampath, binsize, minqual)
    vbins = bins.keys()
    vbins = sorted(vbins, key = lambda x: int(x.split('_')[1]))
    vbins = sorted(vbins, key = lambda x: x.split('_')[0])
    # write header
    head = ['Ref', 'Bin_Position']
    for fn in fileset:
        head.append(fn.split('_aln')[0])
    o.write('\t'.join(head)+'\n')
    # write content
    for item in vbins:
        tup = item.split('_')
        line = [tup[0], int(tup[1])]
        calls = []
        for file in fileset:
            temp = bins[item][file]
            line.append(temp)
            calls.append(temp)
        if sum(calls) >= mincov:
            oline = map(lambda x: str(x), line)
            o.write('\t'.join(oline)+'\n')

## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="Threshold2: find pseudo junction coverage.")
    parser.add_argument("--sampath", type=str, default="sam/", help="Path for sam files (default=sam/).")
    parser.add_argument("--pseudofile", type=str, help="Output pseudo junction coverage file.")
    parser.add_argument("--binsize", type=int, default=10000,  help="Bin size (default=10kb).")
    parser.add_argument("--minqual", type=int, default=10, help="Minimum mapping quality (default = 10).")
    parser.add_argument("--mincov", type=int, default=0, help="Minimum coverage for a bin to be in the output, across all libraries (default=0).")
    return parser.parse_args()

if __name__ in "__main__":
    args = parse_arguments()
    o = open(args.pseudofile, 'w')
    write_output(args.sampath, args.binsize, args.minqual, args.mincov)
