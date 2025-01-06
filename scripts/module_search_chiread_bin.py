# Weier Guo Python
# Create: 01/02/2025
# INTRODUCTION: search bins with chimeric reads. 

## PACKAGES ##
import sys, math, os, time
import argparse 
from collections import defaultdict
import itertools

## GLOBAL ARGUMENTS ##
chromlist = defaultdict()
bins = defaultdict(lambda: defaultdict(int))

## READ ALL SAM FILES ##
def read_sam(sampath):
    li = os.listdir(sampath)
    return li

def search_chiread(sampath, binsize, insert, minqual):
    fileset = read_sam(sampath)
    for fname in fileset:
        print(fname)
        matches = []
        f = open(sampath+fname)
        while 1:
            line = f.readline()
            if line == "":
                break
            if line[0] == '@' and line[1] != "P":
                h = line.split()
                name = h[1].replace("SN:","")
                nlen = int(h[2].replace("LN:",""))
                chromlist[name] = nlen
                continue
            if line[0] == "@" and line[1] == "P":
                continue
            x = line.split("\t")
            tup = []
            tup.append([x[2], int(x[3]), int(x[4])])
            templines = []
            while 1:
                temp = f.tell()
                nxt = f.readline()
                xn = nxt.split('\t')
                if xn[0] == x[0]:
                    tup.append([xn[2], int(xn[3]), int(xn[4])])
                    templines.append(xn)
                else:
                    f.seek(temp)
                    break
            combos = []
            for each in itertools.combinations(tup,2):
                if "Chr" not in each[0][0] or "Chr" not in each[1][0]:
                    continue
                if each[0][0] == each[1][0] and abs(each[0][1]-each[1][1]) < insert:
                    continue
                if each[0][2] < minqual or each[1][2] < minqual:
                    continue
                combos.append(each)
            for item in combos:
                matches.append(list(item))
        f.close()

        for i in matches:
            item = sorted(i, key = lambda x: x[0])
            chr1 = item[0][0]
            pos1 = item[0][1]
            key1 = str(pos1 // binsize)
            chr2 = item[1][0]
            pos2 = item[1][1]
            key2 = str(pos2 // binsize)
            bigkey = '-'.join([chr1,key1,chr2,key2])
            bigkey2 = '-'.join([chr2,key2,chr1,key1])
            bins[bigkey][fname] +=1
            bins[bigkey2][fname] +=1
    return bins 

def write_header(sampath):
    head = ['Ref1', 'Ref1-BinStart', 'Ref2', 'Ref2-BinStart']
    fileset = read_sam(sampath)
    for fn in fileset:
        head.append(fn.split('_aln')[0])
    o.write('\t'.join(head)+'\n')

def write_content(sampath, binsize, insert, minqual, mincov):
    fileset = read_sam(sampath)
    bins = search_chiread(sampath, binsize, insert, minqual)
    vbins = bins.keys()
    vbins = sorted(vbins, key = lambda x: int(x.split('-')[3]))
    vbins = sorted(vbins, key = lambda x: x.split('-')[2])
    vbins = sorted(vbins, key = lambda x: int(x.split('-')[1]))
    vbins = sorted(vbins, key = lambda x: x.split('-')[0])
    same = 0
    for item in vbins:
        tup = item.split('-')
        if int(tup[0][3]) > int(tup[2][3]):
            continue
        elif int(tup[0][3]) == int(tup[2][3]) and int(tup[1]) > int(tup[3]):
            continue 
        elif int(tup[0][3]) == int(tup[2][3]) and int(tup[1]) == int(tup[3]) and same != 0:
            same = 0
            continue
        else:
            if int(tup[0][3]) == int(tup[2][3]) and int(tup[1]) == int(tup[3]) and same == 0:
                same += 1
            line = [tup[0], int(tup[1])*binsize, tup[2], int(tup[3])*binsize]
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
    parser = argparse.ArgumentParser(description="search chireads")
    parser.add_argument("--sampath", type=str, default="sam/", help="Path for sam files (default=sam/).")
    parser.add_argument("--chireadfile", type=str, help="Output chimeric reads bin file.")
    parser.add_argument("--binsize", type=int, default=10000,  help="Bin size (default=10kb).")
    parser.add_argument("--minqual", type=int, default=10, help="Minimum mapping quality (default = 10).")
    parser.add_argument("--insert", type=int, default=2000, help="Maximum estimated insert between forward and reverse (default=2000).")
    parser.add_argument("--mincov", type=int, default=0, help="Minimum coverage for a bin to be in the output, across all libraries (default=0).")
    return parser.parse_args()

if __name__ in "__main__":
    args = parse_arguments()
    o = open(args.chireadfile, 'w')
    write_header(args.sampath)
    write_content(args.sampath, args.binsize, args.insert, args.minqual, args.mincov)
    o.close()