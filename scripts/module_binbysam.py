# Weier Guo Python
# Create: 12/26/2024
# INTRODUCTION: Create relative read coverage for all sam files. 
#### UPDATE: to adjust script to binbysamv8, for including different mapping options. ####

## PACKAGES ##
import sys, math, os, time
import argparse 

## GLOBAL PARAMETERS ## 
remlist = []
alls = []
sizes = []
data = {}
lookup = {}
globalcount = {}
liblist = []

## READ ALL SAM FILES ##
def read_sam(sampath):
    li = os.listdir(sampath)
    #todo = filter(lambda x: x.endswith('.sam'), li)
    return li

## GET CHROM INFO FROM SAM HEADER ##
def get_chrom(sampath, binsize):
    todo = read_sam(sampath)
    f = open(sampath+todo[0])
    for line in f:
        if line[0] != "@":
            break 
        if line[:3] == "@SQ":
            temp = (line[:-1].replace('SN:','').replace('LN:','').split('\t')[1:])
            key2 = temp[0]
            alls.append(key2)
        sizes.append(int(temp[1]))
        lookup[key2] = int(temp[1])
        key1 = temp[0][3:]
        if key1 not in data:
            data[key1] = {}
        if key2 not in data[key1]:
            data[key1][key2] = {}
        key3 = range(int(temp[1])//int(binsize)+1)
        for mod in key3:
            if mod not in data[key1][key2]:
                data[key1][key2][mod] = {}
    f.close()
    return todo 

## ADD BLANK LINES FOR SEPARATE CHROMOSOMES ##
def add_blanks(binsize):
    numblanks = max(sizes)//binsize//10
    return numblanks 

## COUNT READS ##
def count_reads(sampath, binsize):
    todo = get_chrom(sampath, binsize)
    for file in todo:
        f = open(sampath+file)
        print(file)
        libname = file.split('.')[0].replace('_aln','')
        liblist.append(libname)
        globalcount[libname] = 0
        for line in f:
            if line[0] == "@":
                continue
            ln = line[:-1].split('\t')
            if ln[2] == '*':
                continue
            if ln[1] == "16":
                pos = int(ln[3])+len(ln[9])-1
            elif ln[1] in ["0","99","163","83","147"]:
                pos = int(ln[3])
            else:
                continue 
            key1 = ln[2][3:]
            key2 = ln[2]
            key3 = int(pos)//binsize
            try:
                data[key1][key2][key3][libname] += 1
            except:
                if libname not in data[key1][key2][key3]:
                    data[key1][key2][key3][libname] = 1
            globalcount[libname] += 1
        f.close()
    

## CREATE OUTPUT file ##
def write_outfile_header(sampath, controlfile):
    controlpath = sampath+controlfile
    control = controlpath.split('.')[0].replace('_aln','')
    controlkey = control.split('/')[1]
    header = ['Chrom', 'Strt', 'End']
    header += liblist
    header += map(lambda x: x+"/"+controlkey, liblist)
    o.write('\t'.join(header)+'\n') 
    return controlkey

def write_outfile_data(sampath, controlfile, binsize, ploidy, breaks):
    control = write_outfile_header(sampath, controlfile)
    for chrom in alls:
        part = chrom[3:]
        bins = sorted(data[part][chrom].keys())
        for modbin in bins:
            libdata = data[part][chrom][modbin]
            if modbin == bins[-1]:
                line = [chrom, modbin*binsize+1, lookup[chrom]]
            else:
                line = [chrom, modbin*binsize+1, (modbin+1)*binsize]
            pers = {}
            for x in liblist:
                try:
                    line.append(libdata[x])
                    pers[x] = libdata[x]/float(globalcount[x])
                except:
                    line.append(0)
                    pers[x] = 0.0
            sums = sum(pers.values())/float(len(pers.values()))
            for x in liblist:
                if control == "NA":
                    try:
                        line.append(round(pers[x]/sums*ploidy, 3))
                    except ZeroDivisionError:
                        line.append(0.0)
                else:
                    try:
                        line.append(round(pers[x]/pers[control]*ploidy, 3))
                    except ZeroDivisionError:
                        line.append('.')
            fline = map(lambda x: str(x), line)
            o.write('\t'.join(fline)+'\n')
        if breaks == "True":
            numblanks = add_blanks(binsize)
            if chrom != alls[-1]:
                o.write(''.join(map(lambda x: '\n', range(numblanks))))
    o.close()

## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="binbysam")
    parser.add_argument("--sampath", type=str, default="sam/", help="Path for sam files (default=sam/).")
    parser.add_argument("--controlfile", type=str, default="NA", help="Input sam file (default=NA).")
    parser.add_argument("--binbysamfile", type=str, default="binbysam.txt", help="Output bin file.")
    parser.add_argument("--binsize", type=int, default=100000, help="Bin size.")
    parser.add_argument("--breaks", type=str, default = "False", help="Insert breaks.")
    parser.add_argument("--ploidy", type = int, default=2, help="Ploidy multiplier (default=2).")
    return parser.parse_args()

if __name__ in "__main__":
    args = parse_arguments()
    count_reads(args.sampath, args.binsize)
    o = open(args.binbysamfile, 'w')
    write_outfile_data(args.sampath, args.controlfile, args.binsize, args.ploidy, args.breaks)

