# Weier Guo Python
# Create: 01/08/2025
# INTRODUCTION: According to two rounds of threshold selection, get potential chiread bin file, select for seeds for PRICE assembly. 
#### UPDATE: currently only worked for single mapped reads. Should update it to paired-map reads, which can also look for cross-junction mated reads. 

## PACKAGES ##
import sys, math, os, time
import argparse 
import itertools

## GLOBAL PARAMETERS ##
list_name = {}
list_seq = {}
os.system("mkdir seeds")

def get_samplename(input):
    if '/' in input:
        filename = input[input.rindex('/')+1:]
    else:
        filename = input 
    samplename = filename.split("-")[0].split("T2_")[1]
    return samplename 

def extract_chiread(sampath, input, chr1, chr2, pos1, pos2, binsize, insert, minqual):
    matches = []
    samplename = get_samplename(input)
    li = os.listdir(sampath)
    fileset = list(filter(lambda x: x.endswith('.sam'), li))
    for sam in fileset:
        if samplename in sam:
            samfile = sam
            break 
    with open(sampath+samfile, "r") as f:
        while 1:
            l = f.readline()
            if l == '':
                break
            if l[0] == "@":
                continue
            x = l[:-1].split("\t")
            tup = []
            tup.append([x[0], x[2], int(x[3]), int(x[4]), x[9]])
            templines = []
            while 1:
                temp = f.tell()
                nxt = f.readline()
                xn = nxt[:-1].split('\t')
                if xn[0] == x[0]:
                    tup.append([xn[0], xn[2], int(xn[3]), int(xn[4]), xn[9]])
                    templines.append(xn)
                else:
                    f.seek(temp)
                    break
            combos = []
            for each in itertools.combinations(tup,2):
                if "Chr" not in each[0][1] or "Chr" not in each[1][1]:
                    continue
                if each[0][1] == each[1][1] and abs(each[0][2]-each[1][2]) < insert:
                    continue
                if each[0][3] < minqual or each[1][3] < minqual:
                    continue
                if chr1 == chr2:
                    if each[0][1] == each[1][1] == chr1:
                        if (each[0][2] > pos1 and each[0][2] < (pos1 + binsize)) and ((each[1][2] > pos2) and (each[1][2]) < (pos2 + binsize)):
                            combos.append(each)
                        elif (each[1][2] > pos1 and each[1][2] < (pos1 + binsize)) and ((each[0][2] > pos2) and (each[0][2]) < (pos2 + binsize)):
                            combos.append(each)
                elif chr1 != chr2:
                    if each[0][1] == chr1 and each[1][1] == chr2:
                        if (each[0][2] > pos1 and each[0][2] < (pos1 + binsize)) and ((each[1][2] > pos2) and (each[1][2]) < (pos2 + binsize)):
                            combos.append(each)
                    elif each[1][1] == chr1 and each[0][1] == chr2:
                        if (each[1][2] > pos1 and each[1][2] < (pos1 + binsize)) and ((each[0][2] > pos2) and (each[0][2]) < (pos2 + binsize)):
                            combos.append(each)
                else: 
                    continue
            for item in combos:
                matches.append(list(item))
    return matches 

def collect_chiread(sampath, input, binsize, insert, minqual):
    samplename = get_samplename(input)
    with open(input, "r") as chiread:
        head = chiread.readline()
        hd = head[:-1].split('\t')
        spindex = [i for i, value in enumerate(hd) if value == samplename]
        for line in chiread:
            print(line[:-1])
            ln = line[:-1].split('\t')
            chr1 = ln[0]
            pos1 = int(ln[1])
            chr2 = ln[2]
            pos2 = int(ln[3])
            seed = f"{chr1}_{pos1}_{chr2}_{pos2}"
            matches = extract_chiread(sampath, input, chr1, chr2, pos1, pos2, binsize, insert, minqual)
            with open(f"seeds/seeds_{seed}.fasta", "w") as o:
                for item in matches:
                    o.write(">"+item[0][0]+"\n"+item[0][4]+"\n")
                    o.write(">"+item[1][0]+"\n"+item[1][4]+"\n")
                #if seed not in list_seq:
                    #list_name[seed] = [item[0][0]]
                    #list_seq[seed] = [item[0][0],item[0][4],item[1][0],item[1][4]]
                #elif seed in list_seq:
                    #if item[0][0] not in list_name[seed]:
                        #list_name[seed] += [item[0][0]]
                        #list_seq[seed] += [item[0][0],item[0][4],item[1][0],item[1][4]]
    #return list_name, list_seq 

def output_seed(sampath, input, binsize, insert, minqual):
    list_name, list_seq = collect_chiread(sampath, input, binsize, insert, minqual)
    os.system("mkdir seeds")
    for i in list_name:
        with open(f"seeds/seeds_{i}.fasta", "w") as o:
            for a in range(len(list_seq[i])):
                o.write(">"+list_seq[i][a][0]+"\n"+list_seq[i][a][1]+"\n")
                o.write(">"+list_seq[i][a][2]+"\n"+list_seq[i][a][3]+"\n")


## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="Collect seeds.")
    parser.add_argument("--sampath", type=str, default="sam/", help="Path for sam files (default=sam/).")
    parser.add_argument("--input", type=str, help="Input threshold2 filtered chiread file.")
    parser.add_argument("--binsize", type=int, default=10000,  help="Bin size (default=10kb).")
    parser.add_argument("--minqual", type=int, default=10, help="Minimum mapping quality (default = 10).")
    parser.add_argument("--insert", type=int, default=2000, help="Maximum estimated insert between forward and reverse (default=2000).")
    return parser.parse_args()

if __name__ in "__main__":
    args = parse_arguments()
    collect_chiread(args.sampath, args.input, args.binsize, args.insert, args.minqual)
