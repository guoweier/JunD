# Weier Guo Python
# Create: 01/05/2025
# INTRODUCTION: calculate pseudo-junction coverage, apply onto chiread file, filter samples. 
## We assume the breaks only occur on one haplotype. So by getting the pseudo-coverage, we have to divide by ploidy to get single haplotype coverage. 

## PACKAGES ##
import sys, math, os, time
import argparse 
import statistics

## GLOBAL PARAMETERS ##
covs = {}

def get_head(file):
    with open(file, "r") as f:
        head = f.readline()
        head = head.split("\n")[0]
        samples = head.split("\t")
    return samples 

def read_pseudo(pseudojun, sample):
    pseudo_head = get_head(pseudojun)
    samples = [s for s in pseudo_head if s in sample]
    spindex = [i for i, value in enumerate(pseudo_head) if value in samples]
    with open(pseudojun, "r") as pseudo:
        pshead = pseudo.readline()
        for line in pseudo:
            line = line.split("\n")[0]
            ln = line.split("\t")
            for i in spindex:
                if pseudo_head[i] not in covs:
                    covs[pseudo_head[i]] = []
                if int(ln[i]) != 0:
                    covs[pseudo_head[i]] += [int(ln[i])]
    return covs 

def get_coverage(pseudojun, sample, ploidy):
    covs = read_pseudo(pseudojun, sample)
    Scovs = {}
    for sample in covs:
        smode = statistics.mode(covs[sample])
        Scovs[sample] = smode/ploidy 
    return Scovs 

def filter_sample(input, pseudojun, sample, ploidy):
    Scovs = get_coverage(pseudojun, sample, ploidy)
    outdict = {}
    for s in list(Scovs.keys()):
        outdict[s] = []
    allsamples = get_head(input)
    aspindex = [[i,value] for i, value in enumerate(allsamples) if value in Scovs.keys()]
    with open(input, "r") as chiread:
        chihead = chiread.readline()
        for line in chiread:
            line = line.split("\n")[0]
            ln = line.split("\t")
            for item in aspindex:
                t1 = Scovs[item[1]]
                if float(ln[item[0]]) >= t1:
                    outdict[item[1]] += [ln]
    return outdict 

def get_ohead(input):
    with open(input, "r") as chiread:
        ohead = chiread.readline()
    return ohead 

def write_output(input, pseudojun, sample, ploidy):
    outdict = filter_sample(input, pseudojun, sample, ploidy)
    ohead = get_ohead(input)
    for key in outdict:
        with open("T2_"+key+".txt", "w") as o:
            o.write(ohead)
            for line in outdict[key]:
                o.write('\t'.join(line)+'\n')



## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="Threshold2: set samples threshold.")
    parser.add_argument("--input", type=str, help="Input threshold1 filtered chiread file.")
    parser.add_argument("--pseudojun", type=str, help="Input pseudo-junction file.")
    parser.add_argument("--ploidy", type=int, default=2, help="Input ploidy number (default=2).")
    parser.add_argument("--sample", type=str, nargs="+", help="Input sample names.")
    #parser.add_argument("--t2file", type=str, help="Output file with samples>=chosen threshold.")
    return parser.parse_args()

if __name__ in "__main__":
    args = parse_arguments()
    write_output(args.input, args.pseudojun, args.sample, args.ploidy)