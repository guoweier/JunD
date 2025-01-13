## PACKAGES ##
import sys, math, os, time
import argparse
import gzip 

def read_fq():
    li = os.listdir(os.getcwd())
    fqs = list(filter(lambda x: x.endswith(".fq.gz") or x.endswith('.fastq.gz'), li))
    for fq in fqs:
        samplename = fq.split(".fq.gz")[0]
        with gzip.open(fq, "rt") as f, open(samplename+"_new.fq", "w") as o:
            for line in f:
                if line[0] == "@":
                    ln = line.split(" ")
                    x = '_'.join(ln)
                    o.write(x)
                else:
                    o.write(line)
        os.system(f"gzip {samplename}_new.fq")

## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="Mapping step: fastq to bam")
    #parser.add_argument("--database", type=str, help="Input database file for mapping.")
    #parser.add_argument("--mode", type = str, default='u', help="When mapping as uninterleaved, Read type, u = UNINTERLEAVED (default), i = INTERLEAVED")
    #parser.add_argument("--thread", type=int, default=8, help="How many threads to use during alignment.")
    return parser.parse_args()

if __name__ in "__main__":
    args = parse_arguments()
    read_fq()
