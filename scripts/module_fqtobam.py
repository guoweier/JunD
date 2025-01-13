# Weier Guo Python
# Create: 12/26/2024
# Introduction: Mapping reads from fastq to indexed sorted bam. 

## PACKAGES ##
import sys, math, os, time
import argparse
import gzip 

###### TO UPDATE: confirm necessary packages are all correctly loaded ########
###### TO UPDATE: choose mapping steps (whether to generate sam, or bam, or both, etc) ######
## usage: currently this scripts needs to be run in the folder with all fq files. 
#### UPDATE: to specify the location of fq files. #####

def rename_fqs():
    li = os.listdir(os.getcwd())
    fqs = list(filter(lambda x: x.endswith(".fq.gz") or x.endswith('.fastq.gz'), li))
    for fq in fqs:
        samplename = fq.split(".fq.gz")[0]
        print(samplename)
        with gzip.open(fq, "rt") as f, open(samplename+".rename.fq", "w") as o:
            for line in f:
                if line[0] == "@":
                    ln = line.split(" ")
                    x = '_'.join(ln)
                    o.write(x)
                else:
                    o.write(line)
        os.system(f"gzip {samplename}.rename.fq")

## CHECK REF NAME ##
def calibrate_chrname(database):
    # Regard the reference file are downloaded from NCBI, each scaffold name contains "chromosome". 
    if '/' in database:
        dbdir = database[:database.rindex('/')+1]
    else:
        dbdir = ''
    with open(database, "r") as ref, open(dbdir+'Ref_calibrated.fa', "w") as outfile:
        seq_id = 1
        for line in ref:
            if line.startswith(">"):
                if line.split("\n")[0] == ">Chr"+str(seq_id):
                    outfile.write(line)
                elif "chromosome" not in line:
                    break
                elif "chromosome" in line:
                    outfile.write(f">Chr{seq_id}\n")
                    seq_id += 1
            else:
                outfile.write(line)
    return dbdir 

## INDEX REF ##
def ref_indexing(database):
    if '/' in database:
        dbfiles = os.listdir(database[:database.rindex('/')+1])
        dname = database[database.rindex('/')+1:]
    else:
        dbfiles = os.listdir(os.getcwd())
        dname = database
    ind = filter(lambda x: dname in x, dbfiles)
    ref = ['', '.pac', '.ann', '.amb', '.bwt', '.sa']
    check = list((x in ref for x in map(lambda x: x[len(database):],ind)))
    if check.count('False') > 1 or len(check) < len(ref):
        print("bwa index "+database)
        os.system("bwa index "+database)

## MAKE RESULTS DIRECTORIES ##
def create_dict():
    li = os.listdir(os.getcwd())
    dirs = ['sam', 'bai', 'fq', 'sorted_bam']
    if False in list((x in li for x in dirs)):
        os.system("mkdir fq sam bam sorted_bam bai")
    return li

## BASIC MAPPING ##
def fqtobam(database, thread, file, mode):
    if mode == "i":
        name = file.split('.')[0]
        print(name)
        # bwa mem aligning
        print("bwa mem aligning")
        os.system("bwa mem -t "+str(thread)+" "+database+" "+file+" > "+name+"_aln.sam")
    elif mode == "u":
        name = file.split('_1.')[0]
        tail = file.split('_1.')[1]
        print(name)
        # bwa mem aligning
        print("bwa mem aligning")
        os.system("bwa mem -t "+str(threads)+" "+database+" "+name+"_1."+tail+" "+name+"_2."+tail+" > "+name+"_aln.sam")
    # sam to bam
    print("sam to bam")
    os.system("samtools view -bS "+name+"_aln.sam"+" > "+name+"_aln.bam")
    # bam to sorted_bam
    print("bam to sorted_bam")
    os.system("samtools sort "+name+"_aln.bam"+" > "+name+"_aln.sorted.bam")
    # remove bam
    print("remove unsorted bam")
    os.system("rm "+name+"_aln.bam")
    # bam indexing
    print("bam indexing")
    os.system("samtools index "+name+"_aln.sorted.bam"+" > "+name+"_aln.sorted.bam.bai")
    # put files into their directory
    os.system('mv *_aln.sam sam/')
    os.system("mv "+file+" fq/")
    os.system("mv *.sorted.bam sorted_bam/")
    os.system("mv *.bai bai/")

## RUN ##
def mapping(database, mode, thread):
    ref_indexing(database)
    li = create_dict()
    todo = list(filter(lambda x: x.endswith(".rename.fq.gz") or x.endswith(".rename.fastq.gz") or x.endswith(".rename.fq") or x.endswith(".rename.fastq"), li))
    for file in todo:
        fqtobam(database, thread, file, mode)

def arrange_orifq():
    os.system(f"mkdir orifq")
    orifqli = os.listdir(os.getcwd())
    orifqs = list(filter(lambda x: (x.endswith(".fq.gz") or x.endswith(".fastq.gz") or x.endswith(".fq") or x.endswith(".fastq")) and ("rename" not in x), orifqli))
    for fq in orifqs:
        os.system(f"mv {fq} orifq/")
    
## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="Mapping step: fastq to bam")
    parser.add_argument("--database", type=str, help="Input database file for mapping.")
    parser.add_argument("--mode", type = str, default='u', help="When mapping as uninterleaved, Read type, u = UNINTERLEAVED (default), i = INTERLEAVED")
    parser.add_argument("--thread", type=int, default=8, help="How many threads to use during alignment.")

    return parser.parse_args()

if __name__ in "__main__":
    args = parse_arguments()
    rename_fqs()
    arrange_orifq()
    dbdir = calibrate_chrname(args.database)
    database = dbdir+"Ref_calibrated.fa"
    mapping(database, args.mode, args.thread)
    



