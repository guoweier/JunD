# Weier Guo Python
# Create: 12/26/2024
# Introduction: Mapping reads from fastq to indexed sorted bam. 

## PACKAGES ##
import sys, math, os, time
import argparse

###### TO UPDATE: confirm necessary packages are all correctly loaded ########
###### TO UPDATE: choose mapping steps (whether to generate sam, or bam, or both, etc) ######
###### To Update: to automatic ref indexing #######

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
        os.system("mkdir sam bai fq sorted_bam")
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
def mapping(database, index, mode, thread):
    if index == "y":
        ref_indexing(database)
    li = create_dict()
    todo = filter(lambda x: x.endswith(".fq.gz") or x.endswith(".fastq.gz") or x.endswith(".fq") or x.endswith(".fastq"), li)
    for file in todo:
        fqtobam(database, thread, file, mode)
    
## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="Mapping step: fastq to bam")
    parser.add_argument("--database", type=str, required=True, help="Input database file for mapping.")
    ##### To Update: remove this parameter for automatic ref indexing #######
    parser.add_argument("--index", type=str, default ='n', help="whether do reference indexing step. n=NO (default), y=YES.")
    parser.add_argument("--mode", type = str, default='u', help="When mapping as uninterleaved, Read type, u = UNINTERLEAVED (default), i = INTERLEAVED")
    parser.add_argument("--thread", type=int, default=8, help="How many threads to use during alignment.")

    return parser.parse_args()

if __name__ in "__main__":
    args = parse_arguments()
    mapping(args.database, args.index, args.mode, args.thread)



