# Weier Guo Python
# Create: 12/19/2024

## PACKAGES ##
import sys, math, os, time
from optparse import OptionParser

# Please run this script in the folder with all fq files.
usage = "USAGE: python_fqtobam.py -d database_file.fa [-t #threads]"
parser = OptionParser(usage=usage)
parser.add_option("-d", "--database-file", dest="database", help="Input database file for mapping.")
##### To Update: remove this parameter for automatic ref indexing #######
parser.add_option("-i", "--indexing", dest="index", default ='n', help="whether do reference indexing step. n=NO (default), y=YES.")
parser.add_option("-m", "--mode", dest="mode", type = "str", default='u', help="When mapping as uninterleaved, Read type, u = UNINTERLEAVED (default), i = INTERLEAVED")
parser.add_option("-t", "--thread", dest="threads", default="8", help="How many threads to use during alignment.")
parser.add_option("-q", "--trimqual", dest="trimqual", default="20", help="Default mapping quality, for use as the bwa aln -q X during alignment.")
parser.add_option("-X", "--XoutModule", dest="modules", action="store_false", default = True, help="Loading modules for a unix module absed system, defualt = False ")
(opt, args) = parser.parse_args()

database = opt.database
threads = opt.threads
indexstep = opt.index

# load samtool module
#module("load samtools/1.17")

###### Update: to automatic ref indexing #######
## INDEX REF ##
if '/' in database:
   dbfiles = os.listdir(database[:database.rindex('/')+1])
   dname = database[database.rindex('/')+1:]
else:
   dbfiles = os.listdir(os.getcwd())
   dname = database
ind = filter(lambda x: dname in x, dbfiles)

if indexstep == "y":
    ref = ['', '.pac', '.ann', '.amb', '.bwt', '.sa']
    check = list((x in ref for x in map(lambda x: x[len(database):],ind)))
    if check.count('False') > 1 or len(check) < len(ref):
        print("bwa index "+database)
        os.system("bwa index "+database)

## MAKE RESULTS DIRECTORIES ##
li = os.listdir(os.getcwd())
dirs = ['sam', 'bai', 'fq', 'sorted_bam']
if False in list((x in li for x in dirs)):
   os.system("mkdir sam bai fq sorted_bam")

## MAPPING ##
if opt.mode == "i":
    todo = filter(lambda x: x[-6:] == '.fq.gz' or x.endswith("fastq.gz"), li)
    #todo.sort()
    for file in todo:
        name = file.split('.')[0]
        print(name)
        # bwa mem aligning
        print("bwa mem aligning")
        os.system("bwa mem -t "+threads+" "+database+" "+file+" > "+name+"_aln.sam")
        print("bwa mem -t "+threads+" "+database+" "+file+" > "+name+"_aln.sam")
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
elif opt.mode == "u":
    todo = filter(lambda x: x[-8:] == '_1.fq.gz' or x.endswith("_1.fastq.gz"), li)
    #todo.sort()
    for file in todo:
        name = file.split('_1.')[0]
        tail = file.split('_1.')[1]
        print(name)
        # bwa mem aligning
        print("bwa mem aligning")
        os.system("bwa mem -t "+threads+" "+database+" "+name+"_1."+tail+" "+name+"_2."+tail+" > "+name+"_aln.sam")
        print("bwa mem -t "+threads+" "+database+" "+name+"_1."+tail+" "+name+"_2."+tail+" > "+name+"_aln.sam")
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
        os.system("mv "+name+"_*."+tail+" fq/")
        os.system("mv *.sorted.bam sorted_bam/")
        os.system("mv *.bai bai/")
