# Weier Guo Python 
# Create: 12/26/2024

## PACKAGES ##
import os, sys, math, time
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

parser.add_option("-c", "--controlfile", dest="f", default="NA", help="Input sam file.")
parser.add_option("-o", "--out", dest="o", help="Output bin file.")
parser.add_option("--unique", "-u", dest="unique", action="store_true", default = False, help="U/all only U (Off)")
parser.add_option("-m", "--maxsnps", dest="maxsnps", type = "int", default=5, help="Max number of SNPS (sam field 15)")
parser.add_option("-s", "--binsize", dest="binsize", type = "int", default=100000, help="Bin size")
parser.add_option("-b", "--breaks", dest="breaks", action="store_true", default = False, help="Insert breaks")
parser.add_option("-r", "--removefile", dest="r", default=False, help="A sam header file of reference sequences to ignore")
parser.add_option("-p", "--ploidy", dest="ploidy", type = "int", default=2, help="Ploidy multiplier, default is 2 for diploid")
parser.add_option("-C", "--covmode", dest="covmode", action="store_true", default = False, help="Only output coverage columns, not relative percent")
parser.add_option("-P", "--pair", dest="pair", action="store_true", default = False, help="For PE mapped reads take 1 good one")
(opt, args) = parser.parse_args()




