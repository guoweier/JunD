#! /usr/bin/env python

import os, sys, math, time
from optparse import OptionParser

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, Isabelle Henry, 2013
# This work is the property of UC Davis Genome Center - Comai Lab

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------

#This script outputs a read coverage by bin across a reference sequence, using a directory of samtools aligned .sam files as input. 
#It can also output a measure of relative coverage compared to a control dataset. There can be two types of control data: either a 
#control file is indicated or the mean of all files in the directory is calculated and used as the control set. In both cases, the 
#values for relative percentage per bin were calculated by dividing the percentage of reads mapping to that bin for the sample at 
#hand by the mean percentage of reads mapping to that bin for the control set. Finally, all values are multiplied by the ploidy 
#parameter (default 2) such that values for bins present in X copies would oscillate around X.
#
#This script also outputs a second small file containing the number of read processed from each sam file.
#
#Usage: [...] denotes optional parameters, if not indicated, default parameters are used.
#bin-by-sam.py -o output-bin-file.txt -s size-of-bins [-c control .sam file] [-u] [-m number of max snps, default is 5] [-b] [-r] [-p ploidy for relative percent calculation] [-C]
#
#For help
#bin-by-sam.py -h
#
#Input:
#Run in a directory with the input .sam files. If you want to use one of the files as control for the relative coverage, specify the file with the -c option.
#
#Parameters
#
#Required:
#-o, output file name
#-s, bin size (bps)
#
#Optional
#-c, use a control for relative percent coverage calculations, specify the file name here
#-u, use only samtools flagged unique reads (XT:A:U)
#-m, max snps from sam field 15  - default is 5
#-b, inserts empty lines between reference sequences in the result table for easier JMP parsing (do not use if the reference sequence does not contain a few major chromosomes or contigs)
#-r, remove file, a file in sam header format of reference sequences to ignore (example here)
#-p, ploidy, default is 2 (diploid), this is used as the multiplier in the relative coverage calculation
#-C, coverage only mode, this only outputs the read counts for each library, no relative coverage columns. This option cannot be used when a control library is specified
#
#Output:
#One file with a line per bin of each reference sequence and a column for each input .sam library, as well as the relative coverage per input .sam library.

usage = "\n\n%prog"
parser = OptionParser(usage=usage)
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


if opt.f != "NA" and opt.covmode == True:
   parser.error("Cannot specify a contol then supress contol relative coverage percent columns")

#takes in a sam header file of chroms/genes to ignore, must be specified by command line
remlist = []
remcount = {}
remsize = {}
if opt.r != False:
   rem = open(opt.r)
   while 1:
      x = rem.readline()
      if x == '':
         break
      if x[0] != '@':
         break
      if x[:3] == "@SQ":
         temp = (x[:-1].replace('SN:','').replace('LN:','').split('\t')[1:])
         key2 = temp[0]
         remlist.append(key2)
         remcount[key2] = 0
         remsize[key2] = int(temp[1])
   rem.close()

#read in list of sam files, must be in current directory and end in "_aln.sam"
li = os.listdir(os.getcwd())
todo = filter(lambda x: x.endswith('.sam'), li)
todo.sort()
f = open(todo[0])

#read sam header of chrom/genes to use
data = {}
all = []
sizes = []
lookup = {}
while 1:
   x = f.readline()
   if x[0] != '@':
      break
   if x[:3] == "@SQ":
      temp = (x[:-1].replace('SN:','').replace('LN:','').split('\t')[1:])
      key2 = temp[0]
      if key2 in ["ChrUn", "ChrSy"]:
         continue
      if key2 not in remlist:
         all.append(key2)
         sizes.append(int(temp[1]))
         lookup[key2] = int(temp[1])
         key1 = temp[0][-2:]
         if key1 not in data:
            data[key1] = {}
         if key2 not in data[key1]:
            data[key1][key2] = {}
         keys3 = range(int(temp[1])/int(opt.binsize)+1)
         for mod in keys3:
            if mod not in data[key1][key2]:
               data[key1][key2][mod] = {}

#for data export purpose; add blank lines based on size of largest reference
numblanks = max(sizes)/opt.binsize/10
      
f.close()
globalcount = {}
count = 0
liblist = []
if opt.pair == False:
   #per sam file, count reads
   for file in todo:
      f = open(file)
      print file
      libname = file.split('.')[0].replace('_aln','')
      liblist.append(libname)
      globalcount[libname] = 0
      count = 0
      #print time.asctime()
      for x in f:
         count +=1
         if count % 1000000 == 8:
            print count
   
         if x[0] == '@':
            continue       
         
         l = x[:-1].split('\t')   
         if l[2] == '*':
            continue
         if l[2] in ["ChrUn", "ChrSy"]:
            continue
         if l[1] == '0':
            pos = l[3]   
         elif l[1] == '16':   
            pos = int(l[3])+len(l[9])-1
         elif l[1] == '99':
            pos = l[3]
         elif l[1] == '163':
            pos = l[3]
         elif l[1] == '83':
            pos = l[3]
         elif l[1] == '147':
            pos = l[3]
         else:
            continue
    
    
         indexadd = 0
         if 'XC:i:' in x:
            indexadd = 1
            
#         if int(l[15+indexadd].split(':')[-1]) > opt.maxsnps:
#            continue
#         if opt.unique:
#            if l[11+indexadd].split(':')[-1] != 'U':
#               continue
         
         key1 = l[2][-2:]
         key2 = l[2]
         if key2 in remlist:
            remcount[key2] += 1
            continue
         key3 = int(pos) / opt.binsize
   
         
         try:
            data[key1][key2][key3][libname] += 1
         except:
            if libname not in data[key1][key2][key3]:
               data[key1][key2][key3][libname] = 1
         globalcount[libname] +=1
      f.close()
else:
   #per sam file, count reads
   for file in todo:
      f = open(file)
      print file
      libname = file.split('.')[0].replace('_aln','')
      liblist.append(libname)
      globalcount[libname] = 0
      count = 0
      #print time.asctime()
      while 1:
         x1 = f.readline()
         if x1 == '':
            break
         if x1[0] == '@':
            continue 
         x2 = f.readline()
         count +=1
         if count % 1000000 == 8:
            print count
   
         
         
         l1 = x1[:-1].split('\t')
         l2 = x2[:-1].split('\t')
         pic = 0
         if l1[1] not in ['0', '16'] and l2[1] not in ['0', '16']:
            continue
         elif l1[1] not in ['0', '16'] and l2[1] in ['0', '16']:
            pick = 2
         elif l1[1] in ['0', '16'] and l2[1] not in ['0', '16']:
            pick = 1
         elif l1[1] in ['0', '16'] and l2[1] in ['0', '16']:
            val1 = int(l1[5].split('M')[0].split("S")[-1])
            val2 = int(l2[5].split('M')[0].split("S")[-1])
            if val1 >= val2:
               pick = 1
            else:
               pick = 2
         else:
            x[878789]
            
         if pick == 1:      
            l = l1[:]
            x = x1[:]
         elif pick == 2:
            l = l2[:]
            x = x2[:]
         else:
            x[87654567]
           
         if l[2] == '*':
            continue
         if l[2] in ["ChrUn", "ChrSy"]:
            continue
         if l[1] == '0':
            pos = l[3]   
         elif l[1] == '16':   
            pos = int(l[3])+len(l[9])-1
         else:
            continue
    
    
         indexadd = 0
         if 'XC:i:' in x:
            indexadd = 1
            
         if int(l[15+indexadd].split(':')[-1]) > opt.maxsnps:
            continue
         if opt.unique:
            if l[11+indexadd].split(':')[-1] != 'U':
               continue
         
         key1 = l[2][-2:]
         key2 = l[2]
         if key2 in remlist:
            remcount[key2] += 1
            continue
         key3 = int(pos) / opt.binsize
   
         
         try:
            data[key1][key2][key3][libname] += 1
         except:
            if libname not in data[key1][key2][key3]:
               data[key1][key2][key3][libname] = 1
         globalcount[libname] +=1
      f.close()
 
#use control lib if specified, otherwise non applicable
control = opt.f
control = control.split('.')[0].replace('_aln','')

#create header for output file
header = ['Chrom', 'Strt', 'End']
header+= liblist
if opt.covmode == False:
   header+= map(lambda x: x+"/"+control, liblist)
o = open(opt.o, 'w')
o.write('\t'.join(header)+'\n')

#traverse each chromosome by bin
#then in each bin, first tabulate relative %
#then apply to create additional columns for relative%/control(or all)%
#when have one set for cov and other set for %, output line and repeat
#will output blank lines as/if speciifed
for chrom in all:
   part = chrom[-2:]
   bins = data[part][chrom].keys()
   bins.sort()
   for modbin in bins:
      libdata = data[part][chrom][modbin]             
      line = [chrom, modbin*opt.binsize+1, (modbin+1)*opt.binsize]
      if modbin == bins[-1]:
         line = [chrom, modbin*opt.binsize+1, lookup[chrom]]
      pers = {}
      for x in liblist:
         try:
            line.append(libdata[x])
            pers[x] = libdata[x]/float(globalcount[x])
         except:
            line.append(0)
            pers[x] = 0.0
      sums = sum(pers.values())/float(len(pers.values()))
      if opt.covmode == False:                  
         for x in liblist:
            if control == 'NA':
               try:
                  line.append(round(pers[x]/sums*opt.ploidy, 3))
               except ZeroDivisionError:
                  line.append(0.0)
            else:
               try:
                  line.append(round(pers[x]/pers[control]*opt.ploidy, 3))
               except ZeroDivisionError:
                  line.append('.')
      fline = map(lambda x: str(x), line)
      o.write('\t'.join(fline)+'\n')
   if opt.breaks == True:
      o.write(''.join(map(lambda x: '\n', range(numblanks))))

o.close()

#file to hold library statistics
#and output by each lib
o2 = open('readcounts-'+opt.o, 'w')

sechead = ['Lib', 'Reads', 'Reads/MB']
o2.write('\t'.join(sechead)+'\n')
tot = sum(sizes)/1000000.0
for x in liblist:
   o2.write('\t'.join([x, str(globalcount[x]), str(round(globalcount[x]/tot, 2))])+'\n')
if opt.r != False:
   o2.write('\n\nRemoved Reference Counts:\n')
   o2.write('Reference\tReads\n')
   for x in remlist:
      o2.write(x+'\t'+str(remcount[x])+'\n')

o2.close()           
         
         
      
   
   
   
   
   
   
      
      
