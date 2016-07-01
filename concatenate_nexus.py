## script tweaked from http://biopython.org/wiki/Concatenate_nexus by susana freitas, jul 2016
##
## script to concatenate several nexus files, adding the charset at the bottom of the file
## the output will serve as input for SVDQuartets
## it already has at the bottom the charset partitions for each marker

## a list of files have to be in the txt nex_files, present in the same folder as the nexus

from Bio.Nexus import Nexus
# the combine function takes a list of tuples [(name, nexus instance)...],
#if we provide the file names in a list we can use a list comprehension to
# create these tuples

## read all nexus files from a list separated by \n in a txt
with open('nex_files', 'r') as f:
    nexlist = [line.rstrip('\n') for line in f]

print(nexlist)

nexi =  [(fname, Nexus.Nexus(fname)) for fname in nexlist]
combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open('COMBINED.nex', 'w'))
