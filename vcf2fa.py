#### script to convert vcf file to fasta file
#### vcf file only includes SNPs (script does not work with indels) (yet)
## Can be used after the vcflib option vcfallelicprimitives that breaks the multialleles and indels into single snps

## usage >> three input files:

#### INPUT FILE 1
## freebayes vcf file: change filename in line 22

#### INPUT FILE 2
## text file with the names of all individuals ordered as in the vcf file: change filename in line 28

#### INPUT FILE 3
## fasta file with the reference sequences: change filename in line 37



print("################################################################################ \n"
      "################################ SU IS COOL! ################################### \n"
      "################################################################################")

#### INPUT FILE 1
vcf = open("test1000", "r")
# open the vcf file
lines = vcf.readlines()
# read the file line by line

#### INPUT FILE 2
indnames = open('indnames_var15.csv', 'r')
for line in indnames:
    indname = line.split(",")
#print indname[0]
from numpy import matrix

dicref = {}  # make a dictionary with the ref sequences names and

#### INPUT FILE 3
reffa = open('reference_probesgenome.fa', 'r')
seqfa = reffa.readlines()
for line in seqfa:
    if line.startswith('>'):
        seqname = line[1:]
    else:
        seq = line.rstrip() # rstrip() removes the the paragraph character which is the alst character of the line
        seqlen = len(seq)
        dicref[seqname.rstrip()] = seqlen
        # created a dictionary with the name of each marker and respective length
#print dicref
# make a dictionary for the redundant bases (to make consensus seq)
redundant = {'AG': 'R', 'GA': 'R', 'CT': 'Y', 'TC': 'Y', 'AC': 'M', 'CA': 'M', 'GC': 'S', 'CG': 'S', 'AT': 'W',
             'TA': 'W', 'GT': 'K', 'TG': 'K'}
# fmarker = []

# for a in range(len(array)):
#  fmarker[a] = open(array[a] + '.fa', 'w')

# now we will make a dictionary with the contig names and position of variant (key) +
# + all the info in the same line (value)
dici_matrix = {}  # this will make a dictionary for us to build on the positions+marker names and variant calls
for line in lines:
    if line.startswith('#'):
        continue
    else:
        tab = line.split()  # this corresponds to each column in each line
        name = tab[0]
        # previous_item = None # this will help me remembering the previous item in a loop
        pos = tab[1]
        keyname = name + "@&@" + str(pos)
        ref = tab[3]
        alt = tab[4]
        count = alt.count(",")
        #  print(ref) ## it works until here!
        ind = []  # make new array
        for i in range(9, 84, 1):  # individuals are from column 9 to col 84
            ind.append(tab[i])  # this will start reading ind 1 from column 9 until column 160 ind 151
        # print(tab[i]) until here is working!!!!!
        snp = []  # constructing a new array that will look like [A, G, N, ...] and will store the variants for each ind for each pos
        for i in range(0, 74, 1):
            nuc = ind[i].split(':')
            if count == 0:
                if nuc[0] == '0/0':
                    snp.append(ref)
                if nuc[0] == '1/0':
                    com = ref + alt
                    # print(com) # working!!!
                    snp.append(redundant.get(com))  # call the dictionary for heterozygous positions
                if nuc[0] == '0/1':
                    com = alt + ref
                    # print(com) # working!!!
                    snp.append(redundant.get(com))  # call the dictionary for heterozygous positions
                if nuc[0] == '1/1':
                    snp.append(alt)
                if nuc[0] == './.':
                    snp.append('N')
                if nuc[0] == '.':
                    snp.append('N')
                ## if it finds a phased genotype
                if nuc[0] == '0|0':
                    snp.append(ref)
                if nuc[0] == '1|0':
                    com = ref + alt
                    # print(com) # working!!!
                    snp.append(redundant.get(com))  # call the dictionary for heterozygous positions
                if nuc[0] == '0|1':
                    com = alt + ref
                    # print(com) # working!!!
                    snp.append(redundant.get(com))  # call the dictionary for heterozygous positions
                if nuc[0] == '1|1':
                    snp.append(alt)
                if nuc[0] == '.|.':
                    snp.append('N')
            elif count == 1:
                # for a in range(0,1,1):
                    all = alt.split(',')
                    #print(all[0])
                    #print(all[1])
                    if nuc[0] == '0/0':
                        snp.append(ref)
                    if nuc[0] == '1/0':
                        com = ref + all[0]
                        snp.append(redundant.get(com))  # call the dictionary for heterozygous positions
                    if nuc[0] == '0/1':
                        com = all[0] + ref
                        snp.append(redundant.get(com))  # call the dictionary for heterozygous positions
                    if nuc[0] == '1/1':
                        snp.append(all[0])
                    if nuc[0] == './.':
                        snp.append('N')
                    if nuc[0] == '0/2':
                        com = ref + all[1]
                        snp.append(redundant.get(com))
                    if nuc[0] == '1/2':
                        com = all[0] + all[1]
                        snp.append(redundant.get(com))
                    if nuc[0] == '2/0':
                        com = all[1] + ref
                        snp.append(redundant.get(com))
                    if nuc[0] == '2/1':
                        com = all[1] + all[0]
                        snp.append(redundant.get(com))
                    if nuc[0] == '2/2':
                        snp.append(all[1])
                    if nuc[0] == '.':
                        snp.append('N')
                    ## if it finds a phased genotype
                    if nuc[0] == '0|0':
                        snp.append(ref)
                    if nuc[0] == '1|0':
                        com = ref + all[0]
                        snp.append(redundant.get(com))  # call the dictionary for heterozygous positions
                    if nuc[0] == '0|1':
                        com = all[0] + ref
                        snp.append(redundant.get(com))  # call the dictionary for heterozygous positions
                    if nuc[0] == '1|1':
                        snp.append(all[0])
                    if nuc[0] == '.|.':
                        snp.append('N')
                    if nuc[0] == '0|2':
                        com = ref + all[1]
                        snp.append(redundant.get(com))
                    if nuc[0] == '1|2':
                        com = all[0] + all[1]
                        snp.append(redundant.get(com))
                    if nuc[0] == '2|0':
                        com = all[1] + ref
                        snp.append(redundant.get(com))
                    if nuc[0] == '2|1':
                        com = all[1] + all[0]
                        snp.append(redundant.get(com))
                    if nuc[0] == '2|2':
                        snp.append(all[1])
            elif count == 2:
                #for a in range(0,2,1):
                    all = alt.split(',')
                    if nuc[0] == '0/0':
                        snp.append(ref)
                    if nuc[0] == '1/0':
                        com = ref + all[0]
                        snp.append(redundant.get(com))  # call the dictionary for heterozygous positions
                    if nuc[0] == '0/1':
                        com = all[0] + ref
                        snp.append(redundant.get(com))  # call the dictionary for heterozygous positions
                    if nuc[0] == '1/1':
                        snp.append(all[0])
                    if nuc[0] == './.':
                        snp.append('N')
                    if nuc[0] == '0/2':
                        com = ref + all[1]
                        snp.append(redundant.get(com))
                    if nuc[0] == '1/2':
                        com = all[0] + all[1]
                        snp.append(redundant.get(com))
                    if nuc[0] == '2/0':
                        com = all[1] + ref
                        snp.append(redundant.get(com))
                    if nuc[0] == '2/1':
                        com = all[1] + all[0]
                        snp.append(redundant.get(com))
                    if nuc[0] == '2/2':
                        snp.append(all[1])
                    if nuc[0] == '0/3':
                        com = ref + all[2]
                        snp.append(redundant.get(com))
                    if nuc[0] == '1/3':
                        com = all[0] + all[2]
                        snp.append(redundant.get(com))
                    if nuc[0] == '2/3':
                        com = all[1] + all[2]
                        snp.append(redundant.get(com))
                    if nuc[0] == '3/0':
                        com = all[2] + ref
                        snp.append(redundant.get(com))
                    if nuc[0] == '3/1':
                        com = all[2] + all[0]
                        snp.append(redundant.get(com))
                    if nuc[0] == '3/2':
                        com = all[2] + all[1]
                        snp.append(redundant.get(com))
                    if nuc[0] == '3/3':
                        snp.append(all[2])
                    if nuc[0] == '.':
                        snp.append('N')
                    ## if it finds a phased genotype
                    if nuc[0] == '0|0':
                        snp.append(ref)
                    if nuc[0] == '1|0':
                        com = ref + all[0]
                        snp.append(redundant.get(com))  # call the dictionary for heterozygous positions
                    if nuc[0] == '0|1':
                        com = all[0] + ref
                        snp.append(redundant.get(com))  # call the dictionary for heterozygous positions
                    if nuc[0] == '1|1':
                        snp.append(all[0])
                    if nuc[0] == '.|.':
                        snp.append('N')
                    if nuc[0] == '0|2':
                        com = ref + all[1]
                        snp.append(redundant.get(com))
                    if nuc[0] == '1|2':
                        com = all[0] + all[1]
                        snp.append(redundant.get(com))
                    if nuc[0] == '2|0':
                        com = all[1] + ref
                        snp.append(redundant.get(com))
                    if nuc[0] == '2|1':
                        com = all[1] + all[0]
                        snp.append(redundant.get(com))
                    if nuc[0] == '2|2':
                        snp.append(all[1])
                    if nuc[0] == '0|3':
                        com = ref + all[2]
                        snp.append(redundant.get(com))
                    if nuc[0] == '1|3':
                        com = all[0] + all[2]
                        snp.append(redundant.get(com))
                    if nuc[0] == '2|3':
                        com = all[1] + all[2]
                        snp.append(redundant.get(com))
                    if nuc[0] == '3|0':
                        com = all[2] + ref
                        snp.append(redundant.get(com))
                    if nuc[0] == '3|1':
                        com = all[2] + all[0]
                        snp.append(redundant.get(com))
                    if nuc[0] == '3|2':
                        com = all[2] + all[1]
                        snp.append(redundant.get(com))
                    if nuc[0] == '3|3':
                        snp.append(all[2])
            dici_matrix[keyname] = snp
# print dici_matrix

import io
for key in dicref:  # iterate over the keys in the dictionary
    newkey = key.replace('/','_')
    with open(newkey + ".fa", 'w') as markerfile:
   #num += 1
        for n in range(74):
            markerfile.write('>' + indname[n] + '\n')
            for u in range(int(dicref.get(key))):
                keyname2 = key + "@&@" + str(u)
                if dici_matrix.has_key(keyname2):
                    snp2 = dici_matrix.get(keyname2)
                    markerfile.write(str(snp2[n]))
                else:
                    markerfile.write('N')
            markerfile.write('\n')

print "The job is done!"





