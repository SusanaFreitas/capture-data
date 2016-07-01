## made by susana freitas, Jun 2016 
##
### SCRIPT TO REMOVE ENDLINES FROM FASTAFILES
##
## usage:
# python remove_endlines.py filename
## from system (terminal) import the arguments
import sys
filename = sys.argv[1]

fasta = open(filename, 'r+')
lines = fasta.readlines()
with open(filename + ".fa", 'w') as newfasta:
    for line in lines:
        if line.startswith(">"):
            name = line
            newfasta.write('\n' +line + '')
            continue
        else:
            newfasta.write(line.replace("\n", "").replace('\r',''))


#atext = f.read().replace('\n',' ')
print("The job is done!")


