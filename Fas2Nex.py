## script to convert fasta format to nexus
## creates new nexus file from a list of files selected by their extension
# script uses EggLib version 3.0.0b6

import egglib
import os

## to run it recursively in all *.fa files, we will create the step function
def step(ext, dirname, names):
        ext = ext.lower()

    # this will make sure the extension is lowercase
        for name in names:
                if name.lower().endswith(ext):
                        my_align = egglib.io.from_fasta(name)
                        nex = my_align.nexus()
                        nexfile = open(name + '.nexus', 'w')
                        nexfile.write(nex)


exten = '.fa'
topdir = '.'
# Start the walk
os.path.walk(topdir, step, exten)

print("The job is done!")
