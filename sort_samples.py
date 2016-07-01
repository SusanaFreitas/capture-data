## script made by susana freitas, jul 2016

### script to sort sample names according to their species, so that they can be properly identified in the nex input for SVDQuartets
## script uses EggLib version 3.0.0b6 
import egglib
import os

print('###### sorting the samples in the fasta files and converting to nexus #####')
## to run it recursively in all *.fasta files, we will create the step function
def step(ext, dirname, names):
        ext = ext.lower()

    # this will make sure the extension is lowercase

        for name in names:
                if name.lower().endswith(ext):
                        my_align = egglib.io.from_fasta(name)
## index corresponds to the position of each sample in the list -1 (Python!!)
                        brau_list = [0, 47]
                        chlo_list = [1, 2]
                        clar_list = [3, 4]
                        defi_list = [5, 6, 39]
                        der_list = [7]
                        irano_list = [67, 68]
                        mix1_list = [9, 12, 13, 14]
                        mix2_list = [10, 11]
                        par_list = [17, 41]
                        par12_list = [16]
                        por1_list = [19, 21]
                        por2_list = [18, 22, 23]
                        praB_list = [25, 28]
                        praC_list = [24, 26, 27]
                        rad1_list = [64, 15, 20, 29, 30, 31, 32, 33, 35, 36, 37]
                        rad2_list = [34, 38, 40]
                        rud1_list = [45, 46, 54]
                        rud2_list = [42, 43, 44, 59]
                        sax_list = [8, 48]
                        stei_list = [49, 50, 51]
                        val_list = [52, 53, 55, 56, 65, 57, 66, 58, 60, 61, 62, 63]
## subset from my_align each species
                        brau = my_align.subset(brau_list)
                        chlo = my_align.subset(chlo_list)
                        clar = my_align.subset(clar_list)
                        defi = my_align.subset(defi_list)
                        der = my_align.subset(der_list)
                        irano = my_align.subset(irano_list)
                        mix1 = my_align.subset(mix1_list)
                        mix2 = my_align.subset(mix2_list)
                        par = my_align.subset(par_list)
                        par12 = my_align.subset(par12_list)
                        por1 = my_align.subset(por1_list)
                        por2 = my_align.subset(por2_list)
                        praB = my_align.subset(praB_list)
                        praC = my_align.subset(praC_list)
                        rad1 = my_align.subset(rad1_list)
                        rad2 = my_align.subset(rad2_list)
                        rud1 = my_align.subset(rud1_list)
                        rud2 = my_align.subset(rud2_list)
                        sax = my_align.subset(sax_list)
                        stei = my_align.subset(stei_list)
                        val = my_align.subset(val_list)
# to make a new file I will add on the brau samples the following
                        brau.add_samples(chlo)
                        brau.add_samples(clar)
                        brau.add_samples(defi)
                        brau.add_samples(der)
                        brau.add_samples(mix1)
                        brau.add_samples(mix2)
                        brau.add_samples(par)
                        brau.add_samples(par12)
                        brau.add_samples(por1)
                        brau.add_samples(por2)
                        brau.add_samples(praB)
                        brau.add_samples(praC)
                        brau.add_samples(rad1)
                        brau.add_samples(rad2)
                        brau.add_samples(rud1)
                        brau.add_samples(rud2)
                        brau.add_samples(sax)
                        brau.add_samples(stei)
                        brau.add_samples(val)
                        brau.add_samples(irano)
# now the brau subset has all sequences in the order I want
                        brau.to_fasta(fname=name + '.fa')
                        nex = brau.nexus()
                        nexfile = open(name + '.nexus', 'w')
                        nexfile.write(nex)

exten = '.fasta'
topdir = '.'
# Start the walk
os.path.walk(topdir, step, exten)

print("The job is done!")
