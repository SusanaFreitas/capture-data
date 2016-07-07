## script to randomize the individuals selected per species (2)
# script uses EggLib version 3.0.0b6 

import egglib
import os

## to select random items from a list
import random

## to run it recursively in all *.fa files, we will create the step function
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

# to make a new file I will create a new alignment
                        #rand_align = egglib.Align.create(align)
                        rand_align = egglib.Align()

                        # group_of_items = {1, 2, 3, 4}  # a sequence or set will work here.
                        num_of_inds_per_sps = 2  # set the number to select here.
                        #list_of_random_items = random.sample(group_of_items, num_of_inds_per_sps)
                       #  first_random_item = list_of_random_items[0]
                       #  second_random_item = list_of_random_items[1]

                        rand_brau = random.sample(brau_list, num_of_inds_per_sps)
                        rand_chlo = random.sample(chlo_list, num_of_inds_per_sps)
                        rand_clar = random.sample(clar_list, num_of_inds_per_sps)
                        rand_defi = random.sample(defi_list, num_of_inds_per_sps)
                        rand_der = der_list
                        rand_irano = random.sample(irano_list, num_of_inds_per_sps)
                        rand_mix1 = random.sample(mix1_list, num_of_inds_per_sps)
                        rand_mix2 = random.sample(mix2_list, num_of_inds_per_sps)
                        rand_par = random.sample(par_list, num_of_inds_per_sps)
                        rand_par12 = par12_list
                        rand_por1 = random.sample(por1_list, num_of_inds_per_sps)
                        rand_por2 = random.sample(por2_list, num_of_inds_per_sps)
                        rand_praB = random.sample(praB_list, num_of_inds_per_sps)
                        rand_praC = random.sample(praC_list, num_of_inds_per_sps)
                        rand_rad1 = random.sample(rad1_list, num_of_inds_per_sps)
                        rand_rad2 = random.sample(rad2_list, num_of_inds_per_sps)
                        rand_rud1 = random.sample(rud1_list, num_of_inds_per_sps)
                        rand_rud2 = random.sample(rud2_list, num_of_inds_per_sps)
                        rand_sax = random.sample(sax_list, num_of_inds_per_sps)
                        rand_stei = random.sample(stei_list, num_of_inds_per_sps)
                        rand_val = random.sample(val_list, num_of_inds_per_sps)

                        brau = my_align.subset(rand_brau)
                        chlo = my_align.subset(rand_chlo)
                        clar = my_align.subset(rand_clar)
                        defi = my_align.subset(rand_defi)
                        der = my_align.subset(rand_der)
                        irano = my_align.subset(rand_irano)
                        mix1 = my_align.subset(rand_mix1)
                        mix2 = my_align.subset(rand_mix2)
                        par = my_align.subset(rand_par)
                        par12 = my_align.subset(rand_par12)
                        por1 = my_align.subset(rand_por1)
                        por2 = my_align.subset(rand_por2)
                        praB = my_align.subset(rand_praB)
                        praC = my_align.subset(rand_praC)
                        rad1 = my_align.subset(rand_rad1)
                        rad2 = my_align.subset(rand_rad2)
                        rud1 = my_align.subset(rand_rud1)
                        rud2 = my_align.subset(rand_rud2)
                        sax = my_align.subset(rand_sax)
                        stei = my_align.subset(rand_stei)
                        val = my_align.subset(rand_val)

# add on the alignment the random sequences
                        rand_align.add_samples(brau)
                        rand_align.add_samples(chlo)
                        rand_align.add_samples(clar)
                        rand_align.add_samples(defi)
                        rand_align.add_samples(der)
                        rand_align.add_samples(mix1)
                        rand_align.add_samples(mix2)
                        rand_align.add_samples(par)
                        rand_align.add_samples(par12)
                        rand_align.add_samples(por1)
                        rand_align.add_samples(por2)
                        rand_align.add_samples(praB)
                        rand_align.add_samples(praC)
                        rand_align.add_samples(rad1)
                        rand_align.add_samples(rad2)
                        rand_align.add_samples(rud1)
                        rand_align.add_samples(rud2)
                        rand_align.add_samples(sax)
                        rand_align.add_samples(stei)
                        rand_align.add_samples(val)
                        rand_align.add_samples(irano)
# now the brau subset has all sequences in the order I want
                        rand_align.to_fasta(fname=name + '.fa')
                        nex = rand_align.nexus()
                        nexfile = open(name + '.nexus', 'w')
                        nexfile.write(nex)

exten = '.fasta'
topdir = '.'
# Start the walk
os.path.walk(topdir, step, exten)

print("The job is done!")
