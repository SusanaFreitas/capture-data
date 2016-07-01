import egglib
import os
print("###########################################################\n"
      "## extracting columns in fasta files with less than 7 N  ##\n"
      "#######   (corresponding to 10% of missing data)   ########\n"
      "###########################################################\n"
      "\n")

# /home/sanufas/Documents/work/Capture/PYTHON/fb22_filtered
## to run it recursively in all *.fa files, we will create the step function
def step(ext, dirname, names):
    # this will make sure the extension is lowercase
    ext = ext.lower()

    for name in names:
        if name.lower().endswith(ext):
            print(os.path.join(dirname,name))
            my_align = egglib.io.from_fasta(name)
           # my_align1 = my_align2.del_sample(3)
           # my_align = my_align1.del_sample(3)
            length=my_align.ls
            #print(length)
            valid_sites=[]
            #my_align.column(1)
            for i in range(0,int(length),1):
            ## the column returns a list of integers corresponding the the alleles for all samples
                listno = (my_align.column(i, outgroup=False))
                col = map(chr,listno)
            # map(chr, int(my_align.column(i)))
            # col=my_align.column(i)
            # print(col)
                count=col.count('N')
            # print(count)
                if count <= 7:
                   valid_sites.append(i)
            # print(valid_sites)
                new_align=my_align.extract(valid_sites)
                new_align.to_fasta(fname=name + "sta")
                nex = new_align.nexus()
                nexfile = open(name + ".nex", "w")
                nexfile.write(nex)
                phy = new_align.phylip()
                phyfile = open(name + ".phy", "w")
                phyfile.write(phy)
                phyml = new_align.phyml()
                phymlfile = open(name + ".phyml", "w")
                phymlfile.write(phyml)


exten = '.fa'
topdir = '.'
# Start the walk
os.path.walk(topdir, step, exten)

print("The job is done!")

