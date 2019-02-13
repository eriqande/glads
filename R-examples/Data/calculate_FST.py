# -*- coding: utf-8 -*-

# This is a python code which calculates FST, dxy and π using a sliding window 

import genomics
import gzip
import sys

filename = sys.argv[1]
pop1 = sys.argv[2]
pop2 = sys.argv[3]

#Define the input file and create a file handle:

genoFileName = filename
genoFile = gzip.open(genoFileName, "r")


#Provide names for the populations in our dataset, and list the individuals in each.
# We have 9 silvereye  populations

popNames = ["P1", "P2"]

samples = [["P1_" + str(i) for i in xrange(1,int(pop1))], ["P2_" + str(i) for i in xrange(1,int(pop2))]]

#Create a sampleData object, which stores this population information.

sampleData = genomics.SampleData(popInds = samples, popNames = popNames)

#Set the run parameters
windSize = 500000
stepSize = 500000
minSites = 1


#Now we create a windowGenerator object
windowGenerator = genomics.slidingWindows(genoFile, windSize, stepSize, skipDeepcopy = True)


#Define and open the output file.

outFileName = "FST_" + str(genoFileName) + ".csv"
outFile = open(outFileName, "w")


#Write the header line for the output file, which will include “pi” for each population, and “Fst” and “dxy” for each pair of populations.
outFile.write("CHR,BIN_START,BIN_END,N_VARIANTS")

for x in range(len(popNames)-1):
    outFile.write(",pi_" + popNames[x])
    for y in range(x+1,len(popNames)):
        outFile.write(",Fst_" + popNames[x] + "_" + popNames[y])
        outFile.write(",dxy_" + popNames[x] + "_" + popNames[y])

outFile.write("\n")


#Now we run the sliding window analysis. For each window we start by checking that there are enough sites in the window. If so, we make an Alignment object and calculate FST, dXY and π using the function genomics.popDiv. We collect the results as well as window information (scaffold, start, end and number of sites) and then write to the output file.

n=0
for window in windowGenerator:
    if window.seqLen() >= minSites:
        #if there are enough sites, make alignment object
        Aln = genomics.callsToAlignment(window.seqDict(), sampleData, seqType = "pairs")
        #get divergence stats
        statsDict = genomics.popDiv(Aln)
        stats = []
        for x in range(len(popNames)-1):
            #retrieve pi for each population
            stats.append(statsDict["pi_" + popNames[x]])
            for y in range(x+1,len(popNames)):
                #retrieve dxy and Fst for each pair
                stats.append(statsDict["fst_" + popNames[x] + "_" + popNames[y]])
                stats.append(statsDict["dxy_" + popNames[x] + "_" + popNames[y]])
        #add window stats and write to file
        out_data = [window.scaffold, window.start, window.end, window.seqLen()] + [round(s,3) for s in stats]
        outFile.write(",".join([str(x) for x in out_data]) + "\n")
    print n
    n+=1


#Close files and exit.

outFile.close()
genoFile.close()

exit()

