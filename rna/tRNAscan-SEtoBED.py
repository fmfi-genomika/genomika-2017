#! /usr/bin/python3
import sys

"""Script to convert tRNAscan-SE output to a BED file
appropriate for UCSC geneome browser.

tRNAscan-SE output should be on stdin,
BED file is on stdout.

Author: Dominik Bujna
"""


tRNA_type_number = {}
for line in sys.stdin:
    columns = line.split()
    #check if the line contains an entry
    if columns[1].isdigit():
        chromosomeName = columns[0]

        #start & end, smaller first, reposition to start at 0th position
        if int(columns[2]) < int(columns[3]):
            chromStart = int(columns[2]) - 1
            chromEnd =  int(columns[3])
            is_plus_strand = True
        else:
            chromStart = int(columns[3]) - 1
            chromEnd = int(columns[2])
            is_plus_strand = False

        #generate the name
        name = "tRNA-" + columns[4] + "-" + columns[5]
        if name not in tRNA_type_number.keys():
            tRNA_type_number[name] = 1
        else:
            tRNA_type_number[name] += 1
        name += str(tRNA_type_number[name])
        #recalculate the score value from 0-100 to 0-1000
        score = int(float(columns[8])*10)

        #check if there is an intron
        if int(columns[6]) == 0 and int(columns[7]) == 0:
            blockCount = 1
            blockStarts = "0"
            blockSizes = chromEnd - chromStart
        else:
        #tRNAscan-SE seems to give max 1 intron, hence 2 blocks if there is one
            blockCount = 2
        #correct for the strand direction and start at 0
            if is_plus_strand:
                intronStart = int(columns[6]) - 1
                exonStart = int(columns[7])
            else:
                intronStart = int(columns[7]) - 1
                exonStart = int(columns[6])
        #positions with chromStart as 0
            blockStarts = "0," + str(exonStart - chromStart)
            blockSizes = str(intronStart - chromStart) + "," + str(chromEnd - exonStart)

        #put together the processed information
        output = chromosomeName + "\t" + str(chromStart) + "\t" + str(chromEnd) + "\t" + name + "\t" + str(score) + "\t"
        if is_plus_strand:
            output += "+\t"
        else:
            output += "-\t"
        #thickly drawn part: empty (both ends at start); color black
        output += str(chromStart) + "\t" + str(chromStart) + "\t0,0,0\t"
        #add intron/exon info
        output += str(blockCount) + "\t" + str(blockSizes) + "\t" + str(blockStarts)
        print(output)
