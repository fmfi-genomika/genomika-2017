from os import path
import sys

input_path = sys.argv[1]
output_path = path.splitext(input_path)[0] + '_chromosomes.gff'

with open(input_path,'r') as input_file, open(output_path,'w') as output_file:
    chromosome = None
    #iterate through input gff3 file
    for line in input_file:
        #parse comments
        if line[0] == '#':
            output_file.write(line)
            continue

        splitted_line = line.split('\t')
        #parse last info column
        splitted_info = splitted_line[8].split(';')
        #remember current chromosome
        if splitted_line[2] == 'region' and splitted_info[5] == 'genome=chromosome':
            chromosome = splitted_info[3].split('=')[1]
        #write line
        output_file.write("%s\t%s"%(chromosome,"\t".join(splitted_line[1:])))
