from os import path
import sys

input1_path = sys.argv[1]
input2_path = sys.argv[2]
output_path = path.splitext(input2_path)[0] + '.ensembl.bed'

def rnaToLocus(id):
    if not id in rnaToGene:
        return None
    if not rnaToGene[id] in geneToLocus:
        return None
    return geneToLocus[rnaToGene[id]]

geneToLocus = {}
rnaToGene = {}

with open(input1_path,'r') as input_file:
    #iterate through input gff file
    for line in input_file:
        #vyhod komenty
        if line[0] == '#':
            continue

        splitted_line = line.split('\t')
        #parse last info column
        splitted_info = splitted_line[8].split(';')
        if splitted_line[2] == 'gene':
            for i in splitted_info:
                if "locus_tag" in i:
                    geneToLocus[splitted_info[0].split('=')[1]] = i.split('=')[1]
                    break
        if splitted_line[2] == 'mRNA' or splitted_line[2] == 'tRNA':
            rnaToGene[splitted_info[0].split('=')[1]] = splitted_info[1].split('=')[1]

with open(input2_path,'r') as input_file, open(output_path,'w') as output_file:
    #iterate through input bed file
    for line in input_file:
        splitted_line = line.split('\t')
        if 'rna' in splitted_line[3] and rnaToLocus(splitted_line[3]) != None:
            locus = rnaToLocus(splitted_line[3]).rstrip('\n') + '_' + splitted_line[6].split(", ")[1]
        else:
            locus = splitted_line[6].split(", ")[1]
        locus = locus.rstrip('\n')
        output_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s"%(splitted_line[0],splitted_line[1],splitted_line[2],locus,splitted_line[4],splitted_line[5],splitted_line[6]))
