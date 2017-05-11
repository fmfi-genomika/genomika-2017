from os import path
import sys

input1_path = sys.argv[1]
input2_path = sys.argv[2]
output_path = path.splitext(input2_path)[0] + '.ensembl.gtf'

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
        if splitted_line[2] == 'mRNA':
            rnaToGene[splitted_info[0].split('=')[1]] = splitted_info[1].split('=')[1]

with open(input2_path,'r') as input_file, open(output_path,'w') as output_file:
    #iterate through input gtf file
    for line in input_file:
        splitted_line = line.split('\t')
        if splitted_line[2] == "intron":
            continue
        locus = rnaToLocus(splitted_line[8].split(';')[0].split(' ')[1][1:-1])
        if locus != None:
            output_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t gene_id \"transcript:%s\"; transcript_id \"%s\"; PROTEIN_ID \"%s\";\n"%(splitted_line[0],splitted_line[1],splitted_line[2],splitted_line[3],splitted_line[4],splitted_line[5],splitted_line[6],splitted_line[7],locus,locus,locus))
