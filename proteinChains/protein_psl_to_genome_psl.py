from sys import argv
from sys import exit

def load_positions(filename):
    fpos = open(filename, "r")
    
    positions = {}
    line = fpos.readline()
    while True:
        line = fpos.readline()
        if line == "":
            break
        cols = line.strip().split()
        exons = []
        for start, end in zip(cols[9].split(","), cols[10].split(",")):
            if start == "" or end == "":
                continue
            start = int(start)
            end = int(end)
            exons.append((start, end))
        positions[cols[1]] = (cols[2], cols[3], exons)
    return positions

def gene_pos_to_genome_pos(gene_info, pos):
    # DONE: - strand
    if gene_info[1] == "+":
        for exon in gene_info[2]:
            exon_length = exon[1] - exon[0]
            if pos < exon_length:
                return exon[0] + pos
            
            pos -= exon_length
    else:
        for exon in reversed(gene_info[2]):
            exon_length = exon[1] - exon[0]
            if pos < exon_length:
                return exon[1] - pos - 1
            
            pos -= exon_length

    print("nieco sa dodzigalo v gene_pos_to_genome_pos")

def get_all_genome_positions(gene_info, starts, lengths):
    all_positions = []
    for start, length in zip(starts, lengths):
        for protein_pos in range(start, start + length):
            for i in range(3):
                all_positions.append(gene_pos_to_genome_pos(gene_info, protein_pos * 3 + i))
    return all_positions

def get_blocks(positions1, positions2):
    genome_positions1 = []
    genome_positions2 = []
    lengths = []
    current_length = 0
    prev_pos1 = -100
    prev_pos2 = -100

    for pos1, pos2 in zip(positions1, positions2):
        if prev_pos1 + 1 != pos1 or prev_pos2 + 1 != pos2:
            genome_positions1.append(pos1)
            genome_positions2.append(pos2)
            lengths.append(current_length)
            current_length = 0

        current_length += 1
        
        prev_pos1 = pos1
        prev_pos2 = pos2
    lengths = lengths[1:]
    lengths.append(current_length)
    
    return (lengths, genome_positions1, genome_positions2)

def load_chrom_lengths(filename):
    flengths = open(filename, "r")
    lengths = {}
    flengths.readline()
    for line in flengths:
        cols = line.strip().split()
        lengths[cols[0]] = int(cols[1])
    return lengths

def is_ascending(values):
    return values[0] < values[1]

if len(argv) != 6:
    usage = """
    Usage:
    protein_psl_to_genome_psl.py <target gene positions> <query gene positions> <target chromosome lengths> <query chromosome lengths> <protein alighnemts>

    This program creates alignments from genome to genome based on the alignments of proteins. The bases of codons of aligned aminoacids are marked as aligned even if they are different.
    
    You can extract the gene positions from genome browser using command:
    hgsql <genome db> -B -e "SELECT * FROM ensGene"

    To get the chromosome lengths use:
    hgsql <genome db> -B -e "SELECT chrom, size FROM chromInfo"

    """
    print (usage)
    exit(0)

target_gene_positions = load_positions(argv[1])
query_gene_positions = load_positions(argv[2])

target_chrom_lengths = load_chrom_lengths(argv[3])
query_chrom_lengths = load_chrom_lengths(argv[4])

protein_filename = argv[5]

protein_file = open(protein_filename, "r")
while True:
    line = protein_file.readline()
    if line == "":
        break
    cols = line.strip().split()
    target = target_gene_positions[cols[13]]
    query = query_gene_positions[cols[9]]

    lengths = [int(x) for x in cols[18].split(",") if x != ""]
    target_block_starts = [int(x) for x in cols[20].split(",") if x != ""]
    query_block_starts = [int(x) for x in cols[19].split(",") if x != ""]

    # DONE: remove
    # if target[1] == "-" or query[1] == "-":
    #     continue
    
    target_genome_positions = get_all_genome_positions(target, target_block_starts, lengths)
    query_genome_positions = get_all_genome_positions(query, query_block_starts, lengths)

    if not is_ascending(target_genome_positions):
        target_genome_positions.reverse()
        query_genome_positions.reverse()

    qstart = min(query_genome_positions)
    qend = max(query_genome_positions) + 1
    tstart = min(target_genome_positions)
    tend = max(target_genome_positions) + 1
    
    if not is_ascending(query_genome_positions):
        query_genome_positions = [query_chrom_lengths[query[0]] - pos for pos in query_genome_positions]
    

    genome_block_lengths, target_genome_block_starts, query_genome_block_starts = get_blocks(target_genome_positions, query_genome_positions)

    out_cols = []

    #vypisy
    out_cols.append(str(int(cols[0]) * 3))
    out_cols.append(str(int(cols[1]) * 3))
    out_cols.append(str(int(cols[2]) * 3))
    out_cols.append(str(int(cols[3]) * 3))

    # TODO: replace with exact numbers
    out_cols.append(str(int(cols[4]) ))
    out_cols.append(str(int(cols[5]) * 3))
    out_cols.append(str(int(cols[6]) ))
    out_cols.append(str(int(cols[7]) * 3))

    out_cols.append("+" if query[1] == target[1] else "-")

    out_cols.append(query[0])
    out_cols.append(str(query_chrom_lengths[query[0]]))
    out_cols.append(str(qstart))
    out_cols.append(str(qend))

    out_cols.append(target[0])
    out_cols.append(str(target_chrom_lengths[target[0]]))
    out_cols.append(str(tstart))
    out_cols.append(str(tend))

    out_cols.append(str(len(genome_block_lengths)))

    out_cols.append(",".join([str(x) for x in genome_block_lengths]) + ",")
    out_cols.append(",".join([str(x) for x in query_genome_block_starts]) + ",")
    out_cols.append(",".join([str(x) for x in target_genome_block_starts]) + ",")


    print("\t".join(out_cols))


    