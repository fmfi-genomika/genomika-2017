import sys

file = sys.argv[1]

def complement_nt(x):
    d = {'G': 'C', 'A': 'T', 'T': 'A', 'C': 'G', 'N': 'N'}
    return d[x]

def complement(seq):
    return ''.join(list(map(complement_nt, seq)))

repeat = 'GGGTTAGTCA'
repeat_c = complement(repeat)[::-1]

def find_tel(chrom, seq):
    result = []

    for rep, strand in [(repeat, '+'), (repeat_c, '-')]:
        i = 0
        while i < len(seq):
            next = 0
            try:
                next = seq[i:].index(rep)
                next += i
                result.append((next, strand))
            except ValueError:
                break
            i  = next + 1

    result.sort()
    with open('telomeric.bed', 'a') as out:
        for pos, s in result:
            out.write('%s %d %d %c%d 1 %c\n' % (chrom, pos, pos + len(repeat), s, pos+1, s))

with open(file, 'r') as inp:
    chrom = ''
    seq = ''
    for line in inp:
        if line[0] == '>':
            if chrom != '': find_tel(chrom, seq)
            chrom = line[1:].strip()
        else:
            seq += line.strip().upper()
    find_tel(chrom, seq)
