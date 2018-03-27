import re
import sys

'''
Removes `chr` prefix appended by `gff3togtf.pl` to chromosomes
names in the first column, instead of changing `gff3togtf.pl` file.

python3 remove_chr_prefix input_file > output_file
'''

if (len(sys.argv) < 2):
    sys.exit('Please supply file name as the first argument.')

with open(sys.argv[1], 'r') as input_file:
    for line in input_file:
        cols = re.split(r'(\s+)', line.rstrip())
        if (len(cols)):
            cols[0] = re.sub(r'^chr(.*)$',  r'\1', cols[0])
        print(''.join(cols))

