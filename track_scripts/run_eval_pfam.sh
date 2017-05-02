#!/bin/bash

#python3 eval_pfam_positions.py output_sacCer3_short.gff3 sacCer3_gene_positions.txt ##> sacCer3_pfam.bed
python3 eval_pfam_positions.py output_sacCer3.gff3 sacCer3_gene_positions.txt > sacCer3_pfam.bed
python3 eval_pfam_positions.py output_yarLip1.gff3 yarLip1_gene_positions.txt > yarLip1_pfam.bed