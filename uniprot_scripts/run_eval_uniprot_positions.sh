#!/bin/bash

# xml_filename, mapping_filename, gene_positions_filename, output_filename
# yarLip
python3 eval_uniprot_positions.py 'uniprot-proteome%3AUP000001300.xml' mapping_yarLip.txt yarLip1_gene_positions.txt yarLip1_uniprot


# sacCer
python3 eval_uniprot_positions.py 'uniprot-proteome%3AUP000002311.xml' mapping_sacCer.txt sacCer3_gene_positions.txt sacCer3_uniprot