cd /data/sacCer3/similarOrganisms/alligned_data
sudo hgLoadMafSummary sacCer3 sacCerMultAlignLocal aln.maf
sudo cp aln.maf /gbdb/sacCer3/sacCerMultAlignLocal/
cd /gbdb/sacCer3/sacCerMultAlignLocal/
sudo hgLoadMaf sacCer3 sacCerMultAlignLocal
cd /kentsrc/kent/src/hg/makeDb/trackDb/
make alpha

