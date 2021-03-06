modDir="mods"
outDir="wigs"

rm -r $outDir
mkdir $outDir

date

for f in maf/*.maf;
do
	filename=$(basename "$f")
	chromName="${filename%.*}"
	printf "\n\n---------- phyloP: $chromName ----------\n"
	phyloP --method LRT --mode CONACC --wig-scores --chrom $chromName -i MAF "$modDir/$chromName.mod" $f > "$outDir/$chromName.wigFix"
done

date
