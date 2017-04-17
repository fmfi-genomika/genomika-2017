modDir="mods"
outDir="wigs"

rm -r $outDir
mkdir $outDir

date

for f in multiz7way/maf/*.maf;
do
	filename=$(basename "$f")
	chromName="${filename%.*}"
	printf "\n\n---------- phyloP: $chromName ----------\n"
	phyloP --method LRT --mode CONACC --wig-scores --chrom $chromName -i MAF "$modDir/$chromName.mod" $f > "$chromName.wigFix"
done

date
