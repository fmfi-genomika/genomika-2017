outDir="mods"

rm -r $outDir
mkdir $outDir

tree="((((((sacCer3,sacPar),sacMik),sacKud),sacBay),sacCas),sacKlu)"

date

for f in multiz7way/maf/*.maf;
do
	filename=$(basename "$f")
	chromName="${filename%.*}"
	printf "\n\n---------- phyloFit: $chromName ----------\n"
	phyloFit --tree "$tree" --msa-format MAF --out-root "$outDir/$chromName" $f
done

date
