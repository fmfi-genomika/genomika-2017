outDir="mods"

rm -r $outDir
mkdir $outDir

tree="$1"

date

for f in maf/*.maf;
do
	filename=$(basename "$f")
	chromName="${filename%.*}"
	printf "\n\n---------- phyloFit: $chromName ----------\n"
	phyloFit --tree "$tree" --msa-format MAF --out-root "$outDir/$chromName" $f
done

date
