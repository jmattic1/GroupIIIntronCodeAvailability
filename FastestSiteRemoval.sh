	
	###IQ Tree paramter: -wsr
	####Where to put all the new trees
	outdir="/Users/jmattic1/Documents/EukandRiboprotein_ConcatTree/SubsetTree_SR4/"
	rm -rf $outdir
	mkdir $outdir
	###The rate file produced by -wsr
	treefile="/Users/jmattic1/Documents/EukandRiboprotein_ConcatTree/EukRibosomalProteins.concatenated.final.SR4.fas.rate"
	###The fasta file that built the tree
	staticfaafile="/Users/jmattic1/Documents/EukandRiboprotein_ConcatTree/EukRibosomalProteins.concatenated.final.SR4.fas"
	faafile=$(echo $staticfaafile | awk -F '/' -v outdir="$outdir" '{print outdir"/"$NF}')
	###Unwrap FASTA file
	cat $staticfaafile | awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' > $faafile
	cat $treefile | grep -v "#" | tail -n+2 > $outdir/RateFile.tsv
	sitenumber=$(cat $outdir/RateFile.tsv | wc -l)
	echo -e "10\t"$sitenumber > $outdir/RateSites.tsv

	##Run the loop to reduce sites and run the tree
	for i in $(seq 8 -1 1)
	do
		cat $outdir/RateFile.tsv | awk -v limit="$i" '$3 <= limit {print $0}' > $outdir/RateFile.$i.tsv
		sitenumber=$(cat $outdir/RateFile.$i.tsv | wc -l)
		echo -e $i"\t"$sitenumber >> $outdir/RateSites.tsv
		newfaafile=$(echo $faafile | awk -F '/' '{print $NF}' | sed -e "s/.faa/."$i".faa/g" | sed -e "s/.fas/."$i".fas/g")
		positions=$(cat $outdir/RateFile.$i.tsv | awk '{print $1}' | tr '\n' ',' | sed 's/,$//g')
		cat $faafile | awk -v positions="$positions" '{if($0 ~ /^>/) print $0; else {split($0, seqarray,""); split(positions, posarray,","); for(count = 1; count <= length(posarray); ++count){printf seqarray[posarray[count]];}; print ""}}' > $outdir/$newfaafile

		/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir/$newfaafile --alrt 1000 -B 1000

	done
