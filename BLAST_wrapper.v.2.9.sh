###Structure for reading inputs:
###Usage:
#  -blastDB <path to the blast database
#  -fasta <path to the query fasta sequence: can take a multi fasta
#  -eval <number> the e value cut off to apply to BLAST
#  -filter <number> the minimum percentage length of hits to keep
#  -type psiblast OR -type blastp the type of blast to use to search the database
#  -PSSM <path to PSSM file>. In the case of blastp, this has no effect and -fasta is required. In the case of PSIBLAST, this  d
#  replaces -fasta as the input. PSSM should be in the format genename-pssm
#  -outdir Determines the output directory for the run
#  -num_hits Determines how many matches BLAST will keep per gene per organism
#  -num_iters Determines how many times PSIBLAST will keep iterate the BLOSSUM matrix
#  -addorigin Determines whether to add the original fasta files to the aligning file
#  -percID applies a filter to the percent identity of each hit
#  -slenfilter applies a filter to the length of each subject hit, requiring it to be VALUE% of the query length
#  -family_threshold determines the num_hits family threshold to remove splice variants


PSSM=$(echo "notusingit")
num_hits=1
num_iters=6
addorigin=$(echo "yes")
trimal=$(echo "-gappyout")
reciprocal=$(echo "no")
filter=0
percID=0
blasteval=1e-5
type=psiblast
pfam_dom=$(echo "no")
pFAMDB=/Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
slenfilter=0
family_threshold=10
SR4=$(echo "no")
pfam_top=$(echo "no")

for ((var=1; var<=$#; var++))
do
	if [[ "${!var}" = "-blastDB" ]] 
		then
		var=$(($var+1))
		blastDB=${!var}
		echo "Variable blastDB is "$blastDB
	elif [[ "${!var}" = "-fasta" ]] 
		then
		var=$(($var+1))
		fasta=${!var}
		echo "Variable fasta is "$fasta
	elif [[ "${!var}" = "-eval" ]] 
		then
		var=$(($var+1))
		blasteval=${!var}
		echo "Variable blasteval is "$blasteval
	elif [[ "${!var}" = "-filter" ]] 
		then
		var=$(($var+1))
		filter=${!var}
		echo "Variable filter is "$filter
	elif [[ "${!var}" = "-type" ]] 
		then
		var=$(($var+1))
		type=${!var}
		echo "Variable type is "$type
	elif [[ "${!var}" = "-PSSM" ]] 
		then
		var=$(($var+1))
		PSSM=${!var}
		echo "Variable PSSM is "$PSSM
	elif [[ "${!var}" = "-num_hits" ]] 
		then
		var=$(($var+1))
		num_hits=${!var}
		echo "Variable num_hits is "$num_hits
	elif [[ "${!var}" = "-num_iters" ]] 
		then
		var=$(($var+1))
		num_iters=${!var}
		echo "Variable num_iters is "$num_iters
	elif [[ "${!var}" = "-addorigin" ]] 
		then
		var=$(($var+1))
		addorigin=${!var}
		echo "Variable addorigin is "$addorigin
	elif [[ "${!var}" = "-outdir" ]] 
		then
		var=$(($var+1))
		outdir=${!var}
		echo "Variable outdir is "$outdir
	elif [[ "${!var}" = "-trimal" ]] 
		then
		var=$(($var+1))
		trimal=${!var}
		echo "Variable trimal is "$trimal
	elif [[ "${!var}" = "-reciprocal" ]] 
		then
		var=$(($var+1))
		reciprocal=${!var}
		echo "Variable reciprocal is "$reciprocal
	elif [[ "${!var}" = "-pfam_dom" ]] 
		then
		var=$(($var+1))
		pfam_dom=${!var}
		echo "Variable pfam_dom is "$pfam_dom
	elif [[ "${!var}" = "-pfamDB" ]] 
		then
		var=$(($var+1))
		pfamDB=${!var}
		echo "Variable pfamDB is "$pfamDB
	elif [[ "${!var}" = "-pfam_top" ]] 
		then
		var=$(($var+1))
		pfam_top=${!var}
		echo "Variable pfam_top is "$pfam_top
	elif [[ "${!var}" = "-percID" ]] 
		then
		var=$(($var+1))
		percID=${!var}
		echo "Variable percID is "$percID
	elif [[ "${!var}" = "-percID" ]] 
		then
		var=$(($var+1))
		slenfilter=${!var}
		echo "Variable slenfilter is "$slenfilter
	elif [[ "${!var}" = "-family_threshold" ]] 
		then
		var=$(($var+1))
		family_threshold=${!var}
		echo "Variable family_threshold is "$family_threshold
	elif [[ "${!var}" = "-SR4" ]] 
		then
		var=$(($var+1))
		SR4=${!var}
		echo "Variable SR4 is "$SR4
	elif [[ "${!var}" = "-h" ]] 
		then
		echo -e "BLAST Pipeline usage: ./BLAST_wrapper -fasta <inputFASTA> -PSSM notusingit -blastDB <Input database> -eval 1e-5 -filter 0 -type psiblast -outdir <Path to outdir> -addorigin yes -num_hits 1 -num_iters 6 -trimal -gappyout -reciprocal no"
		echo "-fasta Input fasta. Required unless using a PSSM in psiblast mode"
		echo "-PSSM can substitute for the fasta input in psiblast mode"
		echo "-blastDB Input BLASTDB in blast ncbi format"
		echo "-eval The blast cutoff for blast hits"
		echo "-filter Minimum query coverage to be kept in blast hits [0-100]"
		echo "-type Determines which blast to use. Supports blastp and psiblast"
		echo "-outdir The directory where the output directory will be created and written to"
		echo "-addorigin Whether the query sequences should be put into the alignment file: yes/no"
		echo "-num_hits Number of hits to keep per sequence per species in the BLAST search (identical to blast num_hits)"
		echo "-num_iters Number of psiblast iterations to use. Psiblast mode only"
		echo "-trimal Trimal settings. Should have a - in front (ie, -gappyout or -strict)"
		echo "-reciprocal Include a reciprocal best hits filter. If using this, should be a trusted BLAST database of sequences."
		exit
	else
		echo "Unrecognized input, exiting"
		exit 1
	fi
done
rngnumber=$RANDOM
echo $rngnumber
###Use residues.pl to obtain the lengths of each fasta in the file
if [[ $PSSM == "notusingit" ]]
then
	genelength=$(perl /Applications/residues.pl $fasta | awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }')
	echo $genelength
	##Gene name as defined by the name of the fasta file
	genename=$(echo $fasta | awk -F "/" '{print $NF}' | sed 's/.faa//g')
	##Output directory defined from the gene name -- change /Users/jmattick/Documents/ to your own file path
	newoutdir=$(echo $fasta | awk -F "/" '{print $NF}' | sed 's/.faa//g' | sed 's/.fsa//g' | awk -v outdir="$outdir" -v rngnumber="$rngnumber" '{print outdir"/"$0"_BLAST_"rngnumber"_running"}')
else
	genelength=$(cat $PSSM | grep "length" | awk '{print $2}' | sed 's/,//g')
	genename=$(echo $PSSM | awk -F "/" '{print $NF}' | sed 's/-pssm.*//g')
	newoutdir=$(echo $PSSM | awk -F "/" '{print $NF}' | sed 's/-pssm.*//g' | awk -v outdir="$outdir" -v rngnumber="$rngnumber" '{print outdir"/"$0"_BLAST_"rngnumber"_running"}')

fi
##If this directory already exists, remove it and recreate it
rm -rf $newoutdir
mkdir $newoutdir
cd $newoutdir
ls $newoutdir
##If psiblast was specified, search with psiblast using the eval parameter as a cutoff
if [[ $type == "psiblast" ]]
then
	if [[ $PSSM == "notusingit" ]]
	then
		/Applications/miniconda3/bin/psiblast -query $fasta -db $blastDB -outfmt 6 -max_hsps $num_hits -num_threads 8 -out $newoutdir/$genename"."$blasteval".blast" -max_target_seqs 5129560 -evalue $blasteval -num_iterations $num_iters -qcov_hsp_perc $filter
	else
		/Applications/miniconda3/bin/psiblast -in_pssm $PSSM -db $blastDB -outfmt 6 -max_hsps $num_hits -num_threads 8 -out $newoutdir/$genename"."$blasteval".blast" -max_target_seqs 5129560 -evalue $blasteval -num_iterations $num_iters -qcov_hsp_perc $filter
	fi
else

	/Applications/miniconda3/bin/blastp -query $fasta -db $blastDB -outfmt 6 -max_hsps $num_hits -num_threads 8 -out $newoutdir/$genename"."$blasteval".blast" -max_target_seqs 5129560 -evalue $blasteval -qcov_hsp_perc $filter
fi
##Determine unique names of hits as defined by query_hit (to distinguish between hits from multiple query sequences)
paste -d '\t' <(cat $newoutdir/$genename"."$blasteval".blast" | sed '/^$/d' | grep -v "CONVERGED") <(cat $newoutdir/$genename"."$blasteval".blast" | sed '/^$/d' | grep -v "CONVERGED" | awk '{print $2}' | awk -F '_' '{print $1"_"$2}') <(cat $newoutdir/$genename"."$blasteval".blast" | sed '/^$/d' | grep -v "CONVERGED" | awk -v genelength="$genelength" '{print $4/genelength}') > $newoutdir/$genename"."$blasteval".species.blast"
cat $newoutdir/$genename"."$blasteval".species.blast" | awk -v percID="$percID" '$3 >= percID {print $0}' | awk '{print $0"\t"$1"_"$13}' > $newoutdir/$genename"."$blasteval".species.temp.blast"
mv $newoutdir/$genename"."$blasteval".species.temp.blast" $newoutdir/$genename"."$blasteval".species.blast"
#Get unique sequence ids without assessing length filter: output will be genename.blasteval.nolen.filtered.blast
#cat $newoutdir/$genename"."$blasteval".species.blast" | sort -k12,12nr | awk '!arr[$15]++' > $newoutdir/$genename"."$blasteval".filtered.blast"

cat $newoutdir/$genename"."$blasteval".species.blast" | awk '{print $15}' | sort | uniq > $newoutdir/$genename"."$blasteval".species.blast.uniqlist"
while read line
do
	###Editing time
	cat $newoutdir/$genename"."$blasteval".species.blast" | awk -v line="$line" 'line==$15 {print $0}' > $newoutdir/tempblast.eval.blast
	cat $newoutdir/$genename"."$blasteval".species.blast" | awk -v line="$line" 'line==$15 {print $0}' | awk '{print $2}' > $newoutdir/$line.list


	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $newoutdir/$line.list -outfmt %f -out $newoutdir/$line.faa
	cat $newoutdir/$line.faa | sed 's/\*//g' > $newoutdir/$line.nostar.faa
	numseq=$(cat $newoutdir/$line.faa | grep ">" | wc -l)
	echo $line
	if [[ $numseq -gt 1 && $num_hits -gt 1 ]]
	then
		/Applications/miniconda3/bin/mafft --quiet --anysymbol --auto $newoutdir/$line.nostar.faa > $newoutdir/$line.align.faa
		/Applications/miniconda3/bin/trimal -in $newoutdir/$line.align.faa -out $newoutdir/$line.align.trimmed.faa -automated1 -keepseqs
		/Applications/miniconda3/bin/distmat -sequence $newoutdir/$line.align.trimmed.faa -protmethod 2 -nucmethod 2 -outfile $newoutdir/$line.distmat
		cat $newoutdir/$line.distmat  | tail -n "$numseq" | sed 's/\t/,/g' | sed 's/ //g' | awk 'BEGIN{FS=OFS=","}{NF--; print}' | sed 's/^.*,0.00/0.00/g' | sed 's/nan/300/g' | sed 's/,$//g' > $newoutdir/$line.dist
		/usr/local/bin/Rscript /Applications/DistanceMatrix_condenser.R $newoutdir/$line.dist $numseq $family_threshold $newoutdir/tempblast.eval.blast $newoutdir

		cat $newoutdir/tempblast.eval.post.blast | awk -v num_hits="$num_hits" 'NR <= num_hits {print $0}' >> $newoutdir/$genename"."$blasteval".filtered.blast"
	else
		cat $newoutdir/tempblast.eval.blast | awk -v num_hits="$num_hits" 'NR <= num_hits {print $0}' >> $newoutdir/$genename"."$blasteval".filtered.blast"
	fi
	rm $newoutdir/tempblast*
	rm $newoutdir/$line.*
done < $newoutdir/$genename"."$blasteval".species.blast.uniqlist"

#Get unique sequence ids AND assessing length filter: output will be genename.eval.filtered.blast -- Replaced by  -qcov_hsp_perc $filter in the BLAST command
#cat $newoutdir/$genename"."$eval".species.blast" | sort -k12,12nr | awk '!arr[$15]++' | awk -v filter="$filter" '$14 > filter {print $0}' > $newoutdir/$genename"."$eval".filtered.blast"
##Get all unique protein hits in an output file called genename.blast.eval.list
cat $newoutdir/$genename"."$blasteval".filtered.blast" | awk '{print $2}' | sort | uniq > $newoutdir/$genename.blast.$blasteval.list
##Pull all unique protein hits out of the database into a multi fasta file
/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $newoutdir/$genename.blast.$blasteval.list -outfmt %f -out $newoutdir/$genename.blast.$blasteval.faa


####Reciprocal Step
if [[ $reciprocal != "no" ]]
then
	if [[ $type == "psiblast" ]]
	then
		if [[ $PSSM == "notusingit" ]]
		then
			/Applications/miniconda3/bin/psiblast -query $fasta -db $reciprocal -outfmt 6 -max_hsps 1 -num_threads 8 -out $newoutdir/$genename"."$blasteval".reciprocal.query.blast" -max_target_seqs 5129560 -evalue $blasteval -num_iterations $num_iters -qcov_hsp_perc $filter
			/Applications/miniconda3/bin/psiblast -query $newoutdir/$genename.blast.$blasteval.faa -db $reciprocal -outfmt 6 -max_hsps 1 -num_threads 8 -out $newoutdir/$genename"."$blasteval".reciprocal.hit.blast" -max_target_seqs 5129560 -evalue $blasteval -num_iterations $num_iters -qcov_hsp_perc $filter

		else
			/Applications/miniconda3/bin/psiblast -in_pssm $PSSM -db $reciprocal -outfmt 6 -max_hsps 1 -num_threads 8 -out $newoutdir/$genename"."$blasteval".reciprocal.query.blast" -max_target_seqs 5129560 -evalue $blasteval -num_iterations $num_iters -qcov_hsp_perc $filter
			/Applications/miniconda3/bin/psiblast -in_pssm $newoutdir/$genename.blast.$blasteval.faa -db $reciprocal -outfmt 6 -max_hsps 1 -num_threads 8 -out $newoutdir/$genename"."$blasteval".reciprocal.hit.blast" -max_target_seqs 5129560 -evalue $blasteval -num_iterations $num_iters -qcov_hsp_perc $filter

		fi
	else
		/Applications/miniconda3/bin/blastp -query $fasta -db $reciprocal -outfmt 6 -max_hsps 1 -num_threads 8 -out $newoutdir/$genename"."$blasteval".reciprocal.query.blast" -max_target_seqs 5129560 -evalue $blasteval -qcov_hsp_perc $filter
		/Applications/miniconda3/bin/blastp -query $newoutdir/$genename.blast.$blasteval.faa -db $reciprocal -outfmt 6 -max_hsps 1 -num_threads 8 -out $newoutdir/$genename"."$blasteval".reciprocal.hit.blast" -max_target_seqs 5129560 -evalue $blasteval -qcov_hsp_perc $filter
	fi
	##Determine unique names of hits as defined by query_hit (to distinguish between hits from multiple query sequences)
	paste -d '\t' <(cat $newoutdir/$genename"."$blasteval".reciprocal.query.blast" | sed '/^$/d' | grep -v "CONVERGED") <(cat $newoutdir/$genename"."$blasteval".reciprocal.query.blast" | sed '/^$/d' | grep -v "CONVERGED" | awk '{print $2}' | awk -F '_' '{print $1"_"$2}') <(cat $newoutdir/$genename"."$blasteval".reciprocal.query.blast" | sed '/^$/d' | grep -v "CONVERGED" | awk -v genelength="$genelength" '{print $4/genelength}') > $newoutdir/$genename"."$blasteval".reciprocal.query.species.blast"
	cat $newoutdir/$genename"."$blasteval".reciprocal.query.species.blast" | awk '{print $0"\t"$1"_"$13}' > $newoutdir/$genename"."$blasteval".reciprocal.query.species.temp.blast"
	mv $newoutdir/$genename"."$blasteval".reciprocal.query.species.temp.blast" $newoutdir/$genename"."$blasteval".reciprocal.query.species.blast"
	paste -d '\t' <(cat $newoutdir/$genename"."$blasteval".reciprocal.hit.blast" | sed '/^$/d' | grep -v "CONVERGED") <(cat $newoutdir/$genename"."$blasteval".reciprocal.hit.blast" | sed '/^$/d' | grep -v "CONVERGED" | awk '{print $2}' | awk -F '_' '{print $1"_"$2}') <(cat $newoutdir/$genename"."$blasteval".reciprocal.hit.blast" | sed '/^$/d' | grep -v "CONVERGED" | awk -v genelength="$genelength" '{print $4/genelength}') > $newoutdir/$genename"."$blasteval".reciprocal.hit.species.blast"
	cat $newoutdir/$genename"."$blasteval".reciprocal.hit.species.blast" | awk '{print $0"\t"$1"_"$13}' > $newoutdir/$genename"."$blasteval".reciprocal.hit.species.temp.blast"
	mv $newoutdir/$genename"."$blasteval".reciprocal.hit.species.temp.blast" $newoutdir/$genename"."$blasteval".reciprocal.hit.species.blast"

	#Get unique sequence ids from the query and hit blasts -- take only proteins from the hits from the query
	cat $newoutdir/$genename"."$blasteval".reciprocal.query.species.blast" | sort -k12,12nr | awk '!arr[$15]++' > $newoutdir/$genename"."$blasteval".reciprocal.query.nolen.filtered.blast"
	cat $newoutdir/$genename"."$blasteval".reciprocal.hit.species.blast" | sort -k12,12nr | awk '!arr[$15]++' > $newoutdir/$genename"."$blasteval".reciprocal.hit.nolen.filtered.blast"

	cat $newoutdir/$genename"."$blasteval".reciprocal.query.nolen.filtered.blast" | awk '{print $2}' | sort | uniq  > $newoutdir/$genename".blast."$blasteval".query.list"

	while read line
	do
		cat $newoutdir/$genename"."$blasteval".reciprocal.hit.nolen.filtered.blast" | awk -v line="$line" '$2 == line {print $0}' >> $newoutdir/$genename"."$blasteval".reciprocalonly.hit.nolen.filtered.blast"
	done < $newoutdir/$genename".blast."$blasteval".query.list"
	
	cat $newoutdir/$genename"."$blasteval".reciprocalonly.hit.nolen.filtered.blast" | awk '{print $1}' | sort | uniq > $newoutdir/$genename.blast.reciprocal.$blasteval.list
	##Pull all unique protein hits out of the database into a multi fasta file
	mv $newoutdir/$genename.blast.$blasteval.faa $newoutdir/$genename.blast.$blasteval.nonreciprocal.faa
	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $newoutdir/$genename.blast.reciprocal.$blasteval.list -outfmt %f -out $newoutdir/$genename.blast.$blasteval.faa

fi

if [[ $pfam_dom != "no" ]]
then
	/Applications/miniconda3/bin/hmmscan -E 1e-5 -o $newoutdir/$genename.blast.$blasteval.hmm.out --tblout $newoutdir/$genename.blast.$blasteval.hmm.tsv --domtblout $newoutdir/$genename.blast.$blasteval.dom.tsv --cpu 32 $pfamDB $newoutdir/$genename.blast.$blasteval.faa
	if [[ $pfam_top != "no" ]]
	then
		cat $newoutdir/$genename.blast.$blasteval.hmm.tsv | awk '{print $1"\t"$3"\t"$6}' > $newoutdir/$genename.blast.$blasteval.hmm.tab.tsv
		awk '{if ($2 in container && $3 > max[$2]) {max[$2] = $3;line[$2] = $0;container[$2] = $2;} else if (!($2 in container)) { max[$2] = $3;line[$2] = $0;container[$2] = $2;}}END { for (i in container) {print line[i];}}' $newoutdir/$genename.blast.$blasteval.hmm.tab.tsv | grep -E $pfam_dom | awk '{print $2}' > $newoutdir/$genename.blast.$blasteval.pfam.$pfam_dom.list
		/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $newoutdir/$genename.blast.$blasteval.pfam.$pfam_dom.list -outfmt %f -out $newoutdir/$genename.blast.$blasteval.faa

	else
		cat $newoutdir/$genename.blast.$blasteval.hmm.tsv | grep -E $pfam_dom | awk '{print $3}' | sort | uniq > $newoutdir/$genename.blast.reciprocal.$blasteval.pfam.$pfam_dom.list
		mv $newoutdir/$genename.blast.$blasteval.faa $newoutdir/$genename.blast.$blasteval.nonpfam.faa
	fi
	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $newoutdir/$genename.blast.reciprocal.$blasteval.pfam.$pfam_dom.list -outfmt %f -out $newoutdir/$genename.blast.$blasteval.faa
	
fi

mv $newoutdir/$genename.blast.$blasteval.faa $newoutdir/$genename.blast.$blasteval.unfiltered.faa
perl /Applications/residues.pl $newoutdir/$genename.blast.$blasteval.unfiltered.faa | awk -v slenfilter="$slenfilter" -v genelength="$genelength" '$2 > slenfilter*genelength {print $1}' > $newoutdir/$genename.blast.slenfilter.$blasteval.list
/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $newoutdir/$genename.blast.slenfilter.$blasteval.list -outfmt %f -out $newoutdir/$genename.blast.$blasteval.faa



##Add the original query sequences to the multi fasta
if [[ $addorigin == "yes" ]]
then
	cat $newoutdir/$genename.blast.$blasteval.faa $fasta > $newoutdir/$genename.blast.$blasteval.origin.faa
else
	cat $newoutdir/$genename.blast.$blasteval.faa > $newoutdir/$genename.blast.$blasteval.origin.faa
fi

cat $newoutdir/$genename.blast.$blasteval.origin.faa | sed 's/\*//g' > $newoutdir/$genename.blast.$blasteval.origin.nostar.faa
##Align the multi fasta using MAFFT run with automatically determined parameters
/Applications/miniconda3/bin/mafft --anysymbol --reorder --auto $newoutdir/$genename.blast.$blasteval.origin.nostar.faa > $newoutdir/$genename.blast.$blasteval.align.faa
##Trimal to trim the alignment: gappyout parameter errs on the side of keeping gaps
/Applications/miniconda3/bin/trimal -in $newoutdir/$genename.blast.$blasteval.align.faa -out $newoutdir/$genename.blast.$blasteval.align.trimal.faa $trimal
##Build a tree using IQTREE, number of cores is automatically determined as well as appropriate substitution model

if [[ $SR4 == "yes" ]]
then
	/Applications/aa_recoder.py -i $newoutdir/$genename.blast.$blasteval.align.trimal.faa -re SR4 -o $newoutdir/$genename.blast.$blasteval.align.trimal.SR4
	/Applications/miniconda3/bin/iqtree -T AUTO -s $newoutdir/$genename.blast.$blasteval.align.trimal.SR4.fas --alrt 1000 -B 1000 -safe
else
	/Applications/miniconda3/bin/iqtree -T AUTO -s $newoutdir/$genename.blast.$blasteval.align.trimal.faa --alrt 1000 -B 1000 -safe
fi

#IQTree1
##/Users/jmattick/Downloads/iqtree-1.6.12-MacOSX/bin/iqtree -nt AUTO -s $outdir/$genename.blast.$eval.align.trimal.faa
#IQTree2

maxlength=$(perl /Applications/residues.pl $newoutdir/$genename.blast.$blasteval.align.trimal.faa | sort -k2,2nr | awk '{print $2}' | head -1)
numseqs=$(cat $newoutdir/$genename.blast.$blasteval.align.trimal.faa | grep ">" | wc -l | sed 's/ //g')

finaloutdir=$(echo $genename | awk -v outdir="$outdir" -v blasteval="$blasteval" -v filter="$filter" -v num_hits="$num_hits" -v num_iters="$num_iters" -v maxlength="$maxlength" -v numseqs="$numseqs" -v rngnumber="$rngnumber" '{print outdir"/"$0"_e"blasteval"f"filter"nh"num_hits"ni"num_iters"_"numseqs"x"maxlength"_"rngnumber}')

echo $finaloutdir
cd ..
mv $newoutdir $finaloutdir
