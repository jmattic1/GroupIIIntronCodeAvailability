
if [[ 1 -eq 2 ]]
then
	###Running the BLAST searches
	rm -rf /Users/jmattic1/Documents/TERT_reanalysis/
	mkdir /Users/jmattic1/Documents/TERT_reanalysis/

	mkdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout/
	###reciprocal
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/multi_RDP.faa -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout/ -addorigin no -reciprocal /Applications/Databases/ReciprocalDB_archbact/ReciprocalStandard.protein.db.faa -trimal -gappyout
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Introns.faa -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout/ -addorigin no -reciprocal /Applications/Databases/ReciprocalDB_archbact/ReciprocalStandard.protein.db.faa -trimal -gappyout
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/LTR.faa -blastDB /Applications/Databases//EukArchBactViral/EukArchBactViral.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout/ -addorigin no -reciprocal /Applications/Databases/ReciprocalDB_archbact/ReciprocalStandard.protein.db.faa -trimal -gappyout
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/Spliceosome.faa -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout/ -addorigin no -reciprocal /Applications/Databases/ReciprocalDB_archbact/ReciprocalStandard.protein.db.faa -trimal -gappyout
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/TERT.faa -blastDB /Applications/Databases//EukArchBactViral/EukArchBactViral.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout/ -addorigin no -reciprocal /Applications/Databases/ReciprocalDB_archbact/ReciprocalStandard.protein.db.faa -trimal -gappyout

	#####Single taxa gappyout


	mkdir /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout/
	###reciprocal
	samtools faidx /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Introns.faa G2I_Arabidopsis_t > /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Intron_single.faa
	samtools faidx /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/LTR.faa LTR_Daldinia_c > /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/LTR_single.faa
	samtools faidx /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/Spliceosome.faa Splice_Homo_s > /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/Spliceosome_single.faa
	samtools faidx /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/TERT.faa TERT_Arabidopsis_t > /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/TERT_single.faa

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Intron_single.faa -blastDB /Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout/ -addorigin no -reciprocal /Applications/Databases/ReciprocalDB_archbact/ReciprocalStandard.protein.db.faa -trimal -gappyout
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/LTR_single.faa -blastDB /Applications/Databases//EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout/ -addorigin no -reciprocal /Applications/Databases/ReciprocalDB_archbact/ReciprocalStandard.protein.db.faa -trimal -gappyout
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/Spliceosome_single.faa -blastDB /Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout/ -addorigin no -reciprocal /Applications/Databases/ReciprocalDB_archbact/ReciprocalStandard.protein.db.faa -trimal -gappyout
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/TERT_single.faa -blastDB /Applications/Databases//EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout/ -addorigin no -reciprocal /Applications/Databases/ReciprocalDB_archbact/ReciprocalStandard.protein.db.faa -trimal -gappyout

	mkdir /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout_norecip/
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Intron_single.faa -blastDB /Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout_norecip/ -addorigin no -reciprocal no -trimal -gappyout
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/LTR_single.faa -blastDB /Applications/Databases//EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout_norecip/ -addorigin no -reciprocal no -trimal -gappyout
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/Spliceosome_single.faa -blastDB /Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout_norecip/ -addorigin no -reciprocal no -trimal -gappyout
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/TERT_single.faa -blastDB /Applications/Databases//EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa -eval 1e-10 -filter 0.25 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout_norecip/ -addorigin no -reciprocal no -trimal -gappyout


	mkdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/
	###non reciprocal, pfam, low threshold
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/multi_RDP.faa -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 10 -filter 1 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Introns.faa -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 10 -filter 1 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/LTR.faa -blastDB /Applications/Databases//EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 10 -filter 1 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/Spliceosome.faa -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 10 -filter 1 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/TERT.faa -blastDB /Applications/Databases//EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 10 -filter 1 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm


	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Introns.faa -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 1e-5 -filter 1 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/ -addorigin no -trimal -gappyout -pfam_dom RVT -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Intron_single.faa -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 1e-5 -filter 1 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/ -addorigin no -trimal -gappyout -pfam_dom RVT -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm


	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Intron_single.faa -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 1e-5 -filter 1 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/ -addorigin no -trimal -gappyout -pfam_dom RVT -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm


	outdir=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/Asgard_deepsearch_hmm/")
	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa

	pFAMDB=/Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	rm -rf $outdir
	mkdir $outdir
	cat /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout_norecip/*/*trimal.faa | grep ">" | grep -E "Asgard|Thor|Odin|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR" | grep -v "virus" | sort | uniq | sed 's/>//g' > $outdir"/candidateAsgardRDPs.txt"
	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $outdir"/candidateAsgardRDPs.txt" -outfmt %f -out $outdir"/candidateAsgardRDPs.faa"

	/Applications/BLAST_wrapper.sh -fasta $outdir"/candidateAsgardRDPs.faa" -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 1e-5 -filter 1 -type psiblast -outdir $outdir -addorigin no -trimal -gappyout -pfam_dom RVT -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta $outdir"/candidateAsgardRDPs.faa" -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 1e-5 -filter 40 -type psiblast -outdir $outdir -addorigin no -trimal -gappyout -pfam_dom RVT -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	
	cd $outdir
	/Applications/miniconda3/bin/mafft --anysymbol --reorder --auto $outdir"/candidateAsgardRDPs.faa" > $outdir"/candidateAsgardRDPs.align.faa"
	##Trimal to trim the alignment: gappyout parameter errs on the side of keeping gaps
	/Applications/miniconda3/bin/trimal -in $outdir"/candidateAsgardRDPs.align.faa" -out $outdir"/candidateAsgardRDPs.align.trimmed.faa" -gappyout
	##Build a tree using IQTREE, number of cores is automatically determined as well as appropriate substitution model
	#IQTree1
	##/Users/jmattick/Downloads/iqtree-1.6.12-MacOSX/bin/iqtree -nt AUTO -s $outdir/$genename.blast.$eval.align.trimal.faa
	#IQTree2
	/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir"/candidateAsgardRDPs.align.trimmed.faa" --alrt 1000 -B 1000


	cat /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout_norecip/*/*trimal.faa | grep ">" | grep -E "Asgard|Thor|Odin|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR" | grep -v "virus" | grep -v "Asgard_GB" | grep -v "Thor_GBS28" | sort | uniq | sed 's/>//g' > $outdir"/candidateAsgardRDPs_shortrm.txt"
	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $outdir"/candidateAsgardRDPs_shortrm.txt" -outfmt %f -out $outdir"/candidateAsgardRDPs_shortrm.faa"
	/Applications/BLAST_wrapper.sh -fasta $outdir"/candidateAsgardRDPs_shortrm.faa" -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 1e-5 -filter 1 -type psiblast -outdir $outdir -addorigin no -trimal -gappyout -pfam_dom RVT -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta $outdir"/candidateAsgardRDPs_shortrm.faa" -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 1e-5 -filter 40 -type psiblast -outdir $outdir -addorigin no -trimal -gappyout -pfam_dom RVT -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm

	/Applications/BLAST_wrapper.sh -fasta $outdir"/candidateAsgardRDPs_shortrm.faa" -blastDB /Applications/Databases/EukArchBactViral/EukArchBactViral.faa -num_hits 30 -eval 1e-10 -filter 40 -type psiblast -outdir $outdir -addorigin no -trimal -gappyout -pfam_dom RVT -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm

fi

########Take a look at all of the Asgard buddies -- uniques from the non reciprocal -- play with num hits
if [[ 1 -eq 2 ]]
then
	numhits=4
	outdir=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/Asgard_retrolook_"$numhits"/")
	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa
	rm -rf $outdir
	mkdir $outdir
	cat /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout_norecip/*/*trimal.faa | grep ">" | grep -E "Asgard|Thor|Odin|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR" | grep -v "virus" | sort | uniq | sed 's/>//g' > $outdir"/candidateAsgardRDPs.txt"
	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $outdir"/candidateAsgardRDPs.txt" -outfmt %f -out $outdir"/candidateAsgardRDPs.faa"

	while read line
	do
		samtools faidx $outdir"/candidateAsgardRDPs.faa" $line > $outdir/$line.faa
		/Applications/BLAST_wrapper.sh -fasta $outdir/$line.faa -blastDB $blastDB -eval 1e-10 -filter 0.25 -type psiblast -outdir $outdir -addorigin no -reciprocal no -trimal -gappyout -num_hits $numhits
	done < $outdir"/candidateAsgardRDPs.txt"
fi

if [[ 1 -eq 2 ]]
then
	numhits=1
	outdir=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/Asgard_retrolook_strict_"$numhits"/")
	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa
	rm -rf $outdir
	mkdir $outdir
	cat /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout_norecip/*/*trimal.faa | grep ">" | grep -E "Asgard|Thor|Odin|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR" | grep -v "virus" | sort | uniq | sed 's/>//g' > $outdir"/candidateAsgardRDPs.txt"
	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $outdir"/candidateAsgardRDPs.txt" -outfmt %f -out $outdir"/candidateAsgardRDPs.faa"

	while read line
	do
		samtools faidx $outdir"/candidateAsgardRDPs.faa" $line > $outdir/$line.faa
		/Applications/BLAST_wrapper.sh -fasta $outdir/$line.faa -blastDB $blastDB -eval 1e-10 -filter 50 -type psiblast -outdir $outdir -addorigin no -reciprocal no -trimal -gappyout -num_hits $numhits
	done < $outdir"/candidateAsgardRDPs.txt"
fi


if [[ 1 -eq 2 ]]
then
	outdir=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/Asgard_retrolook_dotplot/")
	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa
	rm -rf $outdir
	mkdir $outdir
	cat /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout_norecip/*/*trimal.faa | grep ">" | grep -E "Asgard|Thor|Odin|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR" | grep -v "virus" | sort | uniq | sed 's/>//g' > $outdir"/candidateAsgardRDPs.txt"
	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $outdir"/candidateAsgardRDPs.txt" -outfmt %f -out $outdir"/candidateAsgardRDPs.faa"
	count=1
	while read line
	do
		if [[ $count -eq 1 ]]
		then
			samtools faidx $outdir"/candidateAsgardRDPs.faa" $line > $outdir/$line.query.faa

		else
			query=$(ls $outdir/*.query.faa)
			samtools faidx $outdir"/candidateAsgardRDPs.faa" $line > $outdir/$line.subject.faa
			dotmatcher -asequence $query -bsequence $outdir/$line.subject.faa -graph png -goutfile $outdir/$line.dotplot.png

		fi
		count=$(( $count + 1 ))
	done < $outdir"/candidateAsgardRDPs.txt"
fi
###Scan whole DB........
if [[ 1 -eq 2 ]]
then
	pFAMDB=/Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa

	/Applications/miniconda3/bin/hmmscan -o /Users/jmattic1/Documents/TERT_reanalysis/All.hmm.out --tblout /Users/jmattic1/Documents/TERT_reanalysis/All.hmm.tsv --domtblout /Users/jmattic1/Documents/TERT_reanalysis/All.hmm.dom.tsv --cpu 32 $pFAMDB $blastDB
fi
###First pass HMM Scan
if [[ 1 -eq 2 ]]
then
	numhits=1
	outdir=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/Asgard_retrolook_strict_"$numhits"_hmm/")
	indir=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/Asgard_retrolook_strict_"$numhits"/")
	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa

	pFAMDB=/Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	rm -rf $outdir
	mkdir $outdir
	cat /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout_norecip/*/*trimal.faa | grep ">" | grep -E "Asgard|Thor|Odin|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR" | grep -v "virus" | sort | uniq | sed 's/>//g' > $outdir"/candidateAsgardRDPs.txt"
	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $outdir"/candidateAsgardRDPs.txt" -outfmt %f -out $outdir"/candidateAsgardRDPs.faa"

	while read line
	do
		fasta=$(\ls $indir/$line""*/$line.*.origin.faa)
		/Applications/miniconda3/bin/hmmscan -o $outdir/$line.hmm.out --tblout $outdir/$line.hmm.tsv $pFAMDB $fasta 
	done < $outdir"/candidateAsgardRDPs.txt"
fi

###RVT alignment...
if [[ 1 -eq 2 ]]
then
	lengthfilter=100

	numhits=1
	indir=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/Asgard_retrolook_strict_"$numhits"/")
	outdir=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/Asgard_retrolook_strict_PFAMalign_"$lengthfilter"/")
	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa

	pFAMDB=/Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	rm -rf $outdir
	mkdir $outdir

	cat $indir/*/*origin.faa | grep ">" | sort | uniq | sed 's/>//g' > $outdir"/candidateRDPs.txt"
	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $outdir"/candidateRDPs.txt" -outfmt %f -out $outdir"/candidateRDPs.faa"
	/Applications/miniconda3/bin/hmmscan --cpu 16 -o $outdir/Candidates.hmm.out --tblout $outdir/Candidates.hmm.tsv --domtblout $outdir/Candidates.hmm.dom.tsv $pFAMDB $outdir"/candidateRDPs.faa"
	samtools faidx $outdir"/candidateRDPs.faa"
	for i in $(cat $outdir/Candidates.hmm.tsv | grep "RVT_1" | awk '{print $3}')
	do
		nummatches=$(cat $outdir/Candidates.hmm.dom.tsv | grep $i | grep "RVT_1" | awk '{print $19-$18"\t"$0}' | awk -v lengthfilter="$lengthfilter" '$1 > lengthfilter {print $0}' | wc -l)
		if [[ $nummatches -gt 0 ]]
		then
			start=$(cat $outdir/Candidates.hmm.dom.tsv | grep $i | grep "RVT_1" | awk '{print $19-$18"\t"$0}' | sort -k1,1nr | head -1 | awk '{print $19}')
			end=$(cat $outdir/Candidates.hmm.dom.tsv | grep $i | grep "RVT_1" | awk '{print $19-$18"\t"$0}' | sort -k1,1nr | head -1 | awk '{print $20}')

			samtools faidx $outdir"/candidateRDPs.faa" $i":"$start"-"$end >> $outdir/RVT_1.domains.faa
		fi
	done

	/Applications/miniconda3/bin/mafft --anysymbol --reorder --auto $outdir/RVT_1.domains.faa > $outdir/RVT_1.domains.align.faa
	/Applications/miniconda3/bin/trimal -in $outdir/RVT_1.domains.align.faa -out $outdir/RVT_1.domains.align.trimmed.faa -gappyout
	/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir/RVT_1.domains.align.trimmed.faa --alrt 1000 -B 1000

fi

###Alpha fold -- foldseek
if [[ 1 -eq 2 ]]
then
	indir=$(echo -e "/Users/jmattic1/Documents/RDP_Folds/")
	outdir=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/Foldseek_firstiter/")
	rm -rf $outdir
	mkdir $outdir

	for i in $(ls $indir/ | grep "arch" | grep -v ".zip")
	do
		outfmt=$(echo $i".html")
		infmt=$(echo $i"*.pdb")
		/Applications/foldseek/bin/foldseek easy-search $indir/$i/$infmt /Applications/foldseek/Swiss-Prot $outdir/$outfmt tmp --format-mode 3
	done

	/Applications/foldseek/bin/foldseek easy-search /Users/jmattic1/Documents/Giardia_i_XP_0017084341_putative_Reverse_tra_f42af/*.pdb /Applications/foldseek/Swiss-Prot /Users/jmattic1/Documents/Giardia_i_XP_0017084341_putative_Reverse_tra_f42af/GiardiaAnalysis.html tmp --format-mode 3
	/Applications/foldseek/bin/foldseek easy-search /Users/jmattic1/Documents/Rhodelphis_l_DN2747_c1_g1_i2_e1881/*.pdb /Applications/foldseek/Swiss-Prot /Users/jmattic1/Documents/Rhodelphis_l_DN2747_c1_g1_i2_e1881/RhodelphisAnalysis.html tmp --format-mode 3
fi




####RDP just domain
if [[ 1 -eq 2 ]]
then

	faafile=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/GroupII_Introns_e1e-5f40nh30ni6_1523x/GroupII_Introns.blast.1e-5.origin.faa")
	hmmfile=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/GroupII_Introns_e1e-5f40nh30ni6_1523x/GroupII_Introns.blast.1e-5.hmm.tsv")

	hmmdomfile=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/GroupII_Introns_e1e-5f40nh30ni6_1523x/GroupII_Introns.blast.1e-5.dom.tsv")

	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa
	lengthfilter=75
	outdir=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/GroupII_Intron_Broad_PFAMalign_"$lengthfilter"/")

	rm -rf $outdir
	mkdir $outdir
	samtools faidx $faafile
	for i in $(cat $hmmfile | grep "RVT_1" | awk '{print $3}')
	do
		nummatches=$(cat $hmmdomfile | grep $i | grep "RVT_1" | awk '{print $19-$18"\t"$0}' | awk -v lengthfilter="$lengthfilter" '$1 > lengthfilter {print $0}' | wc -l)
		if [[ $nummatches -gt 0 ]]
		then
			start=$(cat $hmmdomfile | grep $i | grep "RVT_1" | awk '{print $19-$18"\t"$0}' | sort -k1,1nr | head -1 | awk '{print $19}')
			end=$(cat $hmmdomfile | grep $i | grep "RVT_1" | awk '{print $19-$18"\t"$0}' | sort -k1,1nr | head -1 | awk '{print $20}')

			samtools faidx $faafile $i":"$start"-"$end >> $outdir/RVT_1.domains.faa
		fi
	done

	/Applications/miniconda3/bin/mafft --anysymbol --reorder --auto $outdir/RVT_1.domains.faa > $outdir/RVT_1.domains.align.faa
	/Applications/miniconda3/bin/trimal -in $outdir/RVT_1.domains.align.faa -out $outdir/RVT_1.domains.align.trimmed.faa -gappyout
	/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir/RVT_1.domains.align.trimmed.faa --alrt 1000 -B 1000
	/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir/RVT_1.domains.align.trimmed.faa.uniqueseq.phy --alrt 1000 -B 1000
fi

####RDP combine other asgards...
if [[ 1 -eq 2 ]]
then

	oldfaafile=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/gappyout_pfam/GroupII_Introns_e1e-5f40nh30ni6_1523x/GroupII_Introns.blast.1e-5.origin.faa")

	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa
	lengthfilter=75
	outdir=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/GroupII_Intron_Broad_PFAMalign_Asgardadded_"$lengthfilter"/")

	rm -rf $outdir
	mkdir $outdir

	removalgenes=$(cat $oldfaafile | grep ">" | awk '{print $1}' | grep -E "Asgard|Thor|Odin|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR" | grep -v "virus" | sort | uniq | sed 's/>//g' | tr '\n' '|' | sed 's/\|$//g')
	cat /Users/jmattic1/Documents/TERT_reanalysis/singlegene_gappyout_norecip/*/*trimal.faa | grep ">" | awk '{print $1}' | grep -E "Asgard|Thor|Odin|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR" | grep -vE $removalgenes | grep -v "virus" | grep -v "Asgard_GB" | grep -v "Thor_GBS28" | sort | uniq | sed 's/>//g' > $outdir"/candidateAsgardRDPs.txt"
	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $outdir"/candidateAsgardRDPs.txt" -outfmt %f -out $outdir"/candidateAsgardRDPs.faa"

	cat $outdir"/candidateAsgardRDPs.faa" $oldfaafile > $outdir"/candidateCombinedRDPs.faa"

	faafile=$(echo $outdir"/candidateCombinedRDPs.faa")
	samtools faidx $faafile
	/Applications/miniconda3/bin/hmmscan -o $outdir"/candidateCombinedRDPs.hmm.out" --tblout $outdir"/candidateCombinedRDPs.hmm.tsv" --domtblout $outdir"/candidateCombinedRDPs.hmm.dom.tsv" --cpu 32 /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm $faafile

	hmmfile=$(echo -e $outdir"/candidateCombinedRDPs.hmm.tsv")

	hmmdomfile=$(echo -e $outdir"/candidateCombinedRDPs.hmm.dom.tsv")

	for i in $(cat $hmmfile | grep "RVT_1" | awk '{print $3}')
	do
		nummatches=$(cat $hmmdomfile | grep $i | grep "RVT_1" | awk '{print $19-$18"\t"$0}' | awk -v lengthfilter="$lengthfilter" '$1 > lengthfilter {print $0}' | wc -l)
		if [[ $nummatches -gt 0 ]]
		then
			start=$(cat $hmmdomfile | grep $i | grep "RVT_1" | awk '{print $19-$18"\t"$0}' | sort -k1,1nr | head -1 | awk '{print $19}')
			end=$(cat $hmmdomfile | grep $i | grep "RVT_1" | awk '{print $19-$18"\t"$0}' | sort -k1,1nr | head -1 | awk '{print $20}')

			samtools faidx $faafile $i":"$start"-"$end >> $outdir/RVT_1.domains.faa
		fi
	done

	/Applications/miniconda3/bin/mafft --anysymbol --reorder --auto $outdir/RVT_1.domains.faa > $outdir/RVT_1.domains.align.faa
	/Applications/miniconda3/bin/trimal -in $outdir/RVT_1.domains.align.faa -out $outdir/RVT_1.domains.align.trimmed.faa -gappyout
	/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir/RVT_1.domains.align.trimmed.faa --alrt 1000 -B 1000
	/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir/RVT_1.domains.align.trimmed.faa.uniqueseq.phy --alrt 1000 -B 1000
fi


####RDP ALL DOMAINS
if [[ 1 -eq 2 ]]
then
	lengthfilter=75

	outdir=$(echo -e "/Users/jmattic1/Documents/TERT_reanalysis/All_RVT_classification_"$lengthfilter"/")
	rm -rf $outdir
	mkdir $outdir
	cd $outdir
	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa

	cat /Users/jmattic1/Documents/TERT_reanalysis/All.hmm.dom.tsv | awk '$2 == "PF00078.30" {print $0}' | awk '$8 > 100 {print $0}' | awk '{print $19-$18"\t"$0}' | awk -v lengthfilter="$lengthfilter" '$1 > lengthfilter {print $5}' | sort | uniq > $outdir/All.rvts.list

	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch $outdir/All.rvts.list -outfmt %f -out $outdir/All.rvts.faa

	/Applications/miniconda3/bin/mafft --anysymbol --reorder --auto $outdir/All.rvts.faa > $outdir/All.rvts.align.faa
	/Applications/miniconda3/bin/trimal -in $outdir/All.rvts.align.faa -out $outdir/All.rvts.align.trimmed.faa -gappyout
	/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir/All.rvts.align.trimmed.faa --alrt 1000 -B 1000
	/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir/All.rvts.align.trimmed.faa.uniqueseq.phy --alrt 1000 -B 1000
fi


if [[ 1 -eq 2 ]]
then
	numhits=1
	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa

	for i in $(\ls /Users/jmattic1/Documents/NuclearProteins/*/*.faa)
	do
		protein=$(echo $i | awk -F "/" '{print $NF}' | sed 's/.faa//g')
		currentdir=$(echo $i | awk -F "/" 'NF{NF-=1};1'  | sed 's/ /\//g' | awk '{print $0"/"}' )
		outdir=$(echo -e $currentdir"/"$protein"_analysis/")
		rm -rf $outdir
		mkdir $outdir
		cd $outdir
 
		/Applications/BLAST_wrapper.sh -fasta $i -blastDB $blastDB -eval 1e-10 -filter 30 -type psiblast -outdir $outdir -addorigin no -reciprocal no -trimal -gappyout -num_hits $numhits -num_iters 10

	done
	i=/Users/jmattic1/Documents/NuclearProteins/CoatNucleoporinComplex/Nup37.faa
	protein=$(echo $i | awk -F "/" '{print $NF}' | sed 's/.faa//g')
	currentdir=$(echo $i | awk -F "/" 'NF{NF-=1};1'  | sed 's/ /\//g' | awk '{print $0"/"}' )
	outdir=$(echo -e $currentdir"/"$protein"_analysis_reciprocal/")
	rm -rf $outdir
	mkdir $outdir
	cd $outdir
	###Nup37 reciprocal check...
	/Applications/BLAST_wrapper.sh -fasta $i -blastDB $blastDB -eval 1e-10 -filter 30 -type psiblast -outdir $outdir -addorigin no -reciprocal /Applications/Databases/ReciprocalDB/ReciprocalStandard.protein.db.faa -trimal -gappyout -num_hits $numhits -num_iters 10

fi

if [[ 1 -eq 2 ]]
then
	numhits=1
	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa

	for i in $(\ls /Users/jmattic1/Documents/NuclearProteins_notrees/*/*.faa)
	do
		protein=$(echo $i | awk -F "/" '{print $NF}' | sed 's/.faa//g')
		currentdir=$(echo $i | awk -F "/" 'NF{NF-=1};1'  | sed 's/ /\//g' | awk '{print $0"/"}' )
		outdir=$(echo -e $currentdir"/"$protein"_analysis/")
		rm -rf $outdir
		mkdir $outdir
		cd $outdir
 
		/Applications/BLAST_wrapper_notree.sh -fasta $i -blastDB $blastDB -eval 1e-10 -filter 30 -type psiblast -outdir $outdir -addorigin no -reciprocal no -trimal -gappyout -num_hits $numhits -num_iters 10

	done

fi

if [[ 1 -eq 2 ]]
then
	####Make trees in each directory where we find hits to Asgard
	for i in $(cat /Users/jmattic1/Documents/NuclearProteins_notrees/*/*/*/*.filtered.blast | grep -E "Asgard|Thor|Odin|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR" | awk '{print $1}' | sort | uniq)
	do
		fasta=$(grep -l $i /Users/jmattic1/Documents/NuclearProteins_notrees/*/*.faa)
		genename=$(echo $fasta | awk -F '/' '{print $NF}' | sed 's/.faa//g')
		treefile=$(\ls /Users/jmattic1/Documents/NuclearProteins_notrees/*/*/*/ | grep $genename | grep "contree")

		if [[ $treefile == "" ]]
		then
			/Applications/miniconda3/bin/iqtree -T AUTO -s /Users/jmattic1/Documents/NuclearProteins_notrees/*/$genename*/$genename*/$genename.*.align.trimal.faa --alrt 1000 -B 1000
		fi
	done


fi


###Uniprot database download
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/Users/jmattic1/Documents/Dessication_Asgards_notrees/Plasma_Membrane/")
	intable=$(echo "/Users/jmattic1/Downloads/PlasmaMembrane_GO_Archaea.tsv")
	rm -rf $outdir
	mkdir $outdir


	for genename in $(cat $intable | awk '{print $3}' | sort | uniq)
	do
		ProtID=$(cat $intable | awk -v genename="$genename" '$3 == genename {print $2}' | head -1)
		curl -H "Accept: text/plain; format=fasta" "https://rest.uniprot.org/uniprotkb/$ProtID" > $outdir"/"$ProtID".faa"
	done

	numhits=1
	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa

	for i in $(\ls $outdir/*.faa)
	do
		protein=$(echo $i | awk -F "/" '{print $NF}' | sed 's/.faa//g')
		currentdir=$(echo $i | awk -F "/" 'NF{NF-=1};1'  | sed 's/ /\//g' | awk '{print $0"/"}' )
		newoutdir=$(echo -e $currentdir"/"$protein"_analysis/")
		rm -rf $newoutdir
		mkdir $newoutdir
		cd $newoutdir
 
		/Applications/BLAST_wrapper_notree.sh -fasta $i -blastDB $blastDB -eval 1e-10 -filter 30 -type psiblast -outdir $newoutdir -addorigin no -reciprocal no -trimal -gappyout -num_hits $numhits -num_iters 10 &

	done


	for i in $(cat <(cat $outdir/*/*/*.filtered.blast | grep -E "Homo|Arab" | awk '{print $1}' | sort | uniq) <(cat $outdir/*/*/*.filtered.blast | grep -E "Colpo|Guill" | awk '{print $1}' | sort | uniq) <(cat $outdir/*/*/*.filtered.blast | grep -E "Asgard|Thor|Odinarch|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR" | awk '{print $1}' | sort | uniq) | sort | uniq -c | awk '$1 > 2 {print $2}')
	do
		fasta=$(grep -l $i $outdir/*.faa)
		protein=$(echo $fasta | awk -F '/' '{print $NF}' | sed 's/.faa//g')
		currentdir=$(echo $fasta | awk -F "/" 'NF{NF-=1};1'  | sed 's/ /\//g' | awk '{print $0"/"}' )
		newoutdir=$(echo -e $currentdir"/"$protein"_analysis_reciprocal/")
		rm -rf $newoutdir
		mkdir $newoutdir
		cd $newoutdir

		/Applications/BLAST_wrapper.sh -fasta $fasta -blastDB $blastDB -eval 1e-10 -filter 30 -type psiblast -outdir $newoutdir -addorigin no -reciprocal /Applications/Databases/ReciprocalDB/ReciprocalStandard.protein.db.faa -trimal -gappyout -num_hits $numhits -num_iters 10 &


	done

###Non reciprocal tree making
	for i in $(cat <(cat $outdir/*/*/*.filtered.blast | grep -E "Homo|Arab" | awk '{print $1}' | sort | uniq) <(cat $outdir/*/*/*.filtered.blast | grep -E "Colpo|Guill" | awk '{print $1}' | sort | uniq) <(cat $outdir/*/*/*.filtered.blast | grep -E "Asgard|Thor|Odinarch|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR" | awk '{print $1}' | sort | uniq) | sort | uniq -c | awk '$1 > 2 {print $2}')
	do
		fasta=$(grep -l $i $outdir/*.faa)
		protein=$(echo $fasta | awk -F '/' '{print $NF}' | sed 's/.faa//g')
		currentdir=$(echo $fasta | awk -F "/" 'NF{NF-=1};1'  | sed 's/ /\//g' | awk '{print $0"/"}' )
		newoutdir=$(echo -e $currentdir"/"$protein"_analysis_tree_nonreciprocal/")
		rm -rf $newoutdir
		mkdir $newoutdir
		cd $newoutdir

		/Applications/BLAST_wrapper.sh -fasta $fasta -blastDB $blastDB -eval 1e-10 -filter 30 -type psiblast -outdir $newoutdir -addorigin no -reciprocal no -trimal -gappyout -num_hits $numhits -num_iters 10 &


	done


	###Make these trees with full arch tree...
	outdir=$(echo "/Users/jmattic1/Documents/Dessication_Asgards_notrees/Plasma_Membrane_FullArch/")
	rm -rf $outdir
	mkdir $outdir
	cp /Users/jmattic1/Documents/Dessication_Asgards_notrees/Plasma_Membrane/*.faa $outdir
	numhits=1
	blastDB=/Applications/Databases/ArchBactBLAST_DB_v2/CombinedBactArch.faa
	
	for i in $(\ls $outdir/*.faa)
	do
		protein=$(echo $i | awk -F "/" '{print $NF}' | sed 's/.faa//g')
		currentdir=$(echo $i | awk -F "/" 'NF{NF-=1};1'  | sed 's/ /\//g' | awk '{print $0"/"}' )
		newoutdir=$(echo -e $currentdir"/"$protein"_analysis/")
		rm -rf $newoutdir
		mkdir $newoutdir
		cd $newoutdir
 
		/Applications/BLAST_wrapper.sh -fasta $i -blastDB $blastDB -eval 1e-10 -filter 30 -type psiblast -outdir $newoutdir -addorigin no -reciprocal no -trimal -gappyout -num_hits $numhits -num_iters 10 &

	done
	
fi

####Evaluate the hits
if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/Users/jmattic1/Documents/Dessication_Asgards_notrees/Plasma_Membrane/")

	rm -rf /Users/jmattic1/Documents/DeepBranchingMembraneProteins
	mkdir /Users/jmattic1/Documents/DeepBranchingMembraneProteins
	for i in $(cat <(grep -lE "Homo|Arab" $outdir/*_tree_nonreciprocal/*/*.origin.faa) <(grep -lE "Colpo|Guill" $outdir/*_tree_nonreciprocal/*/*.origin.faa) <(grep -lE "Asgard|Thor|Odinarch|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR" $outdir/*_tree_nonreciprocal/*/*.origin.faa) | sort | uniq -c | awk '$1 > 2 {print $2}')
	do
		dir=$(echo $i | awk -F "/" 'NF{NF-=2};1'  | sed 's/ /\//g' | awk '{print $0"/"}')
		cp -r $dir /Users/jmattic1/Documents/DeepBranchingMembraneProteins/
	done

	fastadir="/Users/jmattic1/Documents/Dessication_Asgards_notrees/Plasma_Membrane/"
	for i in $(\ls /Users/jmattic1/Documents/DeepBranchingMembraneProteins/ | grep "x")
	do
		fastaname=$(echo $i | awk -F "_" '{print $1".faa"}')
		genelength=$(perl /Applications/residues.pl $fastadir"/"$fastaname | awk '{print $2}')
		echo $i
		perl /Applications/residues.pl /Users/jmattic1/Documents/DeepBranchingMembraneProteins/$i/*align.trimal.faa

	done
fi

if [[ 1 -eq 2 ]]
then
	####Make trees in each directory where we find hits to Asgard
	numhits=1
	blastDB=/Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa

	for i in $(cat /Users/jmattic1/Documents/NuclearProteins_notrees/*/*/*/*.filtered.blast | grep -E "Asgard|Thor|Odin|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|WOR|Hod|Sif" | awk '{print $1}' | sort | uniq)
	do
		fasta=$(grep -l $i /Users/jmattic1/Documents/NuclearProteins_notrees/*/*.faa)
		protein=$(echo $fasta | awk -F '/' '{print $NF}' | sed 's/.faa//g')
		currentdir=$(echo $fasta | awk -F "/" 'NF{NF-=1};1'  | sed 's/ /\//g' | awk '{print $0"/"}' )
		outdir=$(echo -e $currentdir"/"$protein"_analysis_reciprocal/")
		rm -rf $outdir
		mkdir $outdir
		cd $outdir

		/Applications/BLAST_wrapper_notree.sh -fasta $fasta -blastDB $blastDB -eval 1e-10 -filter 30 -type psiblast -outdir $outdir -addorigin no -reciprocal /Applications/Databases/ReciprocalDB/ReciprocalStandard.protein.db.faa -trimal -gappyout -num_hits $numhits -num_iters 10 &

	done

fi



#####Re-do tree analysis with improved BLAST wrapper and new database
if [[ 1 -eq 2 ]]
then
	conda activate /Applications/miniconda3
	rm -rf /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/
	mkdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/multi_RDP.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 30 -eval 1e-5 -filter 40 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm 
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Introns.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 30 -eval 1e-10 -filter 60 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/LTR.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 30 -eval 1e-5 -filter 40 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/Spliceosome.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 30 -eval 1e-5 -filter 40 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/multi_RDP.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 30 -eval 1e-5 -filter 40 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -SR4 yes
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Introns.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 30 -eval 1e-10 -filter 60 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -SR4 yes
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/LTR.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 30 -eval 1e-5 -filter 40 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -SR4 yes
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/Spliceosome.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 30 -eval 1e-5 -filter 40 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm


	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Introns.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Introns.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -SR4 yes

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Introns.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -SR4 yes

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Introns.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -SR4 yes

fi

#####Re-do tree analysis with improved BLAST wrapper and new database
if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/
	mkdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/

	blastDB="/Applications/Databases/Euk_ArchBact_Condensed_v3.1/Euk_ArchBact_Condensed_v3.1.faa"
	cat /Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/*/*.align.trimal.faa | grep ">" | grep -E "Asgard|Thor|Odin|Hel_|Njord|Jord|Heim|Loki|Kor|Prometh|Baldr|Sif|Hod" | sed 's/>//g' > /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardList.txt
	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardList.txt -outfmt %f -out /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 50 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 50 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -SR4 yes

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1/Euk_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 50 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1/Euk_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 50 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -SR4 yes

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1/Euk_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-10 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1/Euk_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-10 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -SR4 yes

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-10 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-10 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -SR4 yes

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 5 -eval 1e-10 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 5 -eval 1e-10 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -SR4 yes


	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 5 -eval 1e-10 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	
	longestseq=$(perl /Applications/residues.pl /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa | sort -k2,2nr | awk '{print $1}' | head -1)
	heimseq=$(echo "CHeimdall_archaeon-L1_MCY3414279.1")
	samtools faidx /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa $longestseq | sed 's/@/_/g' | sed 's/\|/_/g' | sed 's/\//_/g' > /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.single.faa
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.single.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 5 -eval 1e-10 -filter 50 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	samtools faidx /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa $heimseq | sed 's/@/_/g' | sed 's/\|/_/g' | sed 's/\//_/g' > /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 5 -eval 1e-10 -filter 50 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-10 -filter 50 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1/Euk_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-10 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1/Euk_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-10 -filter 50 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm

fi



###Make a full Asgard DB

if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/TERT_reanalysis/new_db_fullasgard_asgardqueries/
	mkdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_fullasgard_asgardqueries/

	rm -rf /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/
	mkdir /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/

	cat /Applications/Databases/Euk_ArchBact_Condensed_v3.1/Euk_ArchBact_Condensed_v3.1.faa > /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/Euk_ArchBact_Condensed_v3.1_fullasgard.faa
	cat /Users/jmattic1/Documents/Laura_proteomes/*.faa | sed 's/>/>lcl\|/g' | sed 's/Archaeon_WORA1..WORA1/Sif_archaea/g' | sed 's/Archaeon_WORB2..WORB2/Hermod_archaea/g' >> /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/Euk_ArchBact_Condensed_v3.1_fullasgard.faa
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/Euk_ArchBact_Condensed_v3.1_fullasgard.faa -dbtype 'prot' -parse_seqids

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/Euk_ArchBact_Condensed_v3.1_fullasgard.faa -num_hits 1 -eval 1e-10 -filter 50 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm


	/Applications/aa_recoder.py -i /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I_e1e-10f30nh1ni6_408x290_5883/AsgardG2I.blast.1e-10.align.trimal.faa -re SR4 -o /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I_e1e-10f30nh1ni6_408x290_5883/AsgardG2I.blast.1e-10.align.trimal.SR4.faa

	/Applications/miniconda3/bin/iqtree -T AUTO -s /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I_e1e-10f30nh1ni6_408x290_5883/AsgardG2I.blast.1e-10.align.trimal.SR4.faa.fas --alrt 1000 -B 1000 -wsr

fi

####G2I characterization -- protein

if [[ 1 -eq 2 ]]
then
	####new_db_asgardqueries/AsgardG2I_e1e-10f30nh1ni6_408x290_5883
	rm -rf /Users/jmattic1/Documents/TERT_reanalysis/PFAM_G2I_Classification/
	mkdir /Users/jmattic1/Documents/TERT_reanalysis/PFAM_G2I_Classification/
	pFAMDB="/Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm"
	/Applications/miniconda3/bin/hmmscan -o /Users/jmattic1/Documents/TERT_reanalysis/PFAM_G2I_Classification/AsgardG2I_e1e-10f30nh1ni6_408x290_5883.hmm.out --tblout /Users/jmattic1/Documents/TERT_reanalysis/PFAM_G2I_Classification/AsgardG2I_e1e-10f30nh1ni6_408x290_5883.hmm.tsv --domtblout /Users/jmattic1/Documents/TERT_reanalysis/PFAM_G2I_Classification/AsgardG2I_e1e-10f30nh1ni6_408x290_5883.dom.tsv --cpu 32 $pFAMDB /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I_e1e-10f30nh1ni6_408x290_5883//AsgardG2I.blast.1e-10.origin.nostar.faa
fi

####G2I characterization -- RNA secondary structure
###Check if we have RNA for all species

if [[ 1 -eq 2 ]]
then
	fasta=$(\ls /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I_e1e-10f30nh1ni6_408x290_5883//AsgardG2I.blast.1e-10.origin.nostar.faa)

	for i in $(cat $fasta | grep ">" | sed 's/>//g')
	do
		species=$(echo $i | awk -F "_" '{print $1"_"$2}')

		if [[ $(cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.euk.filtered.csv | grep $species | awk -F "," '{print $3}' | wc -l) -gt 0 ]]
		then
			#echo $species
			#cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.euk.filtered.csv | grep $species | awk -F "," '{print $3}'
			#\ls /Users/jmattic1/Documents/Transcriptomes/* | grep $species
			transcriptome1=$(\ls /Users/jmattic1/Documents/Transcriptomes/* | grep $species | wc -l)
			#\ls /Users/jmattic1/Documents/Transcriptomes/*/* | grep $species
			transcriptome2=$(\ls /Users/jmattic1/Documents/Transcriptomes/*/* | grep $species | wc -l)
			if [[ $(cat /Users/jmattic1/Documents/EukProt_included_data_sets.v03.2021_11_22.txt | grep $species | awk '{print $1}' | wc -l) -gt 0 ]]
			then
				eukprot=$(cat /Users/jmattic1/Documents/EukProt_included_data_sets.v03.2021_11_22.txt | grep $species | awk '{print $1}' | wc -l)
			elif [[ $(cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.euk.filtered.csv | grep $species | awk -F "," '{print $3}' | grep -E "GCA|GCF" | wc -l) -gt 0 ]]
			then
				genbank=$(cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.euk.filtered.csv | grep $species | awk -F "," '{print $3}' | grep -E "GCA|GCF" | wc -l)
				#echo "Present in Genbank"
			fi
			total=$(( $transcriptome1 + $transcriptome2 + $eukprot + $genbank))
			if [[ $total -eq 0 ]]
			then
				echo $species" has no home!!!"
			fi
			#echo -e "\n----------------------------------------\n"
		fi
	done

fi

if [[ 1 -eq 2 ]]
then
###Archaeal genbank
	outdir=$(echo "/Users/jmattic1/Documents/ArchaealGenbankCSV_all/")
	rm -rf $outdir
	mkdir $outdir
	fasta=$(\ls /Applications/Databases/ArchBactBLAST_DB_v2_xAsgard/CombinedBactArchxAsgard.faa)
	for i in $(cat $fasta | grep ">" | sed 's/lcl\|//g' | awk -F "_" '{print $1"_"$2}' | sort | uniq | sed 's/>//g')
	do
		species=$(echo $i)
		if [[ $(cat /Users/jmattic1/Documents/FinishedFastas.txt | awk -v species="$species" '$2 == species".faa" {print $1}' | wc -l) -gt 0 ]]
		then
			translatedspecies=$(cat /Users/jmattic1/Documents/FinishedFastas.txt | awk -v species="$species" '$2 == species".faa" {print $1}' | sed 's/.faa//g')
			if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/ | grep "latest_assembly_versions" | wc -l) -gt 0 ]]
			then
				if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/latest_assembly_versions/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
				then
					proteinfile=$(ls -lt /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/latest_assembly_versions/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
				fi
			elif [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/ | grep "all_assembly_versions" | wc -l) -gt 0 ]]
			then
				if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
				then
					proteinfile=$(ls -lt /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
				elif [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/ | grep "suppressed" | wc -l) -gt 0 ]]
				then
					if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/suppressed/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
					then
						proteinfile=$(ls -lt /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/suppressed/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
					fi
				fi
			fi
			accession=$(echo $proteinfile | awk -F "/" '{print $(NF-1)}')
			cat /Users/jmattic1/Documents/FinishedFastas.txt | awk -v accession="$accession" -v species="$species" '$2 == species".faa" {print $1","$2","species","accession}' >> $outdir/Arch.Genbank.translation.csv
		fi
	done
fi

if [[ 1 -eq 2 ]]
then
	#####Adding back in Genbank viruses
	#Pagoda_mosaic-L1_YP_009041481.1
	#mkdir /Users/jmattic1/Documents/ViralRefSeq/viral/Pagoda_yellow_mosaic_virus
	#mkdir /Users/jmattic1/Documents/ViralRefSeq/viral/Pagoda_yellow_mosaic_virus/latest_assembly_versions/
	#mkdir /Users/jmattic1/Documents/ViralRefSeq/viral/Pagoda_yellow_mosaic_virus/latest_assembly_versions/GCF_000923515.1/
	#cp -r /Applications/Databases/EukProtTranscriptomes/ncbi_dataset_\(2\)/ncbi_dataset/data/GCF_000923515.1/ /Users/jmattic1/Documents/ViralRefSeq/viral/Pagoda_yellow_mosaic_virus/latest_assembly_versions/GCF_000923515.1/
	#gzip /Users/jmattic1/Documents/ViralRefSeq/viral/Pagoda_yellow_mosaic_virus/latest_assembly_versions/GCF_000923515.1/*.faa
	#gzip /Users/jmattic1/Documents/ViralRefSeq/viral/Pagoda_yellow_mosaic_virus/latest_assembly_versions/GCF_000923515.1/*.fna
	#echo -e "Pagoda_yellow_mosaic_virus.faa\tPagoda_mosaic-L1.faa" >> /Users/jmattic1/Documents/FinishedFastasVirus_mod.txt

	#Wisteria_-L1_YP_009352866.1
	#sed -i '' 's/Wisteria_mosaic-L1/Wisteria_-L1/g' /Users/jmattic1/Documents/FinishedFastasVirus_mod.txt


	##Adding mitochondria to physcomitrella
	#cat /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/ncbi_dataset_\(36\)/ncbi_dataset/data/GCF_000002425.4/rna.fna /Applications/Databases/EukProtTranscriptomes/sequence.txt > /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/ncbi_dataset_\(36\)/ncbi_dataset/data/GCF_000002425.4/rna.mito.fna
	#mv /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/ncbi_dataset_\(36\)/ncbi_dataset/data/GCF_000002425.4/rna.fna /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/ncbi_dataset_\(36\)/ncbi_dataset/data/GCF_000002425.4/rna.fna.old
	
	#cat /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/ncbi_dataset_\(83\)/ncbi_dataset/data/GCF_000146045.2/rna.fna /Applications/Databases/EukProtTranscriptomes/SCerevisiaeMito.fasta > /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/ncbi_dataset_\(83\)/ncbi_dataset/data/GCF_000146045.2/rna.mito.fna
	#mv /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/ncbi_dataset_\(83\)/ncbi_dataset/data/GCF_000146045.2/rna.fna /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/ncbi_dataset_\(83\)/ncbi_dataset/data/GCF_000146045.2/rna.fna.old

	#cat /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/ncbi_dataset_\(85\)/ncbi_dataset/data/GCF_000002945.1/rna.fna /Applications/Databases/EukProtTranscriptomes/Schizosacch_mito.fasta > /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/ncbi_dataset_\(85\)/ncbi_dataset/data/GCF_000002945.1/rna.mito.fna
	#mv /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/ncbi_dataset_\(85\)/ncbi_dataset/data/GCF_000002945.1/rna.fna /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/ncbi_dataset_\(85\)/ncbi_dataset/data/GCF_000002945.1/rna.fna.old

	#conda activate /Applications/miniconda3/envs/samtools
	#originfasta=$(\ls /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I_e1e-10f30nh1ni6_408x290_5883//AsgardG2I.blast.1e-10.origin.nostar.faa)

	originfasta=$(\ls /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/GroupII_Intron_single_e1e-5f60nh1ni6_35x437_14745/GroupII_Intron_single.blast.1e-5.origin.nostar.faa)

	wd=$(echo $originfasta | awk 'BEGIN{FS=OFS="/"}{NF--; print}')
	outdir=$(echo $wd"/RFAM_Classification/")
	rm -rf $outdir
	mkdir $outdir
	touch $outdir/RNAsequences.fna
	originfasta=$(\ls /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I_e1e-10f30nh1ni6_408x290_5883//AsgardG2I.blast.1e-10.origin.nostar.faa)
	cat $originfasta | sed 's/\//_/g' > $outdir/Blastwrapper.proteins.faa
	fasta=$(echo $outdir/Blastwrapper.proteins.faa)
	samtools faidx $fasta
	for i in $(cat $fasta | grep ">" | sed 's/>//g')
	do
		echo $i
		species=$(echo $i | awk -F "_" '{print $1"_"$2}')
		gene=$(echo $i | cut -f 3- -d '_')
		###Regular Arch
		if [[ $(cat /Users/jmattic1/Documents/FinishedFastas.txt | awk -v species="$species" '$2 == species".faa" {print $1}' | wc -l) -gt 0 ]]
		then
			translatedspecies=$(cat /Users/jmattic1/Documents/FinishedFastas.txt | awk -v species="$species" '$2 == species".faa" {print $1}' | sed 's/.faa//g')
			if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/ | grep "latest_assembly_versions" | wc -l) -gt 0 ]]
			then
				if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/latest_assembly_versions/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
				then
					proteinfile=$(ls -lt /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/latest_assembly_versions/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
				fi
			elif [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/ | grep "all_assembly_versions" | wc -l) -gt 0 ]]
			then
				if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
				then
					proteinfile=$(ls -lt /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
				elif [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/ | grep "suppressed" | wc -l) -gt 0 ]]
				then
					if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/suppressed/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
					then
						proteinfile=$(ls -lt /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/suppressed/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
					fi
				fi
			fi
			accession=$(echo $proteinfile | awk -F "/" '{print $(NF-1)}')
			path=$(echo $proteinfile | awk 'BEGIN{FS=OFS="/"}{NF--; print}')
			###Transcript vs "genomic"
			##RNAfile=$(\ls $path/*cds_from_genomic.fna.gz)
			RNAfile=$(\ls $path/*_genomic.fna.gz | grep -v "from")
			cp $RNAfile $outdir/$species.$gene.rna.fna.gz
			gzip -d $outdir/$species.$gene.rna.fna.gz
			if [[ $(cat $outdir/$species.$gene.rna.fna | grep $gene | awk '{print $1}' | sed 's/>//g' | wc -l) -gt 0 ]]
			then
				rnaname=$(cat $outdir/$species.$gene.rna.fna | grep $gene | awk '{print $1}' | sed 's/>//g')
			else
				/Applications/miniconda3/bin/makeblastdb -in $outdir/$species.$gene.rna.fna -dbtype 'nucl' -parse_seqids
				samtools faidx $fasta $i > $outdir/protein.temp.faa
				/Applications/miniconda3/bin/tblastn -query $outdir/protein.temp.faa -db $outdir/$species.$gene.rna.fna -outfmt 6 -max_hsps 1 -num_threads 8 -out $outdir/$species.$gene.rna.blast -max_target_seqs 5129560 -evalue 1e-5
				if [[ $(cat $outdir/$species.$gene.rna.blast | head -1 | awk '{print $2}' | wc -l) -gt 0 ]]
				then
					rnaname=$(cat $outdir/$species.$gene.rna.blast | head -1 | awk '{print $2}')
				else
					echo "No RNA matches the protein... CHECK ME: "$i
					#exit
				fi
				rm -f $outdir/protein.temp.faa
			fi
			echo ">"$i >> $outdir/RNAsequences.fna
			samtools faidx $outdir/$species.$gene.rna.fna $rnaname | grep -v ">" >> $outdir/RNAsequences.fna
		###All Viral
		elif [[ $(cat /Users/jmattic1/Documents/FinishedFastasVirus_mod.txt | awk -v species="$species" '$2 == species".faa" {print $1}' | wc -l) -gt 0 ]]
		then
			translatedspecies=$(cat /Users/jmattic1/Documents/FinishedFastasVirus_mod.txt | awk -v species="$species" '$2 == species".faa" {print $1}' | sed 's/.faa//g')
			if [[ $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral/$translatedspecies/ | grep "latest_assembly_versions" | wc -l) -gt 0 ]]
			then
				if [[ $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral/$translatedspecies/latest_assembly_versions/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
				then
					proteinfile=$(ls -lt /Users/jmattic1/Documents/ViralRefSeq/viral/$translatedspecies/latest_assembly_versions/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
				fi
			elif [[ $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral/$translatedspecies/ | grep "all_assembly_versions" | wc -l) -gt 0 ]]
			then
				if [[ $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral/$translatedspecies/all_assembly_versions/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
				then
					proteinfile=$(ls -lt /Users/jmattic1/Documents/ViralRefSeq/viral/$translatedspecies/all_assembly_versions/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
				elif [[ $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral/$translatedspecies/all_assembly_versions/ | grep "suppressed" | wc -l) -gt 0 ]]
				then
					if [[ $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral/$translatedspecies/all_assembly_versions/suppressed/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
					then
						proteinfile=$(ls -lt /Users/jmattic1/Documents/ViralRefSeq/viral/$translatedspecies/all_assembly_versions/suppressed/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
					fi
				fi
			fi
			accession=$(echo $proteinfile | awk -F "/" '{print $(NF-1)}')
			path=$(echo $proteinfile | awk 'BEGIN{FS=OFS="/"}{NF--; print}')
			RNAfile=$(\ls $path/*cds_from_genomic.fna.gz)
			cp $RNAfile $outdir/$species.$gene.rna.fna.gz
			gzip -d $outdir/$species.$gene.rna.fna.gz
			rnaname=$(cat $outdir/$species.$gene.rna.fna | grep $gene | awk '{print $1}' | sed 's/>//g')
			echo ">"$i >> $outdir/RNAsequences.fna
			samtools faidx $outdir/$species.$gene.rna.fna $rnaname | grep -v ">" >> $outdir/RNAsequences.fna
		###Eukaryotic evaluation
		elif [[ $(cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.eukbact.filtered.csv | awk 'BEGIN{FS=OFS=","} {gsub(/\./, "", $1)} 1' | sed 's/Porospora_cf_gigantea_B/Porospora_gigantea/g' | grep $species | awk -F "," '{print $3}' | wc -l) -gt 0 ]]
		then
			if [[ $(cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.eukbact.filtered.csv | awk 'BEGIN{FS=OFS=","} {gsub(/\./, "", $1)} 1' | sed 's/Porospora_cf_gigantea_B/Porospora_gigantea/g' | grep $species | awk -F "," '{print $3}' | grep -E "GCA|GCF" | wc -l) -gt 0 ]]
			then
				genbank=$(cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.eukbact.filtered.csv | awk 'BEGIN{FS=OFS=","} {gsub(/\./, "", $1)} 1' | sed 's/Porospora_cf_gigantea_B/Porospora_gigantea/g' | grep $species | awk -F "," '{print $3}' | awk -F "_" '{print $1"_"$2}')
				eukfasta=$(\ls /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/*/*/*/$genbank/*.fna)
				accession=$(echo $eukfasta | awk -F "/" '{print $(NF-1)}')
				cat $eukfasta | sed 's/lcl\|//g' | awk '{if($0 ~ /^>/)  print $1"_"NR; else print $0;}' > $outdir/$species.$gene.rna.fna
				if [[ $(cat $outdir/$species.$gene.rna.fna | grep $gene | awk '{print $1}' | sed 's/>//g' | wc -l) -gt 0 ]]
				then
					rnaname=$(cat $outdir/$species.$gene.rna.fna | grep $gene | awk '{print $1}' | sed 's/>//g')
				else
					/Applications/miniconda3/bin/makeblastdb -in $outdir/$species.$gene.rna.fna -dbtype 'nucl' -parse_seqids
					samtools faidx $fasta $i > $outdir/protein.temp.faa
					/Applications/miniconda3/bin/tblastn -query $outdir/protein.temp.faa -db $outdir/$species.$gene.rna.fna -outfmt 6 -max_hsps 1 -num_threads 8 -out $outdir/$species.$gene.rna.blast -max_target_seqs 5129560 -evalue 1e-5
					if [[ $(cat $outdir/$species.$gene.rna.blast | head -1 | awk '{print $2}' | wc -l) -gt 0 ]]
					then
						rnaname=$(cat $outdir/$species.$gene.rna.blast | head -1 | awk '{print $2}')
					else
						echo "No RNA matches the protein... CHECK ME: "$i
						exit
					fi
					rm -f $outdir/protein.temp.faa
				fi
				echo ">"$i >> $outdir/RNAsequences.fna
				samtools faidx $outdir/$species.$gene.rna.fna $rnaname | grep -v ">" >> $outdir/RNAsequences.fna
				#echo "Present in Genbank"
			elif [[ $(cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.eukbact.filtered.csv | awk 'BEGIN{FS=OFS=","} {gsub(/\./, "", $1)} 1' | sed 's/Porospora_cf_gigantea_B/Porospora_gigantea/g' | grep $species | awk -F "," '{print $3}' | grep "EP" | wc -l) -gt 0 ]]
			then
				if [[ $(cat /Users/jmattic1/Documents/EukProt_included_data_sets.v03.2021_11_22.txt | grep $species | awk '{print $1}' | wc -l) -gt 0 ]]
				then
					if [[ $(cat /Users/jmattic1/Documents/EukProt_included_data_sets.v03.2021_11_22.txt | grep $species | awk '{print $1}' | wc -l) -gt 1 ]]
					then
						EPprotnum=$(echo $i | awk -F "_" '{print $3}')
						eukprot=$(cat /Users/jmattic1/Documents/EukProt_included_data_sets.v03.2021_11_22.txt | grep $species | grep $EPprotnum | awk '{print $1}')
					else
						eukprot=$(cat /Users/jmattic1/Documents/EukProt_included_data_sets.v03.2021_11_22.txt | grep $species | awk '{print $1}')
					fi
					if [[ $(\ls /Applications/Databases/EukProtTranscriptomes/assembled_transcriptomes/ | grep $eukprot | wc -l) -gt 0 ]]
					then
						RNAname=$(\ls /Applications/Databases/EukProtTranscriptomes/assembled_transcriptomes/ | grep $eukprot)
						RNAfile=$(\ls /Applications/Databases/EukProtTranscriptomes/assembled_transcriptomes/$RNAname)
						cat $RNAfile | sed 's/lcl\|//g' | awk '{if($0 ~ /^>/)  print $1"_"NR; else print $0;}' > $outdir/$species.$gene.rna.fna
						#cp $RNAfile $outdir/$species.$gene.rna.fna
						if [[ $(cat $outdir/$species.$gene.rna.fna | grep $gene | awk '{print $1}' | sed 's/>//g' | wc -l) -gt 0 ]]
						then
							rnaname=$(cat $outdir/$species.$gene.rna.fna | grep $gene | awk '{print $1}' | sed 's/>//g')
						else
							/Applications/miniconda3/bin/makeblastdb -in $outdir/$species.$gene.rna.fna -dbtype 'nucl' -parse_seqids
							samtools faidx $fasta $i > $outdir/protein.temp.faa
							/Applications/miniconda3/bin/tblastn -query $outdir/protein.temp.faa -db $outdir/$species.$gene.rna.fna -outfmt 6 -max_hsps 1 -num_threads 8 -out $outdir/$species.$gene.rna.blast -max_target_seqs 5129560 -evalue 1e-5
							if [[ $(cat $outdir/$species.$gene.rna.blast | head -1 | awk '{print $2}' | wc -l) -gt 0 ]]
							then
								rnaname=$(cat $outdir/$species.$gene.rna.blast | head -1 | awk '{print $2}')
							else
								echo "No RNA matches the protein... CHECK ME: "$i
								exit
							fi
							rm -f $outdir/protein.temp.faa
						fi
						echo ">"$i >> $outdir/RNAsequences.fna
						samtools faidx $outdir/$species.$gene.rna.fna $rnaname | grep -v ">" >> $outdir/RNAsequences.fna
					else
						echo "EukProt missing, go find it: "$i >> $outdir/eukprot.neededseqs.txt
						#exit
					fi
				else
					echo "EukProt mismatch: "$i
					exit
				fi
			elif [[ $( \ls /Applications/Databases/EukProtTranscriptomes/OtherEukTranscriptomes/*.fasta | grep $species | wc -l) -gt 0 ]]
			then
				RNAfile=$(\ls /Applications/Databases/EukProtTranscriptomes/OtherEukTranscriptomes/*.fasta | grep $species)
				cat $RNAfile | sed 's/lcl\|//g' | sed 's/\|/_/g' | awk '{if($0 ~ /^>/)  print $1"_"NR; else print $0;}' > $outdir/$species.$gene.rna.fna
				if [[ $(cat $outdir/$species.$gene.rna.fna | grep $gene | awk '{print $1}' | sed 's/>//g' | wc -l) -gt 0 ]]
				then
					rnaname=$(cat $outdir/$species.$gene.rna.fna | grep $gene | awk '{print $1}' | sed 's/>//g')
				else
					/Applications/miniconda3/bin/makeblastdb -in $outdir/$species.$gene.rna.fna -dbtype 'nucl' -parse_seqids
					samtools faidx $fasta $i > $outdir/protein.temp.faa
					/Applications/miniconda3/bin/tblastn -query $outdir/protein.temp.faa -db $outdir/$species.$gene.rna.fna -outfmt 6 -max_hsps 1 -num_threads 8 -out $outdir/$species.$gene.rna.blast -max_target_seqs 5129560 -evalue 1e-5
					if [[ $(cat $outdir/$species.$gene.rna.blast | head -1 | awk '{print $2}' | wc -l) -gt 0 ]]
					then
						rnaname=$(cat $outdir/$species.$gene.rna.blast | head -1 | awk '{print $2}')
					else
						echo "No RNA matches the protein... CHECK ME: "$i
						exit
					fi
					rm -f $outdir/protein.temp.faa
				fi
				echo ">"$i >> $outdir/RNAsequences.fna
				samtools faidx $outdir/$species.$gene.rna.fna $rnaname | grep -v ">" >> $outdir/RNAsequences.fna
			else
				echo "Non standard euks need to be updated: "$i
				exit
			fi
		else
			echo "Not fitting any character state, likely Laura Asgard: "$i
			exit
		fi
	done
	newoutdir=$(echo $wd"/RFAM_Classification/Infernal/")
	rm -rf $newoutdir
	mkdir $newoutdir
	/Applications/bin/cmscan -o $newoutdir/RNAAlignment.Infernal.txt --tblout $newoutdir/RNAAlignment.Infernal.tbl.tsv --cpu 8 --anytrunc /Applications/Rfam/CURRENT/Rfam.cm $outdir/RNAsequences.fna

fi


if [[ 1 -eq 2 ]]
then
	outdir=$(echo "/Users/jmattic1/Documents/TERT_reanalysis/RFAM_G2I_Classification/Infernal/")
	rm -rf $outdir
	mkdir $outdir
	cd $outdir
	cp /Users/jmattic1/Documents/ArchaeaGenbank/archaea/Candidatus_Heimdallarchaeota_archaeon_LC_3/latest_assembly_versions/GCA_001940645.1_ASM194064v1/GCA_001940645.1_ASM194064v1_rna_from_genomic.fna.gz $outdir
	cp /Users/jmattic1/Documents/ArchaeaGenbank/archaea/Candidatus_Heimdallarchaeota_archaeon_LC_3/latest_assembly_versions/GCA_001940645.1_ASM194064v1/GCA_001940645.1_ASM194064v1_genomic.fna.gz $outdir

	gzip -d GCA_001940645.1_ASM194064v1_rna_from_genomic.fna.gz
	gzip -d GCA_001940645.1_ASM194064v1_genomic.fna.gz
	/Applications/bin/cmpress /Applications/Rfam/CURRENT/Rfam.cm
	/Applications/bin/cmscan -o $outdir/RNAAlignment.Infernal.txt --tblout $outdir/RNAAlignment.Infernal.tbl.tsv --cpu 8 --anytrunc /Applications/Rfam/CURRENT/Rfam.cm /Users/jmattic1/Documents/TERT_reanalysis/RFAM_G2I_Classification/RNAsequences.fna
	/Applications/bin/cmscan -o $outdir/Heimdall.Infernal.genomic.txt --tblout $outdir/Heimdall.Infernal.genomic.tbl.tsv --cpu 8 --anytrunc /Applications/Rfam/CURRENT/Rfam.cm $outdir/GCA_001940645.1_ASM194064v1_genomic.fna


	samtools faidx /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I_e1e-10f30nh1ni6_408x290_5883//AsgardG2I.blast.1e-10.origin.nostar.faa CHeimdall_archaeon-L1_MCY3414279.1 > Heimdall.G2I.protein.faa
	/Applications/miniconda3/bin/makeblastdb -in GCA_001940645.1_ASM194064v1_genomic.fna -dbtype 'nucl' -parse_seqids
	/Applications/miniconda3/bin/tblastn -query Heimdall.G2I.protein.faa -db GCA_001940645.1_ASM194064v1_genomic.fna -outfmt 6 -max_hsps 1 -num_threads 8 -out Heimdall.G2I.protein.blast -max_target_seqs 5129560 -evalue 1e-5
fi


###Time to start over ... too many cross linked pathways
if [[ 1 -eq 2 ]]
then
	mkdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Introns.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 60 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Introns.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 30 -eval 1e-5 -filter 60 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 60 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 30 -eval 1e-5 -filter 60 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm


	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Intron_single.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 60 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Intron_single.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 30 -eval 1e-5 -filter 60 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_broad/RDP_Family_Separate_Analysis/GroupII_Intron_single.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm

	rm -rf /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/
	mkdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/Euk_ArchBact_Condensed_v3.1_fullasgard.faa -num_hits 1 -eval 1e-5 -filter 60 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/Euk_ArchBact_Condensed_v3.1_fullasgard.faa -num_hits 30 -eval 1e-5 -filter 60 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/Euk_ArchBact_Condensed_v3.1_fullasgard.faa -num_hits 1 -eval 1e-5 -filter 30 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom RVT_1 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/Euk_ArchBact_Condensed_v3.1_fullasgard.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom Intron_maturas2 -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	pfamdomain=$(echo "Intron_maturas2|GIIM")
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/Euk_ArchBact_Condensed_v3.1_fullasgard.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa -blastDB /Applications/Databases/ArchBactBLAST_DB_v2/CombinedBactArch.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	
	#Chloroplasts
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/ArabidopsisMatK.faa -blastDB /Applications/Databases/ArchBactBLAST_DB_v2/CombinedBactArch.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/ArabidopsisMatK.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/Euk_ArchBact_Condensed_v3.1_fullasgard.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm

	cat /Users/jmattic1/Documents/TERT_reanalysis/ArabidopsisMatK.faa /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa > /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.chloro.faa
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.chloro.faa -blastDB /Applications/Databases/ArchBactBLAST_DB_v2/CombinedBactArch.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	
	##Here we are.
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.chloro.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/Euk_ArchBact_Condensed_v3.1_fullasgard.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.chloro.faa -blastDB /Applications/Databases/Euk_ArchBact_Condensed_v3.1_fullasgard/Euk_ArchBact_Condensed_v3.1_fullasgard.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.chloro.faa -blastDB /Users/jmattic1/Documents/Orc_analysis/All_Asgard.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm


	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.chloro.faa -blastDB /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1/Viral_Euk_ArchBact_Condensed_v4.1.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	mkdir /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard/
	cat /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1/Viral_Euk_ArchBact_Condensed_v4.1.faa /Users/jmattic1/Documents/Orc_analysis/All_Asgard.faa > /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard.faa
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard.faa -dbtype 'prot' -parse_seqids
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.chloro.faa -blastDB /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.chloro.faa -blastDB /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -pfam_eval 1e-10
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.chloro.faa -blastDB /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard.faa -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm -pfam_eval 1e-20

	###Make Prp8 Tree
	rm -rf /Users/jmattic1/Documents/TERT_reanalysis/Prp8concat/
	mkdir /Users/jmattic1/Documents/TERT_reanalysis/Prp8concat/
	cat /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim.chloro_e1e-1f20nh20ni6_242x241_17305/AsgardG2I.heim.chloro.blast.1e-1.origin.nostar.faa /Users/jmattic1/Documents/Prp8.four.faa > /Users/jmattic1/Documents/TERT_reanalysis/Prp8concat/PFAM.1e-5Filter.allasgard.allbact.faa
	cat /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim.chloro_e1e-1f20nh20ni6_207x348_5531/AsgardG2I.heim.chloro.blast.1e-1.origin.nostar.faa /Users/jmattic1/Documents/Prp8.four.faa > /Users/jmattic1/Documents/TERT_reanalysis/Prp8concat/PFAM.1e-10Filter.allasgard.allbact.faa
	cat /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim.chloro_e1e-1f20nh20ni6_74x489_3230/AsgardG2I.heim.chloro.blast.1e-1.origin.nostar.faa /Users/jmattic1/Documents/Prp8.four.faa > /Users/jmattic1/Documents/TERT_reanalysis/Prp8concat/PFAM.1e-20Filter.allasgard.allbact.faa

	for i in $(ls /Users/jmattic1/Documents/TERT_reanalysis/Prp8concat/*.faa | grep "faa" | grep -v "align")
	do
		mafft=$(echo $i | sed 's/.faa/.align.faa/g')
		trim=$(echo $i | sed 's/.faa/.align.trim.faa/g')
		/Applications/miniconda3/bin/mafft --anysymbol --reorder --auto $i > $mafft
		/Applications/miniconda3/bin/trimal -in $mafft -out $trim -gappyout
		/Applications/miniconda3/bin/iqtree -T AUTO -s $trim --alrt 1000 -B 1000 -safe

	done

	##Bacterial G2I
	DB=/Users/jmattic1/Documents/BactBLAST_DB/Bact.faa
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/TERT_reanalysis/new_db_asgardqueries/AsgardG2I.heim.faa -blastDB $DB -num_hits 20 -eval 1e-1 -filter 20 -type psiblast -outdir /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/ -addorigin no -trimal -gappyout -pfam_dom $pfamdomain -pfamDB /Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm


	####PRP8 Single and multi Gene Tree
	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/Prp8AF/PRP8.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 1 -eval 1e-5 -filter 40 -type psiblast -outdir /Users/jmattic1/Documents/Prp8AF/ -addorigin no -trimal -gappyout 

	/Applications/BLAST_wrapper.sh -fasta /Users/jmattic1/Documents/Prp8AF/PRP8.faa -blastDB /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -num_hits 20 -eval 1e-5 -filter 40 -type psiblast -outdir /Users/jmattic1/Documents/Prp8AF/ -addorigin no -trimal -gappyout

fi

###Asgard DB RVTs
if [[ 1 -eq 1 ]]
then

	originfasta=$(\ls /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardDB_AsgardG2I.heim.chloro_e1e-1f20nh20ni6_81x383_17417/*origin.nostar.faa)

	wd=$(echo $originfasta | awk 'BEGIN{FS=OFS="/"}{NF--; print}')
	outdir=$(echo $wd"/RFAM_Classification/")
	rm -rf $outdir
	mkdir $outdir
	touch $outdir/RNAsequences.fna
	originfasta=$(\ls /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardDB_AsgardG2I.heim.chloro_e1e-1f20nh20ni6_81x383_17417/*origin.nostar.faa)
	cat $originfasta | sed 's/\//_/g' > $outdir/Blastwrapper.proteins.faa
	fasta=$(echo $outdir/Blastwrapper.proteins.faa)
	samtools faidx $fasta
	for i in $(cat $fasta | grep ">" | sed 's/>//g' | awk '{print $1}')
	do
		echo $i
		species=$(echo $i | awk -F "_" '{print $1"_"$2}')
		GCA=$(echo $i | awk -F "_" '{print $2}' | sed 's/-/_/g')

		gene=$(echo $i | cut -f 3- -d '_')
		###Regular Arch
			if [[ $(\ls /Users/jmattic1/Documents/Orc_analysis/AsgardGenomes/ncbi_dataset/data/$GCA/$GCA*.fna | wc -l) -gt 0 ]]
			then
				RNAfile=$(\ls /Users/jmattic1/Documents/Orc_analysis/AsgardGenomes/ncbi_dataset/data/$GCA/$GCA*.fna | head -1)

			else
				echo $species" is messed up, check it."
				break
			fi
			cp $RNAfile $outdir/$species.$gene.rna.fna
			if [[ $(cat $outdir/$species.$gene.rna.fna | grep $gene | awk '{print $1}' | sed 's/>//g' | wc -l) -gt 0 ]]
			then
				rnaname=$(cat $outdir/$species.$gene.rna.fna | grep $gene | awk '{print $1}' | sed 's/>//g')
			else
				/Applications/miniconda3/bin/makeblastdb -in $outdir/$species.$gene.rna.fna -dbtype 'nucl' -parse_seqids
				samtools faidx $fasta $i > $outdir/protein.temp.faa
				/Applications/miniconda3/bin/tblastn -query $outdir/protein.temp.faa -db $outdir/$species.$gene.rna.fna -outfmt 6 -max_hsps 1 -num_threads 8 -out $outdir/$species.$gene.rna.blast -max_target_seqs 5129560 -evalue 1e-5
				if [[ $(cat $outdir/$species.$gene.rna.blast | head -1 | awk '{print $2}' | wc -l) -gt 0 ]]
				then
					rnaname=$(cat $outdir/$species.$gene.rna.blast | head -1 | awk '{print $2}')
				else
					echo "No RNA matches the protein... CHECK ME: "$i
					#exit
				fi
				rm -f $outdir/protein.temp.faa
			fi
			echo ">"$i >> $outdir/RNAsequences.fna
			samtools faidx $outdir/$species.$gene.rna.fna $rnaname | grep -v ">" >> $outdir/RNAsequences.fna
		
	done
	newoutdir=$(echo $wd"/RFAM_Classification/Infernal/")
	rm -rf $newoutdir
	mkdir $newoutdir
	/Applications/bin/cmscan -o $newoutdir/RNAAlignment.Infernal.txt --tblout $newoutdir/RNAAlignment.Infernal.tbl.tsv --cpu 8 --anytrunc /Applications/Rfam/CURRENT/Rfam.cm $outdir/RNAsequences.fna

fi



######Intein stuff
if [[ 1 -eq 2 ]]
then
	for i in $(\ls /Users/jmattic1/Documents/InteinSearches/*.faa | grep "faa" | grep -v "DS2")
	do
		/Applications/BLAST_wrapper.sh -fasta $i -blastDB /Applications/Databases/ArchBactBLAST_DB_v2/CombinedBactArch.faa -num_hits 20 -eval 1e-1 -filter 95 -type psiblast -outdir /Users/jmattic1/Documents/InteinSearches/ArchBact/ -addorigin yes -trimal -automated1

	done

fi


if [[ 1 -eq 2 ]]
then
	endonucleasedomains=$(cat /Users/jmattic1/Documents/InteinSearches/*/InteinAnalysis/FullProt.dom.tsv | grep -E "endonuclease|LAGLIDADG" | awk '{print $1}' | sort | uniq | tr '\n' '|' | sed 's/\|$//g')
	for i in $(\ls /Users/jmattic1/Documents/InteinSearches/*/InteinAnalysis/FullProt.dom.tsv | grep "dom")
	do
		output=$(echo $i | sed 's/FullProt/InteinOnly/g')
		inteinlist=$(cat $i | awk '$1 == "Intein_splicing" {print $4}' | sort | uniq | tr '\n' '|' | sed 's/\|$//g')
		cat $i | grep -E "$inteinlist" | grep -v "#" > $output
		output2=$(echo $i | sed 's/FullProt/Intein_and_HED/g')
		output3=$(echo $i | sed 's/FullProt/Intein_and_HED.domonly/g')

		inteinlist2=$(cat $output | grep -E "$endonucleasedomains" | awk '{print $4}' | sort | uniq | tr '\n' '|' | sed 's/\|$//g')

		cat $output | sed 's/,//g' | grep -E "$inteinlist2" | awk '{print $1","$2","$3","$4","$5","$6","$7","$8","$9","$10","$11","$12","$2","$13","$14","$15","$16","$17","$18","$19","$20","$21","$22","$23"_"$24"_"$25"_"$26"_"$27"_"$28"_"$29"_"$30}' | sed 's/_*$//g' > $output2
		cat $output2 | grep -E "endonuclease|LAGLIDADG|Intein_splicing" | awk -F ',' '{print $1","$2","$3","$4","$5","$6","$7","$8","$9","$10","$11","$12","$2","$13","$14","$15","$16","$17","$18","$19","$20","$21","$22","$23"_"$24"_"$25"_"$26"_"$27"_"$28"_"$29"_"$30}' | sed 's/_*$//g' > $output3

	done


	for i in $(\ls /Users/jmattic1/Documents/InteinSearches/ | grep "t6b" | grep -v "faa")
	do
		perl /Applications/residues.pl /Users/jmattic1/Documents/InteinSearches/$i/*.origin.nostar.faa > /Users/jmattic1/Documents/InteinSearches/$i/InteinAnalysis/ProteinLengths.tsv
	done
fi

if [[ 1 -eq 2 ]]
then
	####new_db_asgardqueries/AsgardG2I_e1e-10f30nh1ni6_408x290_5883
	faafile=$(echo "/Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim_e1e-5f60nh1ni6_182x283_18474/AsgardG2I.heim.blast.1e-5.origin.nostar.faa")

	RNAfile=$(echo "/Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim_e1e-5f60nh1ni6_182x283_18474/RFAM_Classification/RNAsequences.fna")
	outdir=$(echo "/Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim_e1e-5f60nh1ni6_182x283_18474/")
	rm -rf $outdir/PFAM_Classification/
	mkdir $outdir/PFAM_Classification/
	pFAMDB="/Users/jmattic1/Documents/PFAM_DB/Pfam-A.hmm"
	/Applications/miniconda3/bin/hmmscan -o $outdir/PFAM_Classification/Prot.hmm.out --tblout $outdir/PFAM_Classification/Prot.hmm.tsv --domtblout $outdir/PFAM_Classification/Prot.dom.tsv --cpu 16 $pFAMDB $faafile
	#/Applications/miniconda3/bin/nhmmscan -o $outdir/PFAM_Classification/RNA.hmm.out --tblout $outdir/PFAM_Classification/RNA.hmm.tsv --dfamtblout $outdir/PFAM_Classification/RNA.dom.tsv --cpu 16 $pFAMDB $RNAfile
fi


####Internal G2I evaluation -- synteny
##Haloferax proteome against the genome -- blastx
if [[ 1 -eq 2 ]]
then
	####new_db_asgardqueries/AsgardG2I_e1e-10f30nh1ni6_408x290_5883
	faafile=$(echo "/Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim_e1e-1f20nh20ni6_138x500_10053/RFAM_Classification/RNAsequences.fna")
	wd=$(echo $faafile | awk 'BEGIN{FS=OFS="/"}{NF--; print}')
	grouptwointronregex=$(cat $wd/Infernal/RNAAlignment.Infernal.tbl.tsv | grep "intron" | awk '{print $1}' | sort | uniq | grep -E "group|Intron" | tr '\n' '|' | sed 's/\|$//g')
	cat $wd/Infernal/RNAAlignment.Infernal.tbl.tsv | grep -v "#" | sed 's/,/_/g' | awk '{print $1","$2","$3","$4","$5","$6","$7","$8","$9","$10","$11","$12","$13","$14","$15","$16","$17}' > $wd/Infernal/RNAAlignment.Infernal.tbl.csv
	rm -rf $wd/G2I_Synteny/
	mkdir $wd/G2I_Synteny/
	for i in $(cat $faafile | grep ">" | sed 's/>//g')
	do
		species=$(echo $i | awk -F "_" '{print $1"_"$2}')
		gene=$(echo $i | cut -f 3- -d '_')
		translatedspecies=$(cat /Users/jmattic1/Documents/FinishedFastas.txt | awk -v species="$species" '$2 == species".faa" {print $1}' | sed 's/.faa//g')
		if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/ | grep "latest_assembly_versions" | wc -l) -gt 0 ]]
		then
			if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/latest_assembly_versions/*/ | grep "gff.gz" | wc -l) -gt 0 ]]
			then
				proteinfile=$(ls -lt /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/latest_assembly_versions/*/*gff.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
			fi
			cp $proteinfile $wd/G2I_Synteny/$i.gff.gz
			gzip -d $wd/G2I_Synteny/$i.gff.gz
			contigname=$(cat $wd/$species.$gene.rna.blast | head -1 | awk '{print $2}')
			cat $wd/G2I_Synteny/$i.gff | grep -v "#" | awk -v filename="$i" -v contigname="$contigname" '$1 == contigname {print $0"\t"filename}' | sed 's/,/_/g' | sed "s/\'/prime/g" | tr '\t' ',' >> $wd/G2I_Synteny/All.region.gff.csv
		elif [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/ | grep "all_assembly_versions" | wc -l) -gt 0 ]]
		then
			if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/*/ | grep "gff.gz" | wc -l) -gt 0 ]]
			then
				proteinfile=$(ls -lt /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/*/*gff.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
			elif [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/ | grep "suppressed" | wc -l) -gt 0 ]]
			then
				if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/suppressed/*/ | grep "gff.gz" | wc -l) -gt 0 ]]
				then
					proteinfile=$(ls -lt /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$translatedspecies/all_assembly_versions/suppressed/*/*gff.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
				fi
			fi
			cp $proteinfile $wd/G2I_Synteny/$i.gff.gz
			gzip -d $wd/G2I_Synteny/$i.gff.gz
			contigname=$(cat $wd/$species.$gene.rna.blast | head -1 | awk '{print $2}')
			cat $wd/G2I_Synteny/$i.gff | grep -v "#" | awk -v filename="$i" -v contigname="$contigname" '$1 == contigname {print $0"\t"filename}' | sed 's/,/_/g' | sed "s/\'/prime/g" | tr '\t' ',' >> $wd/G2I_Synteny/All.region.gff.csv
		else

			echo "Doesn't have a gff file: "$i 
		fi
	done
fi

####PRP8 Alignments with G2Is

if [[ 1 -eq 2 ]]
then
	conda activate /Applications/miniconda3/envs/samtools
	cd /Users/jmattic1/Documents/Prp8AF/
	cat Group2Intron.hits.16504.faa > Group2Intron.Prp8.combined.faa
	for i in $(cat PRP8_e1e-5f40nh1ni6_182x2293_25360/PRP8.blast.1e-5.origin.nostar.faa | grep -E "Homo|Nebula|Arab|Symb" | sed 's/>//g')
	do
		samtools faidx PRP8_e1e-5f40nh1ni6_182x2293_25360/PRP8.blast.1e-5.origin.nostar.faa $i >> Group2Intron.Prp8.combined.faa
	done

	for i in $(cat /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim_e1e-1f20nh20ni6_138x500_10053/AsgardG2I.heim.blast.1e-1.origin.nostar.faa | grep -E "Heim|Loki" | grep -v "CLokiarchaeum_75" | sed 's/>//g')
	do
		samtools faidx /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim_e1e-1f20nh20ni6_138x500_10053/AsgardG2I.heim.blast.1e-1.origin.nostar.faa $i >> Group2Intron.Prp8.combined.faa
	done

	/Applications/miniconda3/bin/mafft --anysymbol --reorder --auto Group2Intron.Prp8.combined.faa > Group2Intron.Prp8.combined.align.faa
	/Applications/miniconda3/bin/trimal -in Group2Intron.Prp8.combined.align.faa -out Group2Intron.Prp8.combined.align.trimal.faa -gappyout
	/Applications/miniconda3/bin/iqtree -T AUTO -s Group2Intron.Prp8.combined.align.trimal.faa --alrt 1000 -B 1000
fi
if [[ 1 -eq 2 ]]
then
	####RNA Alignments
	Qfasta=$(\ls /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim_e1e-1f20nh20ni6_138x500_10053/RFAM_Classification//RNAsequences.fna)
	outdir=$(echo "/Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim_e1e-1f20nh20ni6_138x500_10053/RFAM_Classification/Infernal/")
	#conda activate /Applications/miniconda3/envs/samtools
	samtools faidx $Qfasta
	Rfamqueries=$(cat /Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim_e1e-1f20nh20ni6_138x500_10053/RFAM_Classification/Infernal/RNAAlignment.Infernal.tbl.tsv | grep -E "group-II|gpII" | awk '{print $1}' | sort | uniq | tr '\n' '|')
	rm -rf $outdir/Alignment/
	mkdir $outdir/Alignment/
	cat $outdir/RNAAlignment.Infernal.tbl.tsv | grep -E "group-II|gpII" > $outdir/Alignment/G2I_Subset_RFAM.tsv
	#/Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim_e1e-1f20nh20ni6_138x500_10053/RFAM_Classification/Infernal/RNAAlignment.Infernal.tbl.tsv 
	while read line
	do
		g2Iclass=$(echo $line | awk '{print $1}')
		ori=$(echo $line | awk '{print $10}')
		contig=$(echo $line | awk '{print $3}')
		
		if [[ $ori == "-" ]]
		then
			start=$(echo $line | awk '{print $9}')
			end=$(echo $line | awk '{print $8}')
			samtools faidx $Qfasta $contig:$start-$end > $outdir/Alignment/subFasta.fna
			cat $outdir/Alignment/subFasta.fna | grep ">" > $outdir/Alignment/subFasta.rev.fna
			cat $outdir/Alignment/subFasta.fna | grep -v ">" | rev | tr 'actgACTG' 'tgacTGAC' >> $outdir/Alignment/subFasta.rev.fna
			cat $outdir/Alignment/subFasta.rev.fna > $outdir/Alignment/subFasta.fna
		else
			start=$(echo $line | awk '{print $8}')
			end=$(echo $line | awk '{print $9}')
			samtools faidx $Qfasta $contig:$start-$end > $outdir/Alignment/subFasta.fna
		fi
		cat $outdir/Alignment/subFasta.fna >> $outdir/Alignment/$g2Iclass.fasta
		rm $outdir/Alignment/subFasta*

	done < $outdir/Alignment/G2I_Subset_RFAM.tsv
	cat $outdir/Alignment/*.fasta > $outdir/Alignment/All.fasta
	rm -rf ~/Documents/TurboFold/
	mkdir ~/Documents/TurboFold/
	for i in $(\ls $outdir/Alignment/*.fasta)
	do
		outfile=$(echo $i | awk -F '/' '{print $NF}' | sed 's/.fasta//g')
		cat <(cat ~/Downloads/RF00020.fasta | awk '{print $1}' | sed 's/\//_/g') <(cat $i | sed 's/:/_/g') | sed -E '/>/!s/R/A/g' | sed -E '/>/!s/Y/T/g' | sed -E '/>/!s/W/A/g' > ~/Documents/TurboFold/$outfile.U5.dup.fasta
		seqkit rmdup -s < ~/Documents/TurboFold/$outfile.U5.dup.fasta > ~/Documents/TurboFold/$outfile.U5.nodupfasta
		cat ~/Documents/TurboFold/$outfile.U5.nodupfasta | awk '{if($0 ~ /^>/)  print $0"_"NR; else print $0;}' >  ~/Documents/TurboFold/$outfile.U5.fasta
		/Applications/miniconda3/bin/mafft --anysymbol --reorder --auto ~/Documents/TurboFold/$outfile.U5.fasta > ~/Documents/TurboFold/$outfile.U5.align.fasta
		/Applications/miniconda3/bin/trimal -in ~/Documents/TurboFold/$outfile.U5.align.fasta -out ~/Documents/TurboFold/$outfile.U5.align.trimal.fasta -gappyout
		/Applications/miniconda3/bin/iqtree -T AUTO -s ~/Documents/TurboFold/$outfile.U5.align.trimal.fasta --alrt 1000 -B 1000
	done

	#####Excel sheet generator
	cat //Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/LauraIncluded/AsgardG2I.heim_e1e-1f20nh20ni6_138x500_10053/RFAM_Classification/Infernal/Alignment/G2I_Subset_RFAM.tsv | grep "Intron_gpII" | awk '$15 > 40 {print $3}' | sort | uniq -c | awk '{print $2"\t"$1}' > ./RDP_Paper_Nature/RFAM_incidence.v0.tsv
fi


