if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/ArchRepresentativeTaxa/
	mkdir /Users/jmattic1/Documents/ArchRepresentativeTaxa/
	grep -l "Esch" /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_automated/*/*trimal.faa | awk -F '/' '{print $6}' > /Users/jmattic1/Documents/ArchRepresentativeTaxa/FullBactProteins.txt

	while read line
	do
		/Applications/Alignment_refiner.sh -blastoutput /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_automated/$line -cutoff 0.95 -db /Applications/Databases/ArchBactBLAST_DB/CombinedBactArch.faa -trimal -automated1

	done < /Users/jmattic1/Documents/ArchRepresentativeTaxa/FullBactProteins.txt

	rm -rf /Users/jmattic1/Documents/ArchRepresentativeTaxa/RPL8_subsettree/
	mkdir /Users/jmattic1/Documents/ArchRepresentativeTaxa/RPL8_subsettree/

	cat <(cat /Users/jmattic1/Documents/Archaealfamilies.0.5.tsv | awk '{print $1}') <(cat /Users/jmattic1/Documents/RPL8.HVolc.blast.1e-5.align.trimal.faa.mldist.tab.txt | awk '{print $1}' | grep -E "Asgard|Heim|Thor|Odin|Loki|Hel") | grep -v "Representative" | sed 's/Bacillus_s_NP_388000.2_ribosomal_protein_L2__$/Bacillus_s_NP_388000.2_ribosomal_protein_L2_( /g' | sort | uniq > /Users/jmattic1/Documents/ArchRepresentativeTaxa/RPL8_subsettree/proteinlist.txt

	while read line
	do
		proteinname=$(cat /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_automated//RPL8.HVolc_e1e-5f0.25nh1ni6_1054x/RPL8.HVolc_BLAST_27377_running/RPL8.HVolc.blast.1e-5.origin.faa | grep $line | sed 's/>//g')
		echo $proteinname >> /Users/jmattic1/Documents/ArchRepresentativeTaxa/RPL8_subsettree/actualproteinlist.txt
		echo $line
		echo $proteinname
	done < /Users/jmattic1/Documents/ArchRepresentativeTaxa/RPL8_subsettree/proteinlist.txt

	faidxproteins=$(cat /Users/jmattic1/Documents/ArchRepresentativeTaxa/RPL8_subsettree/actualproteinlist.txt | tr '\n' ' ')

	#samtools faidx /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_automated//RPL8.HVolc_e1e-5f0.25nh1ni6_1054x/RPL8.HVolc_BLAST_27377_running/RPL8.HVolc.blast.1e-5.origin.faa $faidxproteins  /Users/jmattic1/Documents/ArchRepresentativeTaxa/RPL8_subsettree/RPL8.subset.faa>
	/Applications/miniconda3/bin/mafft --reorder --auto /Users/jmattic1/Documents/ArchRepresentativeTaxa/RPL8_subsettree/RPL8.subset.faa > /Users/jmattic1/Documents/ArchRepresentativeTaxa/RPL8_subsettree/RPL8.subset.align.faa
	/Applications/miniconda3/bin/trimal -in /Users/jmattic1/Documents/ArchRepresentativeTaxa/RPL8_subsettree/RPL8.subset.align.faa -out /Users/jmattic1/Documents/ArchRepresentativeTaxa/RPL8_subsettree/RPL8.subset.align.trimal.faa -automated1
	/Applications/miniconda3/bin/iqtree -T AUTO -s /Users/jmattic1/Documents/ArchRepresentativeTaxa/RPL8_subsettree/RPL8.subset.align.trimal.faa --alrt 1000 -B 1000
fi


##################Using MLDist to subset. 
if [[ 1 -eq 2 ]]
then
	##Make the master fasta file
	rm /Users/jmattic1/Documents/AllEcoliContainingRBPs.faa
	while read line
	do
		cat /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_automated/$line/*trimal.faa >> /Users/jmattic1/Documents/AllEcoliContainingRBPs.faa
	done < /Users/jmattic1/Documents/ArchRepresentativeTaxa/FullBactProteins.txt

	###Check all files for family members......
	cat /Users/jmattic1/Documents/Archaealfamilies.0.3.tsv | tail -n+2 > /Users/jmattic1/Documents/Archaealfamilies.0.5.nohead.tsv
	rm -f /Users/jmattic1/Documents/ArchBactCondensed.txt
	rm -f /Users/jmattic1/Documents/ArchBactCondensed.faa
	rm -f /Users/jmattic1/Documents/ArchBactCondensed.faa.gz
	rm -f /Users/jmattic1/Documents/RemovedTaxa.txt
	rm -f /Users/jmattic1/Documents/KeptTaxa.txt
	touch /Users/jmattic1/Documents/KeptTaxa.txt

	while read line
	do
		halopresence=$(echo $line | grep "Haloferax")
		ecolipresence=$(echo $line | grep "Escherichia_c")
		sulfpresence=$(echo $line | grep "acidocaldarius")

		if [[ $halopresence != "" ]]
		then
			representative=$(echo "Haloferax_volcanii")
		else
			representative=$(echo $line | awk '{print $1}' | awk -F "_" '{print $1"_"$2}')
		fi
		if [[ $ecolipresence != "" ]]
		then
			representative=$(echo "Escherichia_c")
		fi
		if [[ $sulfpresence != "" ]]
		then
			representative=$(echo "Sulfolobus_acidocaldarius")
		fi
		
		numhits=$(cat /Users/jmattic1/Documents/AllEcoliContainingRBPs.faa | grep $representative | wc -l)

		if [[ $numhits -ne 27 ]]
		then
			count=1
			totalcount=$(echo $line | awk '{print $2}' | awk -F ',' '{print NF}')
			while [ "$count" -le $totalcount ] && [ "$numhits" -ne 27 ]
			do
				currentcheck=$(echo $line | awk '{print $2}' | awk -v count="$count" -F ',' '{print $count}')
				currentrepresentative=$(echo $currentcheck | awk '{print $1}' | awk -F "_" '{print $1"_"$2}')
				numhits=$(cat /Users/jmattic1/Documents/AllEcoliContainingRBPs.faa | grep $currentrepresentative | wc -l)
				count=$(( $count + 1 ))
			done
			if [[ $numhits -ne 27 ]]
			then
				echo "No hits found, removing"
				echo $representative >> /Users/jmattic1/Documents/RemovedTaxa.txt
			else
				representative=$(echo $currentrepresentative)
				presence=$(cat /Users/jmattic1/Documents/KeptTaxa.txt | grep $representative)
				if [[ $presence == "" ]]
				then
					cat /Applications/Databases/ArchBactBLAST_DB_v2/CombinedBactArch.faa | grep $representative >> /Users/jmattic1/Documents/ArchBactCondensed.txt
					faidxgrab=$(cat /Applications/Databases/ArchBactBLAST_DB_v2/CombinedBactArch.faa | grep $representative | sed 's/>//g' | tr '\n' ' ')
					samtools faidx /Applications/Databases/ArchBactBLAST_DB_v2/CombinedBactArch.faa $faidxgrab >> /Users/jmattic1/Documents/ArchBactCondensed.faa
					echo $representative >> /Users/jmattic1/Documents/KeptTaxa.txt
				else
					echo $representative" was duplicated"
				fi
			fi
		else
			presence=$(cat /Users/jmattic1/Documents/KeptTaxa.txt | grep $representative)
			if [[ $presence == "" ]]
			then
				cat /Applications/Databases/ArchBactBLAST_DB_v2/CombinedBactArch.faa | grep $representative >> /Users/jmattic1/Documents/ArchBactCondensed.txt
				faidxgrab=$(cat /Applications/Databases/ArchBactBLAST_DB_v2/CombinedBactArch.faa | grep $representative | sed 's/>//g' | tr '\n' ' ')
				samtools faidx /Applications/Databases/ArchBactBLAST_DB_v2/CombinedBactArch.faa $faidxgrab >> /Users/jmattic1/Documents/ArchBactCondensed.faa
				echo $representative >> /Users/jmattic1/Documents/KeptTaxa.txt
			else
				echo $representative" was duplicated"
			fi
		fi


	done < /Users/jmattic1/Documents/Archaealfamilies.0.5.nohead.tsv

	rm -rf /Applications/Databases/ArchBact_Condensed/
	mkdir /Applications/Databases/ArchBact_Condensed/
	cp /Users/jmattic1/Documents/ArchBactCondensed.faa /Applications/Databases/ArchBact_Condensed/
	gzip /Users/jmattic1/Documents/ArchBactCondensed.faa
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/ArchBact_Condensed/ArchBactCondensed.faa -dbtype 'prot' -parse_seqids

fi


##################Unsubsetted concatenation w/ Lauras Asgards -- include Hod and E. coli
if [[ 1 -eq 2 ]]
then
	##Make the master fasta file
	#conda activate samtools
	rm -rf /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/
	mkdir /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/
	cat <(grep -l "Deinococcus_r" /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_combinedarchbact_extraasgard/*/*trimal.faa | awk -F '/' '{print $6}') <(grep -l "Bdellovibrio" /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_combinedarchbact_extraasgard/*/*trimal.faa | awk -F '/' '{print $6}') <(grep -l "Esch" /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_combinedarchbact_extraasgard/*/*trimal.faa | awk -F '/' '{print $6}') <(grep -l "Hod" /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_combinedarchbact_extraasgard/*/*trimal.faa | awk -F '/' '{print $6}') | sort | uniq -c | awk '$1 > 3 {print $2}' > /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/FullBactProteins.txt

	mkdir /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/HodContainingRPs/
	for i in $(cat /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/FullBactProteins.txt | grep "RP")
	do
		cp -r /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_combinedarchbact_extraasgard/$i /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/HodContainingRPs/
	done

	cat /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/HodContainingRPs/*/*.trimal.faa | grep ">" | sed 's/>//g' | awk -F "_" '{print $1"_"$2}' | sort | uniq -c | awk '$1 == 19 {print $2}' > /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/IncludedSpecies.txt

	searchcriteria=$(cat /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/IncludedSpecies.txt | tr '\n' '|' | sed 's/|$//g')


	while read line
	do
		for i in $(cat /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/FullBactProteins.txt | grep "RP")
		do
			fasta=$(\ls /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/HodContainingRPs/$i/*.align.trimal.faa)

			extractiongene=$(cat $fasta | grep $line"_" | sed 's/>//g')
			samtools faidx $fasta $extractiongene >> /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/Unconcat.faa
		done
		echo ">"$line > /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/IndvConcat.faa
		cat /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/Unconcat.faa | grep -v ">" | tr -d "\n" >> /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/IndvConcat.faa
		rm -f /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/Unconcat.faa
		fold -w 80 /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/IndvConcat.faa >> /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/RibosomalProteins.concatenated.uncondensed.faa
	done < /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/IncludedSpecies.txt
	cat /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/RibosomalProteins.concatenated.uncondensed.faa | sed 's/>/\n>/g' | sed '/^$/d' > /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/RibosomalProteins.concatenated.uncondensed.final.faa
	/Applications/miniconda3/bin/iqtree -T AUTO -s /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/RibosomalProteins.concatenated.uncondensed.final.faa --alrt 1000 -B 1000 -wsr
	while read line
	do
		i=RPL8.HVolc_e1e-5f25nh1ni6_1551x
		fasta=$(\ls /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/HodContainingRPs/$i/*.align.trimal.faa)

		extractiongene=$(cat $fasta | grep $line"_" | sed 's/>//g')
		samtools faidx $fasta $extractiongene >> /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/Unconcat.faa
		echo ">"$line > /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/IndvConcat.faa
		cat /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/Unconcat.faa | grep -v ">" | tr -d "\n" >> /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/IndvConcat.faa
		rm -f /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/Unconcat.faa
		fold -w 80 /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/IndvConcat.faa >> /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/RPL8.concatenated.uncondensed.faa
	done < /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/IncludedSpecies.txt
	cat /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/RPL8.concatenated.uncondensed.faa | sed 's/>/\n>/g' | sed '/^$/d' > /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/RPL8.concatenated.uncondensed.final.faa
	/Applications/miniconda3/bin/iqtree -T AUTO -s /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/RPL8.concatenated.uncondensed.final.faa --alrt 1000 -B 1000 -wsr

	##Run R script
	mamba activate fisher
	/Applications/aa_recoder.py -i RibosomalProteins.concatenated.uncondensed.final.faa -re SR4
	/Applications/miniconda3/bin/iqtree -T AUTO -s /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/RibosomalProteins.concatenated.uncondensed.SR4recode.final.faa --alrt 1000 -B 1000 -wsr



	finaltable=$(echo "RibosomalProteinConcat.basalnode.1.tsv")
	withdrawallist=$(tail -n+2 $finaltable | awk '{print $1}' | tr '\n' ' ')
	samtools faidx /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/RibosomalProteins.concatenated.uncondensed.final.faa $withdrawallist > /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/RibosomalProteins.concatenated.subset.final.faa
	/Applications/miniconda3/bin/iqtree -T AUTO -s /Users/jmattic1/Documents/ArchRepresentativeTaxa_v2/RibosomalProteins.concatenated.subset.final.faa --alrt 1000 -B 1000 -wsr

	rm -rf /Applications/Databases/ArchBact_Condensed_v2/
	mkdir /Applications/Databases/ArchBact_Condensed_v2/
	for i in $(tail -n+2 $finaltable | awk '{print $1}')
	do
		cat /Applications/Databases/ArchBactBLAST_DB_v2_xAsgard/CombinedBactArchxAsgard.faa | grep $i"_" >> /Applications/Databases/ArchBact_Condensed_v2/ArchBact_Condensed_v2.faa

	done
fi

##################Make the paper tree vvvvvvvvvvvvvvvvvvvv
if [[ 1 -eq 2 ]]
then
	##Make the master fasta file
	#conda activate samtools
	rm -rf /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/
	mkdir /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/
	cat <(grep -l "Thermus_a" /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_reducedset_v3.1/*/*trimal.faa | awk -F '/' '{print $6}') <(grep -l "Bdellovibrio" /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_reducedset_v3.1/*/*trimal.faa | awk -F '/' '{print $6}') <(grep -l "Esch" /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_reducedset_v3.1/*/*trimal.faa | awk -F '/' '{print $6}') <(grep -l "Hod" /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_reducedset_v3.1/*/*trimal.faa | awk -F '/' '{print $6}') | sort | uniq -c | awk '$1 > 3 {print $2}' > /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/FullBactProteins.txt

	mkdir /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/HodContainingRPs/
	for i in $(cat /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/FullBactProteins.txt | grep "RP")
	do
		cp -r /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_reducedset_v3.1/$i /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/HodContainingRPs/
	done

	cat /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/HodContainingRPs/*/*.trimal.faa | grep ">" | sed 's/>//g' | awk -F "_" '{print $1"_"$2}' | sort | uniq -c | awk '$1 == 20 {print $2}' > /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/IncludedSpecies.txt

	searchcriteria=$(cat /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/IncludedSpecies.txt | tr '\n' '|' | sed 's/|$//g')


	while read line
	do
		for i in $(cat /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/FullBactProteins.txt | grep "RP")
		do
			fasta=$(\ls /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/HodContainingRPs/$i/*.align.trimal.faa)

			extractiongene=$(cat $fasta | grep $line"_" | sed 's/>//g')
			samtools faidx $fasta $extractiongene >> /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/Unconcat.faa
		done
		echo ">"$line > /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/IndvConcat.faa
		cat /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/Unconcat.faa | grep -v ">" | tr -d "\n" >> /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/IndvConcat.faa
		rm -f /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/Unconcat.faa
		fold -w 80 /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/IndvConcat.faa >> /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/RibosomalProteins.concatenated.uncondensed.faa
	done < /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/IncludedSpecies.txt
	cat /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/RibosomalProteins.concatenated.uncondensed.faa | sed 's/>/\n>/g' | sed '/^$/d' > /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/RibosomalProteins.concatenated.final.faa
	/Applications/miniconda3/bin/iqtree -T AUTO -s /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/RibosomalProteins.concatenated.final.faa --alrt 1000 -B 1000 -wsr
	
	/Applications/aa_recoder.py -i /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/RibosomalProteins.concatenated.final.faa -re SR4 -o /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/RibosomalProteins.concatenated.final.SR4
	/Applications/miniconda3/bin/iqtree -T AUTO -s /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/RibosomalProteins.concatenated.final.SR4.fas --alrt 1000 -B 1000 -wsr



	finaltable=$(echo "RibosomalProteinConcat.basalnode.1.tsv")
	withdrawallist=$(tail -n+2 $finaltable | awk '{print $1}' | tr '\n' ' ')
	samtools faidx /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/RibosomalProteins.concatenated.uncondensed.final.faa $withdrawallist > /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/RibosomalProteins.concatenated.subset.final.faa
	/Applications/miniconda3/bin/iqtree -T AUTO -s /Users/jmattic1/Documents/RibosomalProtein_ConcatTree/RibosomalProteins.concatenated.subset.final.faa --alrt 1000 -B 1000 -wsr

	rm -rf /Applications/Databases/ArchBact_Condensed_v2/
	mkdir /Applications/Databases/ArchBact_Condensed_v2/
	for i in $(tail -n+2 $finaltable | awk '{print $1}')
	do
		cat /Applications/Databases/ArchBactBLAST_DB_v2_xAsgard/CombinedBactArchxAsgard.faa | grep $i"_" >> /Applications/Databases/ArchBact_Condensed_v2/ArchBact_Condensed_v2.faa

	done
fi



###Make each individual alignment
if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/
	mkdir /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/
	while read gene
	do
		name=$(echo $gene | awk -F "_" '{print $1}')
		greptaxa=$(cat /Users/jmattic1/Documents/KeptTaxa.txt | tr '\n' '|' | sed 's/\|$//g')
		faidxgenes=$(cat /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_automated/$gene/*origin.faa | grep -E $greptaxa | sed 's/>//g' | tr '\n' ' ')
		samtools faidx /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_automated/$gene/*origin.faa $faidxgenes > /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/$name.faa
		/Applications/miniconda3/bin/mafft --auto /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/$name.faa > /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/$name.align.faa
		/Applications/miniconda3/bin/trimal -in /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/$name.align.faa -out /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/$name.align.trimmed.faa -automated1

	done < /Users/jmattic1/Documents/ArchRepresentativeTaxa/FullBactProteins.txt


fi

###Make concatenated alignment
if [[ 1 -eq 2 ]]
then
	rm /Users/jmattic1/Documents/ArchRepresentativeTaxa/ConcatenatedRibosomeProteins.aligned.trimmed.faa
	cat /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/*.align.trimmed.faa > /Users/jmattic1/Documents/ArchRepresentativeTaxa/AllTaxa.faa
	samtools faidx /Users/jmattic1/Documents/ArchRepresentativeTaxa/AllTaxa.faa
	rm -f /Users/jmattic1/Documents/ArchRepresentativeTaxa/ConcatenatedRibosomeProteins*
	while read species
	do
		seqnumber=$(cat /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/*.align.trimmed.faa | grep $species | sort | uniq | wc -l)
		if [[ $seqnumber -ne 27 ]]
		then
			echo $species" has less than 27 sequences"
			break
		fi
		genefaidx=$(cat /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/*.align.trimmed.faa | grep $species | sed 's/>//g' | tr '\n' ' ')
		echo ">"$species >> /Users/jmattic1/Documents/ArchRepresentativeTaxa/ConcatenatedRibosomeProteins.aligned.trimmed.faa
		samtools faidx /Users/jmattic1/Documents/ArchRepresentativeTaxa/AllTaxa.faa $genefaidx | grep -v ">" | tr -d '\n' | awk '{print $0}' >> /Users/jmattic1/Documents/ArchRepresentativeTaxa/ConcatenatedRibosomeProteins.aligned.trimmed.faa
	done < /Users/jmattic1/Documents/KeptTaxa.txt

	fold -w 80 /Users/jmattic1/Documents/ArchRepresentativeTaxa/ConcatenatedRibosomeProteins.aligned.trimmed.faa > /Users/jmattic1/Documents/ArchRepresentativeTaxa/ConcatenatedRibosomeProteins.aligned.trimmed.folded.faa
	/Applications/miniconda3/bin/iqtree -T AUTO -s /Users/jmattic1/Documents/ArchRepresentativeTaxa/ConcatenatedRibosomeProteins.aligned.trimmed.folded.faa --alrt 1000 -B 1000

fi



###Make concatenated alignment of untrimmed
if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/ArchRepresentativeTaxa/Untrimmed/
	mkdir /Users/jmattic1/Documents/ArchRepresentativeTaxa/Untrimmed/

	cat /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/*.align.faa > /Users/jmattic1/Documents/ArchRepresentativeTaxa/Untrimmed/AllTaxa.faa
	samtools faidx /Users/jmattic1/Documents/ArchRepresentativeTaxa/Untrimmed/AllTaxa.faa
	while read species
	do
		seqnumber=$(cat /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/*.align.faa | grep $species | sort | uniq | wc -l)
		if [[ $seqnumber -ne 27 ]]
		then
			echo $species" has less than 27 sequences"
			break
		fi
		genefaidx=$(cat /Users/jmattic1/Documents/ArchRepresentativeTaxa/FinalAlignments/*.align.faa | grep $species | sed 's/>//g' | tr '\n' ' ')
		echo ">"$species >> /Users/jmattic1/Documents/ArchRepresentativeTaxa/Untrimmed/ConcatenatedRibosomeProteins.aligned.faa
		samtools faidx /Users/jmattic1/Documents/ArchRepresentativeTaxa/Untrimmed/AllTaxa.faa $genefaidx | grep -v ">" | tr -d '\n' | awk '{print $0}' >> /Users/jmattic1/Documents/ArchRepresentativeTaxa/Untrimmed/ConcatenatedRibosomeProteins.aligned.faa
	done < /Users/jmattic1/Documents/KeptTaxa.txt

	fold -w 80 /Users/jmattic1/Documents/ArchRepresentativeTaxa/Untrimmed/ConcatenatedRibosomeProteins.aligned.faa > /Users/jmattic1/Documents/ArchRepresentativeTaxa/Untrimmed/ConcatenatedRibosomeProteins.aligned.folded.faa
	#/Applications/miniconda3/bin/iqtree -T AUTO -s /Users/jmattic1/Documents/ArchRepresentativeTaxa/ConcatenatedRibosomeProteins.aligned.trimmed.folded.faa --alrt 1000 -B 1000

fi

###Make the euk arch bact tree
if [[ 1 -eq 2 ]]
then
	#outdir="/Users/jmattic1/Documents/EukArchBactTree/"
	outdir="/Users/jmattic1/Documents/EukArchBactTree_automatedxAsgard/"
	rm -rf $outdir
	mkdir $outdir
	#aligndir="/Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_eukarchbact_automated/"
	#aligndir="/Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_eukarchbact_strict/"
	aligndir="/Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_eukarchbact_automated_extraAsgard/"

	cat <(grep -l "Esch" $aligndir/*/*trimal.faa | awk -F '/' '{print $NF}') <(grep -l "Homo_s" $aligndir/*/*trimal.faa | awk -F '/' '{print $NF}') | sort | uniq -c | awk '$1 > 1 {print $2}' > $outdir/FullBactEukProteins.txt

	mkdir $outdir/Alignments
	while read line
	do
		cat $aligndir/*/$line >> $outdir/AllEcoliHomoContainingRBPs.faa
		cp $aligndir/*/$line $outdir/Alignments/
	done < $outdir/FullBactEukProteins.txt

	cat $outdir/AllEcoliHomoContainingRBPs.faa | grep ">" | sed 's/>//g' | awk -F "_" '{print $1"_"$2}' | sort | uniq -c | awk '$1 == 28 {print $2}' > $outdir/IncludedSpecies.txt
	cat $outdir/AllEcoliHomoContainingRBPs.faa | grep ">" | sed 's/>//g' | awk -F "_" '{print $1"_"$2}' | sort | uniq -c | awk '$1 < 28 {print $2}' > $outdir/NonIncludedSpecies.txt
	samtools faidx $outdir/AllEcoliHomoContainingRBPs.faa
	mkdir $outdir/ConcatTrimmed/
	while read species
	do
		seqnumber=$(cat $outdir/Alignments/*.align.trimal.faa | grep $species | sort | uniq | wc -l)
		if [[ $seqnumber -ne 28 ]]
		then
			echo $species" has less than 28 sequences"
			break
		fi
		genefaidx=$(cat $outdir/Alignments/*.align.trimal.faa | grep $species | sed 's/>//g' | tr '\n' ' ')
		echo ">"$species >> $outdir/ConcatTrimmed/ConcatenatedRibosomeProteins.aligned.trimmed.faa
		samtools faidx $outdir/AllEcoliHomoContainingRBPs.faa $genefaidx | grep -v ">" | tr -d '\n' | awk '{print $0}' >> $outdir/ConcatTrimmed/ConcatenatedRibosomeProteins.aligned.trimmed.faa
	done < $outdir/IncludedSpecies.txt

	fold -w 80 $outdir/ConcatTrimmed/ConcatenatedRibosomeProteins.aligned.trimmed.faa > $outdir/ConcatTrimmed/ConcatenatedRibosomeProteins.aligned.trimmed.folded.faa
	/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir/ConcatTrimmed/ConcatenatedRibosomeProteins.aligned.trimmed.folded.faa --alrt 1000 -B 1000 -wsr


fi


###Add lauras asgard and build new DB
if [[ 1 -eq 2 ]]
then
	rm -rf /Applications/Databases/EukArchBact_xAsgard/
	mkdir /Applications/Databases/EukArchBact_xAsgard/

	rm -rf /Applications/Databases/EukArchBactViral_xAsgard/
	mkdir /Applications/Databases/EukArchBactViral_xAsgard/

	rm -rf /Applications/Databases/EukArchBactViral_xAsgard/
	mkdir /Applications/Databases/ArchBactBLAST_DB_v2_xAsgard/

	cat /Applications/Databases/EukArchBact/EukArchBact.faa > /Applications/Databases/EukArchBact_xAsgard/EukArchBact_xAsgard.faa
	cat /Applications/Databases/EukArchBactViral/EukArchBactViral.faa > /Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa
	cat /Applications/Databases/ArchBactBLAST_DB_v2/CombinedBactArch.faa > /Applications/Databases/ArchBactBLAST_DB_v2_xAsgard/CombinedBactArchxAsgard.faa

	##Old
	for i in $(\ls /Users/jmattic1/Documents/Laura_proteomes/*.faa | grep "faa")
	do

		cat $i | sed 's/>/>lcl\|/g' >> /Applications/Databases/EukArchBact_xAsgard/EukArchBact_xAsgard.faa
		cat $i | sed 's/>/>lcl\|/g' >> /Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa

	done

	##New
	for i in $(\ls /Users/jmattic1/Documents/Laura_proteomes/*.faa | grep "faa")
	do
		name=$(echo $i | awk -F '/' '{print $NF}' | sed 's/.faa//g')
		cat $i | sed "s/>.*\.\./>lcl\|$name\_/g" >> /Applications/Databases/EukArchBact_xAsgard/EukArchBact_xAsgard.faa
		cat $i | sed "s/>.*\.\./>lcl\|$name\_/g" >> /Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa
		cat $i | sed "s/>.*\.\./>lcl\|$name\_/g" >>  /Applications/Databases/ArchBactBLAST_DB_v2_xAsgard/CombinedBactArchxAsgard.faa
	done
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/EukArchBact_xAsgard/EukArchBact_xAsgard.faa -dbtype 'prot' -parse_seqids
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/EukArchBactViral_xAsgard/EukArchBactViral_xAsgard.faa -dbtype 'prot' -parse_seqids
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/ArchBactBLAST_DB_v2_xAsgard/CombinedBactArchxAsgard.faa -dbtype 'prot' -parse_seqids


fi

###Subset the RBP alignment
if [[ 1 -eq 2 ]]
then
	outdir="/Users/jmattic1/Documents/EukArchBactTree_automatedxAsgard/SubsetTrees/"
	rm -rf $outdir
	mkdir $outdir
	treefile="/Users/jmattic1/Documents/EukArchBactTree_automatedxAsgard/ConcatTrimmed/ConcatenatedRibosomeProteins.aligned.trimmed.folded.faa.rate"
	faafile="/Users/jmattic1/Documents/EukArchBactTree_automatedxAsgard/ConcatTrimmed/ConcatenatedRibosomeProteins.aligned.trimmed.faa"

	cat $treefile | grep -v "#" | tail -n+2 > $outdir/RateFile.tsv
	sitenumber=$(cat $outdir/RateFile.tsv | wc -l)
	echo -e "10\t"$sitenumber > $outdir/RateSites.tsv
	for i in $(seq 9 -1 1)
	do
		cat $outdir/RateFile.tsv | awk -v limit="$i" '$3 <= limit {print $0}' > $outdir/RateFile.$i.tsv
		sitenumber=$(cat $outdir/RateFile.$i.tsv | wc -l)
		echo -e $i"\t"$sitenumber >> $outdir/RateSites.tsv
		newfaafile=$(echo $faafile | awk -F '/' '{print $NF}' | sed -e "s/.faa/."$i".faa/g")
		positions=$(cat $outdir/RateFile.$i.tsv | awk '{print $1}' | tr '\n' ',' | sed 's/,$//g')
		cat $faafile | awk -v positions="$positions" '{if($0 ~ /^>/) print $0; else {split($0, seqarray,""); split(positions, posarray,","); for(count = 1; count <= length(posarray); ++count){printf seqarray[posarray[count]];}; print ""}}' > $outdir/$newfaafile

		/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir/$newfaafile --alrt 1000 -B 1000

	done
fi


###Banoos candidates
if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/Banoo_Candidate_UPs/
	mkdir /Users/jmattic1/Documents/Banoo_Candidate_UPs/
	conda activate samtools
	cat /Users/jmattic1/Documents/banoosqueries.faa | sed 's/\|/_/g' >  /Users/jmattic1/Documents/banoosqueries.nopipes.faa
	for i in $(cat /Users/jmattic1/Documents/banoosqueries.nopipes.faa | grep ">" | sed 's/>//g' | awk '{print $1}')
	do
		outfile=$(echo $i".faa")
		samtools faidx /Users/jmattic1/Documents/banoosqueries.nopipes.faa $i > /Users/jmattic1/Documents/Banoo_Candidate_UPs/$outfile
	done
fi


###Euk Candidates in general
if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/EukTree_UPs/
	mkdir /Users/jmattic1/Documents/EukTree_UPs/
	conda activate /Applications/miniconda3/envs/samtools
	#Already done riboproteins, no need to do them here
	for i in $(\ls /Users/jmattic1/Documents/EukaryoteTreeProteins/ | grep ".fa" | grep -v "rps" | grep -v "rpl" | grep -v ".fai")
	do
		longestseq=$(perl /Applications/residues.pl /Users/jmattic1/Documents/EukaryoteTreeProteins/$i | sort -k2,2nr | awk '{print $1}' | head -1)
		outfile=$(echo $i | sed 's/fa/faa/g')
		samtools faidx /Users/jmattic1/Documents/EukaryoteTreeProteins/$i $longestseq | sed 's/@/_/g' | sed 's/\|/_/g' | sed 's/\//_/g' > /Users/jmattic1/Documents/EukTree_UPs/$outfile
	done
fi


####Combine the Euk candidates with the RPs

###Updated for v4 with all asgard
if [[ 1 -eq 1 ]]
then
	#outdir="/Users/jmattic1/Documents/EukandRiboprotein_ConcatTree/"
	#EukProteinDir="/Users/jmattic1/Documents/EukTree_UPs_blast_v3.1/"
	outdir="/Users/jmattic1/Documents/EukandRiboprotein_ConcatTree_60perc_v4_AllAsgard/"
	EukProteinDir="/Users/jmattic1/Documents/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard_Proteome_Output/"
	RiboProteinDir="/Users/jmattic1/Documents/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard_RiboProteins_blast_v4_60perc/"
	rm -rf $outdir
	mkdir $outdir
	specieslist=$(echo "/Users/jmattic1/Documents/IncludedTaxa.txt")
	speciesnumber=$(cat $specieslist | wc -l)
	###Make this list-able
	rm -f $outdir/IncludedProteins.txt
	while read line
	do
		grep -l $line $RiboProteinDir/*/*trimal.faa | awk -F '/' '{print $(NF-1)}' >> $outdir/IncludedProteins.txt

	done < $specieslist
	cat $outdir/IncludedProteins.txt | sort | uniq -c | awk -v speciesnumber="$speciesnumber" '$1 == speciesnumber {print $2}' > $outdir/FullBactProteins.txt 
	
	rm -f $outdir/IncludedProteins.txt
	while read line
	do
		grep -l $line $EukProteinDir/*/*trimal.faa | awk -F '/' '{print $(NF-1)}' >> $outdir/IncludedProteins.txt

	done < $specieslist
	cat $outdir/IncludedProteins.txt | sort | uniq -c | awk -v speciesnumber="$speciesnumber" '$1 == speciesnumber {print $2}' > $outdir/EukFullBactProteins.txt 

	rm -rf  $outdir/HodContainingRPs/
	mkdir $outdir/HodContainingRPs/
	for i in $(cat $outdir/FullBactProteins.txt | grep "RP")
	do
		cp -r $RiboProteinDir/$i $outdir/HodContainingRPs/
	done

	for i in $(cat $outdir/EukFullBactProteins.txt | grep -E "e|running")
	do
		cp -r $EukProteinDir/$i $outdir/HodContainingRPs/
	done

	genenumber=$(ls $outdir/HodContainingRPs/ | wc -l)
	cat $outdir/HodContainingRPs/*/*.trimal.faa | grep ">" | sed 's/>//g' | awk -F "_" '{print $1"_"$2}' | sort | uniq -c | awk -v genenumber="$genenumber" '$1 == genenumber {print $2}' > $outdir/IncludedSpecies.txt

	cat $outdir/FullBactProteins.txt $outdir/EukFullBactProteins.txt > $outdir/AllProteins.txt

	while read line
	do
		for i in $(cat $outdir/AllProteins.txt)
		do
			fasta=$(\ls $outdir/HodContainingRPs/$i/*.align.trimal.faa)

			extractiongene=$(cat $fasta | grep $line"_" | sed 's/>//g')
			samtools faidx $fasta $extractiongene >> $outdir/Unconcat.faa
		done
		echo ">"$line > $outdir/IndvConcat.faa
		cat $outdir/Unconcat.faa | grep -v ">" | tr -d "\n" >> $outdir/IndvConcat.faa
		rm -f $outdir/Unconcat.faa
		fold -w 80 $outdir/IndvConcat.faa >> $outdir/RibosomalProteins.concatenated.uncondensed.faa
	done < $outdir/IncludedSpecies.txt
	cat $outdir/RibosomalProteins.concatenated.uncondensed.faa | sed 's/>/\n>/g' | sed '/^$/d' > $outdir/RibosomalProteins.concatenated.final.faa
	mv $outdir/RibosomalProteins.concatenated.final.faa $outdir/EukRibosomalProteins.concatenated.final.faa
	/Applications/miniconda3/bin/iqtree -T 20 -s $outdir/EukRibosomalProteins.concatenated.final.faa --alrt 1000 -B 1000 -wsr
	
	/Applications/aa_recoder.py -i $outdir/EukRibosomalProteins.concatenated.final.faa -re SR4 -o $outdir/EukRibosomalProteins.concatenated.final.SR4
	/Applications/miniconda3/bin/iqtree -T 10 -s $outdir/EukRibosomalProteins.concatenated.final.SR4.fas --alrt 1000 -B 1000 -wsr

fi


####Combine the Euk candidates with the RPs -- UNTRIMMED AS THE PAPER INDICATES IS BETTER
if [[ 1 -eq 2 ]]
then
	#outdir="/Users/jmattic1/Documents/EukandRiboprotein_ConcatTree/"
	#EukProteinDir="/Users/jmattic1/Documents/EukTree_UPs_blast_v3.1/"
	outdir="/Users/jmattic1/Documents/EukandRiboprotein_ConcatTree_60perc_untrimmed/"
	EukProteinDir="/Users/jmattic1/Documents/EukTree_UPs_blast_v3.1_60perc/"
	rm -rf $outdir
	mkdir $outdir
	specieslist=$(echo "/Users/jmattic1/Documents/IncludedTaxa.txt")
	speciesnumber=$(cat $specieslist | wc -l)
	###Make this list-able
	rm -f $outdir/IncludedProteins.txt
	while read line
	do
		grep -l $line /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_reducedset_v3.1/*/*align.faa | awk -F '/' '{print $(NF-1)}' >> $outdir/IncludedProteins.txt

	done < $specieslist
	cat $outdir/IncludedProteins.txt | sort | uniq -c | awk -v speciesnumber="$speciesnumber" '$1 == speciesnumber {print $2}' > $outdir/FullBactProteins.txt 
	
	rm -f $outdir/IncludedProteins.txt
	while read line
	do
		grep -l $line $EukProteinDir/*/*align.faa | awk -F '/' '{print $(NF-1)}' >> $outdir/IncludedProteins.txt

	done < $specieslist
	cat $outdir/IncludedProteins.txt | sort | uniq -c | awk -v speciesnumber="$speciesnumber" '$1 == speciesnumber {print $2}' > $outdir/EukFullBactProteins.txt 

	rm -rf  $outdir/HodContainingRPs/
	mkdir $outdir/HodContainingRPs/
	for i in $(cat $outdir/FullBactProteins.txt | grep "RP")
	do
		cp -r /Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_reducedset_v3.1/$i $outdir/HodContainingRPs/
	done

	for i in $(cat $outdir/EukFullBactProteins.txt | grep "e")
	do
		cp -r $EukProteinDir/$i $outdir/HodContainingRPs/
	done

	genenumber=$(ls $outdir/HodContainingRPs/ | wc -l)
	cat $outdir/HodContainingRPs/*/*.align.faa | grep ">" | sed 's/>//g' | awk -F "_" '{print $1"_"$2}' | sort | uniq -c | awk -v genenumber="$genenumber" '$1 == genenumber {print $2}' > $outdir/IncludedSpecies.txt

	cat $outdir/FullBactProteins.txt $outdir/EukFullBactProteins.txt > $outdir/AllProteins.txt

	while read line
	do
		for i in $(cat $outdir/AllProteins.txt)
		do
			fasta=$(\ls $outdir/HodContainingRPs/$i/*.align.faa)

			extractiongene=$(cat $fasta | grep $line"_" | sed 's/>//g')
			samtools faidx $fasta $extractiongene >> $outdir/Unconcat.faa
		done
		echo ">"$line > $outdir/IndvConcat.faa
		cat $outdir/Unconcat.faa | grep -v ">" | tr -d "\n" >> $outdir/IndvConcat.faa
		rm -f $outdir/Unconcat.faa
		fold -w 80 $outdir/IndvConcat.faa >> $outdir/RibosomalProteins.concatenated.uncondensed.faa
	done < $outdir/IncludedSpecies.txt
	cat $outdir/RibosomalProteins.concatenated.uncondensed.faa | sed 's/>/\n>/g' | sed '/^$/d' > $outdir/RibosomalProteins.concatenated.final.faa
	mv $outdir/RibosomalProteins.concatenated.final.faa $outdir/EukRibosomalProteins.concatenated.final.faa
	/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir/EukRibosomalProteins.concatenated.final.faa --alrt 1000 -B 1000 -wsr
	
	/Applications/aa_recoder.py -i $outdir/EukRibosomalProteins.concatenated.final.faa -re SR4 -o $outdir/EukRibosomalProteins.concatenated.final.SR4
	/Applications/miniconda3/bin/iqtree -T AUTO -s $outdir/EukRibosomalProteins.concatenated.final.SR4.fas --alrt 1000 -B 1000 -wsr



fi

###Subsetting SR4
if [[ 1 -eq 2 ]]
then
	outdir="/Users/jmattic1/Documents/EukandRiboprotein_ConcatTree/SubsetTree_SR4/"
	rm -rf $outdir
	mkdir $outdir
	treefile="/Users/jmattic1/Documents/EukandRiboprotein_ConcatTree/EukRibosomalProteins.concatenated.final.SR4.fas.rate"
	staticfaafile="/Users/jmattic1/Documents/EukandRiboprotein_ConcatTree/EukRibosomalProteins.concatenated.final.SR4.fas"
	faafile=$(echo $staticfaafile | awk -F '/' -v outdir="$outdir" '{print outdir"/"$NF}')
	cat $staticfaafile | awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' > $faafile
	cat $treefile | grep -v "#" | tail -n+2 > $outdir/RateFile.tsv
	sitenumber=$(cat $outdir/RateFile.tsv | wc -l)
	echo -e "10\t"$sitenumber > $outdir/RateSites.tsv
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
fi

###Converting to Phylobayes' stupid phylip format jfc
if [[ 1 -eq 2 ]]
then
	conda activate /Applications/miniconda3/envs/samtools
	inputfasta=$(\ls /Users/jmattic1/Documents/EukandRiboprotein_ConcatTree/EukRibosomalProteins.concatenated.final.faa)
	echo "220 25546" > /Users/jmattic1/Documents/EukRibosomalProteins.concatenated.final.unfolded.phy
	cat $inputfasta | awk '{if($0 ~ /^>/) print $0"\t"; else print $0}' | tr -d '\n' | tr '>' '\n' >> /Users/jmattic1/Documents/EukRibosomalProteins.concatenated.final.unfolded.phy
	sed '/^$/d' /Users/jmattic1/Documents/EukRibosomalProteins.concatenated.final.unfolded.phy > /Users/jmattic1/Documents/EukRibosomalProteins.concatenated.final.unfolded.lines.phy
fi
