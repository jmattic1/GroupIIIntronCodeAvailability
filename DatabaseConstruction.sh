if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/BactRefSeq
	mkdir /Users/jmattic1/Documents/BactRefSeq
	rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria /Users/jmattic1/Documents/BactRefSeq
	rclone copy ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria /Users/jmattic1/Documents/BactRefSeq
	rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/*/latest_assembly_versions/*/ /Users/jmattic1/Documents/BactRefSeq

	rm -rf /Users/jmattic1/Documents/BactGenbank
	mkdir /Users/jmattic1/Documents/BactGenbank
	rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria /Users/jmattic1/Documents/BactGenbank
	
	rm -rf /Users/jmattic1/Documents/BactRefSeq2
	mkdir /Users/jmattic1/Documents/BactRefSeq2

	rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/*.faa.gz /Users/jmattic1/Documents/BactRefSeq2/
fi

if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/ArchaeaRefSeq
	mkdir /Users/jmattic1/Documents/ArchaeaRefSeq
	rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea /Users/jmattic1/Documents/ArchaeaRefSeq

	rm -rf /Users/jmattic1/Documents/ArchaeaGenbank
	rm /Users/jmattic1/Documents/ArchaeaGenbank.zip
	mkdir /Users/jmattic1/Documents/ArchaeaGenbank
	rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea /Users/jmattic1/Documents/ArchaeaGenbank
	

fi

###Get viral sequences from refseq/genbank
if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/ViralRefSeq
	mkdir /Users/jmattic1/Documents/ViralRefSeq
	rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral /Users/jmattic1/Documents/ViralRefSeq

	rm -rf /Users/jmattic1/Documents/ViralGenbank
	mkdir /Users/jmattic1/Documents/ViralGenbank
	rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral /Users/jmattic1/Documents/ViralGenbank
	

fi

#####Build the Archaeal BlastDB
if [[ 1 -eq 2 ]]
then
	###Get all of the usable FASTAs (then remove the archaeal genbank directory -- too large)
	if [[ 1 -eq 2 ]]
	then
		rm -rf /Users/jmattic1/Documents/ArchaealBLAST_FASTA/
		mkdir /Users/jmattic1/Documents/ArchaealBLAST_FASTA/
		cd /Users/jmattic1/Documents/ArchaealBLAST_FASTA/

		for i in $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea)
		do
			if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$i/ | grep "latest_assembly_versions" | wc -l) -gt 0 ]]
			then
				if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$i/latest_assembly_versions/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
				then
					proteinfile=$(ls -lt /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$i/latest_assembly_versions/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
					cp $proteinfile /Users/jmattic1/Documents/ArchaealBLAST_FASTA/$i.faa.gz
				fi
			elif [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$i/ | grep "all_assembly_versions" | wc -l) -gt 0 ]]
			then
				if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$i/all_assembly_versions/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
				then
					proteinfile=$(ls -lt /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$i/all_assembly_versions/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
					cp $proteinfile /Users/jmattic1/Documents/ArchaealBLAST_FASTA/$i.faa.gz
				elif [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$i/all_assembly_versions/ | grep "suppressed" | wc -l) -gt 0 ]]
				then
					if [[ $(\ls /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$i/all_assembly_versions/suppressed/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
					then
						proteinfile=$(ls -lt /Users/jmattic1/Documents/ArchaeaGenbank/archaea/$i/all_assembly_versions/suppressed/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
						cp $proteinfile /Users/jmattic1/Documents/ArchaealBLAST_FASTA/$i.faa.gz
					fi
				fi
			fi

		done

		gzip -d /Users/jmattic1/Documents/ArchaealBLAST_FASTA/*
	fi
	###Obselete
	if [[ 1 -eq 2 ]]
	then
		##Build the actual BLAST DB
		rm -rf /Users/jmattic1/Documents/ArchaealBLAST_DB/
		mkdir /Users/jmattic1/Documents/ArchaealBLAST_DB/
		cd /Users/jmattic1/Documents/ArchaealBLAST_DB/
		for i in $(\ls /Users/jmattic1/Documents/ArchaealBLAST_FASTA/ | grep ".faa")
		do
			species=$(echo $i | sed 's/\.faa//g' | sed 's/_archaeon//g' | sed 's/Candidatus_//g' | sed 's/miscellaneous_//g' | sed 's/BIG_FIL_POST_REV_//g' | sed 's/big_fil_rev_//g' | sed 's/um_filter//g' | sed 's/Marine_Group_//g' | sed 's/Methano.*_//g' | sed 's/land_//g' | sed 's/et_al.*//g' | sed 's/BACL13_MAG-//g' | sed 's/CG//g' | sed 's/Nitrosopumilus_sp\._ex/Nitrosopumilus/g')
			cat /Users/jmattic1/Documents/ArchaealBLAST_FASTA/$i | awk '{print $1}' | sed "s/>/>lcl\|$species/g" >> /Users/jmattic1/Documents/ArchaealBLAST_DB/Archaea.faa
		done

		/Users/jmattic1/miniconda3/bin/makeblastdb -in /Users/jmattic1/Documents/ArchaealBLAST_DB/Archaea.faa -dbtype 'prot' -parse_seqids
	fi
fi

#####Get Viral FASTAs
if [[ 1 -eq 2 ]]
then
	###Get all of the usable FASTAs (then remove the archaeal genbank directory -- too large)
	if [[ 1 -eq 1 ]]
	then
		rm -rf /Users/jmattic1/Documents/ViralBLAST_FASTA/
		mkdir /Users/jmattic1/Documents/ViralBLAST_FASTA/
		cd /Users/jmattic1/Documents/ViralBLAST_FASTA/

		for i in $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral)
		do
			if [[ $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral/$i/ | grep "latest_assembly_versions" | wc -l) -gt 0 ]]
			then
				if [[ $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral/$i/latest_assembly_versions/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
				then
					proteinfile=$(ls -lt /Users/jmattic1/Documents/ViralRefSeq/viral/$i/latest_assembly_versions/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
					cp $proteinfile /Users/jmattic1/Documents/ViralBLAST_FASTA/$i.faa.gz
					gzip -d /Users/jmattic1/Documents/ViralBLAST_FASTA/$i.faa.gz
				fi
			elif [[ $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral/$i/ | grep "all_assembly_versions" | wc -l) -gt 0 ]]
			then
				if [[ $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral/$i/all_assembly_versions/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
				then
					proteinfile=$(ls -lt /Users/jmattic1/Documents/ViralRefSeq/viral/$i/all_assembly_versions/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
					cp $proteinfile /Users/jmattic1/Documents/ViralBLAST_FASTA/$i.faa.gz
					gzip -d /Users/jmattic1/Documents/ViralBLAST_FASTA/$i.faa.gz
				elif [[ $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral/$i/all_assembly_versions/ | grep "suppressed" | wc -l) -gt 0 ]]
				then
					if [[ $(\ls /Users/jmattic1/Documents/ViralRefSeq/viral/$i/all_assembly_versions/suppressed/*/ | grep "protein.faa.gz" | wc -l) -gt 0 ]]
					then
						proteinfile=$(ls -lt /Users/jmattic1/Documents/ViralRefSeq/viral/$i/all_assembly_versions/suppressed/*/*protein.faa.gz | sort -k5,5nr | awk '{print $NF}' | head -1)
						cp $proteinfile /Users/jmattic1/Documents/ViralBLAST_FASTA/$i.faa.gz
						gzip -d /Users/jmattic1/Documents/ViralBLAST_FASTA/$i.faa.gz
					fi
				fi
			fi

		done

		#gzip -d /Users/jmattic1/Documents/ViralBLAST_FASTA/*
	fi

	##Build the actual BLAST DB

fi
###Euk Genbank transcriptome download
if [[ 1 -eq 2 ]]
then
	#/Users/jmattic1/Documents/EukGenbankTranscriptomes/
	cp -r /Users/jmattic1/Documents/EukGenbankTranscriptomes/ /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/
	cd /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/
	for f in *\ *; do mv "$f" "${f// /_}"; done
	for i in $(\ls /Users/jmattic1/Documents/EukGenbankTranscriptomesMod/*.zip | grep "zip")
	do
		newdir=$(echo $i | sed 's/.zip//g')
		mkdir $newdir
		cd $newdir
		cp $i $newdir
		unzip ./*.zip
	done

fi


####Archaea + Bacteria
if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/ArchBactBLAST_DB/
	mkdir /Users/jmattic1/Documents/ArchBactBLAST_DB/
	
	cd /Users/jmattic1/Documents/EukProtProteinDatabase/

	cat /Users/jmattic1/Documents/EukProteomeRefFile.v3.loop.filtered.csv | grep "Bacteria" > /Users/jmattic1/Documents/EukProteomeRefFile.v3.loop.bacterial.filtered.csv

	while read line
	do
		speciesunedited=$(echo $line | awk -F ',' '{print $1}' | sed 's/ /_/g')
		if [[ $(echo $speciesunedited | grep "_sp" | wc -l) -eq 0 && $(echo $speciesunedited | grep "_cf" | wc -l) -eq 0 && $(echo $speciesunedited | awk -F "_" '{print NF}') -eq 2 ]]
		then
			species=$(echo $speciesunedited | awk -F '_' '{print $1"_"substr($2,1,1)}')
		else
			species=$(echo $speciesunedited)
		fi
		echo $speciesunedited
		echo $species
		faafile=$(echo $line | awk -F ',' '{print $3".faa"}')
		if [[ $(echo $faafile | grep "GCA_" | wc -l) -gt 0 || $(echo $faafile | grep "GCF_" | wc -l) -gt 0 ]]
		then
			cat $faafile | sed -e "s/>/>"$species"_/g" | sed 's/ /_/g' | sed 's/\_\[.*//g' >> /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa
		elif [[ $(echo $faafile | grep "_Gene" | wc -l) -gt 0 ]]
		then
			cat $faafile | sed -e "s/>jgi/>"$species"_/g" | sed 's/\|/_/g' >> /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa
		elif [[ $(echo $faafile | grep "EP" | wc -l) -gt 0 ]]
		then
			cat $faafile | sed 's/ .*//g' | awk -F "_" '{if($0 ~ /^>/)  print $1"_"$NF; else print $0;}' | sed 's/ .*//g' | sed -e "s/"$species"_//g" | sed -e "s/>/>"$species"_/g" >> /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa
		elif [[ $(echo $faafile | grep "Bigna" | wc -l) -gt 0 ]]
		then
			cat $faafile | sed -e "s/>/>"$species"_/g" | sed 's/ .*gene:/_/g' | sed 's/ .*//g' >> /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa
		elif [[ $(echo $faafile | grep "_prot_cdhit95" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F "::" '{if($0 ~ /^>/)  print ">"$2; else print $0;}' | sed -e "s/>/>"$species"_/g" >> /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa
		elif [[ $(echo $faafile | grep "Nibbleromonas" | wc -l) -gt 0 || $(echo $faafile | grep "Ubysesseya" | wc -l) -gt 0 || $(echo $faafile | grep "Nebulamonas" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F "::" '{if($0 ~ /^>/)  print ">"$2"_"$3; else print $0;}' | sed 's/NEBULAMONAS_//g' | sed 's/>.*DN.*_c/>c/g' | sed -e "s/>/>"$species"_/g" >> /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa
		elif [[ $(echo $faafile | grep "Nephromyces" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F " " '{if($0 ~ /^>/)  print ">"$1; else print $0;}' | sed 's/>.*DN.*_c/>c/g' | sed 's/\|/_/g' | sed -e "s/>/>"$species"_/g" >> /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa
		elif [[ $(echo $faafile | grep "Malawimonas_californiana" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F " " '{if($0 ~ /^>/)  print ">"$1; else print $0;}' | awk -F "_" '{if($0 ~ /^>/)  print ">"$3; else print $0;}' | sed -e "s/>/>"$species"_/g" >> /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa
		elif [[ $(echo $faafile | grep "Selenidium_serpulae" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F "@" '{if($0 ~ /^>/)  print $1"_"$2; else print $0;}' | sed -e "s/Selendium_serpulae/Selenidium_s/g" >> /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa
		elif [[ $(echo $faafile | grep "Lankesteria_metandrocarpae" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F "@" '{if($0 ~ /^>/)  print $1"_"$2; else print $0;}' | sed -e "s/Lankesteria_metandrocarpae/Lankesteria_m/g" >> /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa
		else
			cat $faafile | awk -F "@" '{if($0 ~ /^>/)  print $1"_"$2; else print $0;}' | sed -e "s/>/>"$species"_/g" >> /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa
		fi

	done < /Users/jmattic1/Documents/EukProteomeRefFile.v3.loop.bacterial.filtered.csv

	cat /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa | sed '/^$/d' | sed 's/SV\*>Cafeteria_burkhardae_EP00749/SV\*\n>C_burkhardae_EP00749/g' | sed 's/Ochro1393_1_4_//g' | awk '{print $1}' | sed 's/,.*//g' | sed 's/_TRINITY//g' | sed 's/\|/_/g' | sed 's/>/>lcl\|/g' | awk '{if($0 ~ /^>/)  print substr($0,1,50); else print $0;}' > /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.mod.faa
	mv /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.mod.faa /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa

	cat /Users/jmattic1/Documents/ArchBactBLAST_DB/Bact.protein.db.faa > /Users/jmattic1/Documents/ArchBactBLAST_DB/CombinedBactArch.faa

	rm -rf /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/
	cp -r /Users/jmattic1/Documents/ArchaealBLAST_FASTA/ /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/

	paste -d "\t" <(\ls /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/ | grep ".faa") <(\ls /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/ | grep ".faa"  | sed 's/\.faa//g'| sed 's/Candidatus_/C/g' | sed 's/candidate_division/c/g' | sed 's/_archaeon_//g' | sed 's/_big_fil_rev_//g' | sed 's/archaeota//g' | sed 's/miscellaneous_//g' | awk '{print substr($0,1,27)}' | awk -F "_" '{if($1 == $NF){print $1"_arch"} else {print $1"_"$NF}}') | awk '{print $1"\t"$2}' > /Users/jmattic1/Documents/ArchaeaToBeFixed.tsv

	rm -f /Users/jmattic1/Documents/FinishedFastas.txt
	touch /Users/jmattic1/Documents/FinishedFastas.txt
	while read line
	do
		fasta=$(echo $line | awk '{print $1}')
		presence=$(cat /Users/jmattic1/Documents/FinishedFastas.txt | awk -v fasta="$fasta" '$1 == fasta {print $0}')
		if [[ $presence == "" ]]
		then
			genus=$(echo $line | awk '{print $2}')
			cat /Users/jmattic1/Documents/ArchaeaToBeFixed.tsv | awk -v genus="$genus" '$2 == genus {print $0}' > /Users/jmattic1/Documents/subArchaeaToBeFixed.tsv
			duplicates=$(cat /Users/jmattic1/Documents/subArchaeaToBeFixed.tsv | wc -l)
			if [[ $duplicates -gt 1 ]]
			then
				count=1
				while read change
				do
					newgenus=$(echo $change | awk '{print $2}')
					newfasta=$(echo $change | awk '{print $1}')
					echo -e $newfasta"\t"$newgenus"-L"$count".faa" >> /Users/jmattic1/Documents/FinishedFastas.txt

					mv /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/$newfasta /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/$newgenus"-L"$count".faa"
					count=$(( $count + 1 ))
				done < /Users/jmattic1/Documents/subArchaeaToBeFixed.tsv
			else
				echo -e $fasta"\t"$genus"-L1.faa" >> /Users/jmattic1/Documents/FinishedFastas.txt
				mv /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/$fasta /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/$genus"-L1.faa"

			fi
		fi
	done < /Users/jmattic1/Documents/ArchaeaToBeFixed.tsv

	for i in $(\ls /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/ | grep ".faa")
	do
		species=$(echo $i | sed 's/\.faa//g')
		#if [[ $(cat /Users/jmattic1/Documents/ArchBactBLAST_DB/CombinedBactArch.faa | grep $species | wc -l) -eq 0 ]]
		#then
			cat /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/$i | awk '{print $1}' | sed "s/>/>lcl\|$species\_/g" | awk '{if($0 ~ /^>/)  print substr($0,1,50); else print $0;}' >> /Users/jmattic1/Documents/ArchBactBLAST_DB/CombinedBactArch.faa
		#fi
	done
	/Applications/miniconda3/bin/makeblastdb -in /Users/jmattic1/Documents/ArchBactBLAST_DB/CombinedBactArch.faa -dbtype 'prot' -parse_seqids

fi



####Eukaryota

#########FIX 
if [[ 1 -eq 2 ]]
then
	rm /Users/jmattic1/Documents/EukSpeciesKey.tsv
	rm -rf /Users/jmattic1/Documents/EukBLASTDB/
	mkdir /Users/jmattic1/Documents/EukBLASTDB/
	cd /Users/jmattic1/Documents/EukProtProteinDatabase/
	cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.csv | sed 's/.faa//g' | grep -v "ProteomeRef" | grep -v "No genome" | grep -v "No proteome" | grep -v "^,," | grep -v "PRJNA" | grep -v "Genus-Species" > /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.filtered.csv
	cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.filtered.csv | awk -F ',' '$4 == "Eukaryota" {print $0}' | grep -v "NCBI_TaxID" > /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.euk.filtered.csv
	cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.filtered.csv | awk -F ',' '$4 == "Bacteria" {print $0}' | grep -v "NCBI_TaxID" > /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.bact.filtered.csv

	cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.filtered.csv | awk -F ',' '$4 == "Eukaryota" {print $0}' | grep -v "NCBI_TaxID" > /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.eukbact.filtered.csv
	cat /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.filtered.csv | awk -F ',' '$4 == "Bacteria" {print $0}' | grep -v "NCBI_TaxID" >> /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.eukbact.filtered.csv
	##EukProk Protein Database
	while read line
	do
		truespeciesunedited=$(echo $line | awk -F ',' '{print $1}')
		speciesunedited=$(echo $line | awk -F ',' '{print $1}' | sed 's/ /_/g' | sed 's/_cf//g' | sed 's/\.//g')
		if [[ $(echo $speciesunedited | grep -v "specki" | grep "_sp" | wc -l) -eq 0 && $(echo $speciesunedited | grep "_cf" | wc -l) -eq 0 && $(echo $speciesunedited | awk -F "_" '{print NF}') -eq 2 ]]
		then
			species=$(echo $speciesunedited | awk -F '_' '{print $1"_"substr($2,1,1)}')
		elif [[ $(echo $speciesunedited | awk -F "_" '{print NF}') -lt 2 ]]
		then
			species=$(echo $speciesunedited | awk -F "_" '{print $1"_sp"}')
		else
			species=$(echo $speciesunedited | awk -F "_" '{print $1"_"$2}')
		fi

		echo $truespeciesunedited
		echo $species
		echo -e $truespeciesunedited"\t"$species >> /Users/jmattic1/Documents/EukSpeciesKey.tsv
		faafile=$(echo $line | awk -F ',' '{print $3".faa"}')
		if [[ $(echo $faafile | grep "GCA_" | wc -l) -gt 0 || $(echo $faafile | grep "GCF_" | wc -l) -gt 0 ]]
		then
			cat $faafile | sed 's/Mantamonas_s//g' | sed -e "s/>/>"$species"_/g" | sed 's/ /_/g' | sed 's/\_\[.*//g' >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		elif [[ $(echo $faafile | grep "_Gene" | wc -l) -gt 0 ]]
		then
			cat $faafile | sed -e "s/>jgi/>"$species"_/g" | sed 's/\|/_/g' >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		elif [[ $(echo $faafile | grep "EP" | wc -l) -gt 0 ]]
		then
			cat $faafile | sed 's/ .*//g' | awk -F "_" '{if($0 ~ /^>/)  print $1"_"$NF; else print $0;}' | sed 's/ .*//g' | sed -e "s/"$species"_//g" | sed -e "s/>/>"$species"_/g" >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		elif [[ $(echo $faafile | grep "Bigna" | wc -l) -gt 0 ]]
		then
			cat $faafile | sed -e "s/>/>"$species"_/g" | sed 's/ .*gene:/_/g' | sed 's/ .*//g' >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		elif [[ $(echo $faafile | grep "_prot_cdhit95" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F "::" '{if($0 ~ /^>/)  print ">"$2; else print $0;}' | sed -e "s/>/>"$species"_/g" >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		elif [[ $(echo $faafile | grep "Nibbleromonas" | wc -l) -gt 0 || $(echo $faafile | grep "Ubysesseya" | wc -l) -gt 0 || $(echo $faafile | grep "Nebulamonas" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F "::" '{if($0 ~ /^>/)  print ">"$2"_"$3; else print $0;}' | sed 's/NEBULAMONAS_//g' | sed 's/>.*DN.*_c/>c/g' |  sed -e "s/>/>"$species"_/g" >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		elif [[ $(echo $faafile | grep "Nephromyces" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F " " '{if($0 ~ /^>/)  print ">"$1; else print $0;}' | sed 's/>.*DN.*_c/>c/g' | sed 's/\|/_/g' | sed -e "s/>/>"$species"_/g" >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		elif [[ $(echo $faafile | grep "Malawimonas_californiana" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F " " '{if($0 ~ /^>/)  print ">"$1; else print $0;}' | awk -F "_" '{if($0 ~ /^>/)  print ">"$3; else print $0;}' | sed -e "s/>/>"$species"_/g" >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		elif [[ $(echo $faafile | grep "Selenidium_serpulae" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F "@" '{if($0 ~ /^>/)  print $1"_"$2; else print $0;}' | sed -e "s/Selendium_serpulae/Selenidium_s/g" >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		elif [[ $(echo $faafile | grep "Lankesteria_metandrocarpae" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F "@" '{if($0 ~ /^>/)  print $1"_"$2; else print $0;}' | sed -e "s/Lankesteria_metandrocarpae/Lankesteria_m/g" >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		elif [[ $(echo $faafile | grep "Filipodium_phascolosomae" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F "@" '{if($0 ~ /^>/)  print $1"_"$2; else print $0;}' | sed -e "s/Filipodium_phascolosomae/Filipodium_p/g" >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		elif [[ $(echo $faafile | grep "Platyproteum_sp" | wc -l) -gt 0 ]]
		then
			cat $faafile | awk -F "@" '{if($0 ~ /^>/)  print $1"_"$2; else print $0;}' | sed -e "s/Platyproteum_1/"$species"/g" | sed -e "s/Platyproteum_2/"$species"/g" >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		else
			cat $faafile | awk -F "@" '{if($0 ~ /^>/)  print $1"_"$2; else print $0;}' | sed "s/$speciesunedited//g" | sed -e "s/>/>"$species"_/g" >> /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
		fi

	done < /Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.euk.filtered.csv

	cat /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa | sed '/^$/d' | sed 's/SV\*>Cafeteria_burkhardae_EP00749/SV\*\n>C_burkhardae_EP00749/g' | sed 's/Ochro1393_1_4_//g' | awk '{print $1}' | sed 's/,.*//g' | sed 's/_TRINITY//g' | sed 's/\|/_/g' | sed 's/>/>lcl\|/g' | awk '{if($0 ~ /^>/)  print substr($0,1,50); else print $0;}' > /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.mod.faa
	mv /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.mod.faa /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa
	
	cat /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa | grep ">" | awk -F "_" '{print $1"_"$2}' | sort | uniq -c
	/Applications/miniconda3/bin/makeblastdb -in /Users/jmattic1/Documents/EukBLASTDB/Euk.protein.db.faa -dbtype 'prot' -parse_seqids
	

	rm -rf /Applications/Databases/EukBLASTDB/
	mv /Users/jmattic1/Documents/EukBLASTDB/ /Applications/Databases/EukBLASTDB/
fi


####Viruses
if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/Virus_DB/
	mkdir /Users/jmattic1/Documents/Virus_DB/

	rm -rf /Users/jmattic1/Documents/ViralBLAST_FASTA_modded/
	cp -r /Users/jmattic1/Documents/ViralBLAST_FASTA/ /Users/jmattic1/Documents/ViralBLAST_FASTA_modded/

	paste -d "\t" <(\ls /Users/jmattic1/Documents/ViralBLAST_FASTA_modded/ | grep ".faa") <(\ls /Users/jmattic1/Documents/ViralBLAST_FASTA_modded/ | grep ".faa"  | sed 's/\.faa//g'| sed 's/-like/C/g' | sed 's/_virus/c/g' | sed 's/associated_/c/g' | sed 's/strain//g' | awk '{print substr($0,1,20)}' | awk -F "_" '{if($1 == $NF){print $1"_virus"} else {print $1"_"$NF}}') | awk '{print $1"\t"$2}' > /Users/jmattic1/Documents/VirusToBeFixed.tsv

	rm -f /Users/jmattic1/Documents/FinishedFastasVirus.txt
	touch /Users/jmattic1/Documents/FinishedFastasVirus.txt
	while read line
	do
		fasta=$(echo $line | awk '{print $1}')
		presence=$(cat /Users/jmattic1/Documents/FinishedFastasVirus.txt | awk -v fasta="$fasta" '$1 == fasta {print $0}')
		if [[ $presence == "" ]]
		then
			genus=$(echo $line | awk '{print $2}')
			cat /Users/jmattic1/Documents/VirusToBeFixed.tsv | awk -v genus="$genus" '$2 == genus {print $0}' > /Users/jmattic1/Documents/subVirusToBeFixed.tsv
			duplicates=$(cat /Users/jmattic1/Documents/subVirusToBeFixed.tsv | wc -l)
			if [[ $duplicates -gt 1 ]]
			then
				count=1
				while read change
				do
					newgenus=$(echo $change | awk '{print $2}')
					newfasta=$(echo $change | awk '{print $1}')
					echo -e $newfasta"\t"$newgenus"-L"$count".faa" >> /Users/jmattic1/Documents/FinishedFastasVirus.txt

					mv /Users/jmattic1/Documents/ViralBLAST_FASTA_modded/$newfasta /Users/jmattic1/Documents/ViralBLAST_FASTA_modded/$newgenus"-L"$count".faa"
					count=$(( $count + 1 ))
				done < /Users/jmattic1/Documents/subVirusToBeFixed.tsv
			else
				echo -e $fasta"\t"$genus"-L1.faa" >> /Users/jmattic1/Documents/FinishedFastasVirus.txt
				mv /Users/jmattic1/Documents/ViralBLAST_FASTA_modded/$fasta /Users/jmattic1/Documents/ViralBLAST_FASTA_modded/$genus"-L1.faa"

			fi
		fi
	done < /Users/jmattic1/Documents/VirusToBeFixed.tsv

	for i in $(\ls /Users/jmattic1/Documents/ViralBLAST_FASTA_modded/ | grep ".faa")
	do
		species=$(echo $i | sed 's/\.faa//g')
		#if [[ $(cat /Users/jmattic1/Documents/ArchBactBLAST_DB/CombinedBactArch.faa | grep $species | wc -l) -eq 0 ]]
		#then
			cat /Users/jmattic1/Documents/ViralBLAST_FASTA_modded/$i | awk '{print $1}' | sed "s/>/>lcl\|$species\_/g" | awk '{if($0 ~ /^>/)  print substr($0,1,50); else print $0;}' >> /Users/jmattic1/Documents/Virus_DB/Viral.faa
		#fi
	done
	/Users/jmattic1/miniconda3/bin/makeblastdb -in /Users/jmattic1/Documents/Virus_DB/Viral.faa -dbtype 'prot' -parse_seqids



fi

####JUST Archaea
if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/ArchBLAST_DB/
	mkdir /Users/jmattic1/Documents/ArchBLAST_DB/
	
	cd /Users/jmattic1/Documents/EukProtProteinDatabase/


	touch /Users/jmattic1/Documents/ArchBactBLAST_DB/Arch.faa

	rm -rf /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/
	cp -r /Users/jmattic1/Documents/ArchaealBLAST_FASTA/ /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/

	paste -d "\t" <(\ls /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/ | grep ".faa") <(\ls /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/ | grep ".faa"  | sed 's/\.faa//g'| sed 's/Candidatus_/C/g' | sed 's/candidate_division/c/g' | sed 's/_archaeon_//g' | sed 's/_big_fil_rev_//g' | sed 's/archaeota//g' | sed 's/miscellaneous_//g' | awk '{print substr($0,1,27)}' | awk -F "_" '{if($1 == $NF){print $1"_arch"} else {print $1"_"$NF}}') | awk '{print $1"\t"$2}' > /Users/jmattic1/Documents/ArchaeaToBeFixed.tsv

	rm -f /Users/jmattic1/Documents/FinishedFastas.txt
	touch /Users/jmattic1/Documents/FinishedFastas.txt
	while read line
	do
		fasta=$(echo $line | awk '{print $1}')
		presence=$(cat /Users/jmattic1/Documents/FinishedFastas.txt | awk -v fasta="$fasta" '$1 == fasta {print $0}')
		if [[ $presence == "" ]]
		then
			genus=$(echo $line | awk '{print $2}')
			cat /Users/jmattic1/Documents/ArchaeaToBeFixed.tsv | awk -v genus="$genus" '$2 == genus {print $0}' > /Users/jmattic1/Documents/subArchaeaToBeFixed.tsv
			duplicates=$(cat /Users/jmattic1/Documents/subArchaeaToBeFixed.tsv | wc -l)
			if [[ $duplicates -gt 1 ]]
			then
				count=1
				while read change
				do
					newgenus=$(echo $change | awk '{print $2}')
					newfasta=$(echo $change | awk '{print $1}')
					echo -e $newfasta"\t"$newgenus"-L"$count".faa" >> /Users/jmattic1/Documents/FinishedFastas.txt

					mv /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/$newfasta /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/$newgenus"-L"$count".faa"
					count=$(( $count + 1 ))
				done < /Users/jmattic1/Documents/subArchaeaToBeFixed.tsv
			else
				echo -e $fasta"\t"$genus"-L1.faa" >> /Users/jmattic1/Documents/FinishedFastas.txt
				mv /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/$fasta /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/$genus"-L1.faa"

			fi
		fi
	done < /Users/jmattic1/Documents/ArchaeaToBeFixed.tsv

	for i in $(\ls /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/ | grep ".faa")
	do
		species=$(echo $i | sed 's/\.faa//g')
		#if [[ $(cat /Users/jmattic1/Documents/ArchBactBLAST_DB/CombinedBactArch.faa | grep $species | wc -l) -eq 0 ]]
		#then
			cat /Users/jmattic1/Documents/ArchaealBLAST_FASTA_modded/$i | awk '{print $1}' | sed "s/>/>lcl\|$species\_/g" | awk '{if($0 ~ /^>/)  print substr($0,1,75); else print $0;}' >> /Users/jmattic1/Documents/ArchBLAST_DB/Arch.faa
		#fi
	done
	/Applications/miniconda3/bin/makeblastdb -in /Users/jmattic1/Documents/ArchBLAST_DB/Arch.faa -dbtype 'prot' -parse_seqids

fi

####All Bacteria
if [[ 1 -eq 2 ]]
then
	rm -rf /Users/jmattic1/Documents/BactBLAST_DB/
	mkdir /Users/jmattic1/Documents/BactBLAST_DB/
	


	###Problem bacteria -- multispecies and sp. Solved, hopefully
	#cat /Users/jmattic1/Documents/BactRefSeq2/bacteria.wp_protein.1.protein.faa | sed 's/\[/yYyYyYy/g' | sed 's/\]/zZzZzZz/g' | awk '{if($0 ~ /^>/ && $0 ~ /MULTISPECIES/) { print ">"$NF"_multispecies_"$1;} else if($0 ~ /^>/) { for (i=1; i<=NF; ++i) { if ($i ~ "yYyYyYy") first=i }; for (i=1; i<=NF; ++i) { if ($i ~ "zZzZzZz") second=i } print ">"$first"_"$second"_"$1;} else {print $1}}' | sed 's/yYyYyYy//g' | sed 's/zZzZzZz//g' | sed 's/_>/_/g'
#cat /Users/jmattic1/Documents/BactRefSeq2/bacteria.wp_protein.*.protein.faa | awk '{if($0 ~ /^>/)  print ">"$(NF-1)"_"$NF"_"$1; else print $1}'


	cat /Users/jmattic1/Documents/BactRefSeq2/bacteria.wp_protein.*.protein.faa | sed 's/\[/yYyYyYy/g' | sed 's/\]/zZzZzZz/g' | awk '{if($0 ~ /^>/ && $0 ~ /MULTISPECIES/) { for (i=1; i<=NF; ++i) { if ($i ~ "yYyYyYy") first=i }; print ">"$first"_multispecies_"$1;} else if($0 ~ /^>/ && $0 ~ /endosymbiont/) { for (i=1; i<=NF; ++i) { if ($i ~ "zZzZzZz") second=i }; print ">"$second"_endosymbiont_"$1;} else if($0 ~ /^>/) { for (i=1; i<=NF; ++i) { if ($i ~ "yYyYyYy") first=i }; for (i=1; i<=NF; ++i) { if ($i ~ "zZzZzZz") second=i } print ">"$first"_"$second"_"$1;} else {print $1}}' | sed 's/yYyYyYy//g' | sed 's/zZzZzZz//g' | sed 's/_>/_/g' | sed 's/\//-/g' | sed 's/,/-/g' | sed 's/>/>lcl\|/g' > /Users/jmattic1/Documents/BactBLAST_DB/Bact.faa
	/Applications/miniconda3/bin/makeblastdb -in /Users/jmattic1/Documents/BactBLAST_DB/Bact.faa -dbtype 'prot' -parse_seqids

	cat /Users/jmattic1/Documents/BactBLAST_DB/Bact.faa | sed "s/\'//g" | sed 's/(.*)//g' > /Users/jmattic1/Documents/BactBLAST_DB/Bact.noquotes.faa
	/Applications/miniconda3/bin/makeblastdb -in /Users/jmattic1/Documents/BactBLAST_DB/Bact.noquotes.faa -dbtype 'prot' -parse_seqids

	faadir="/Users/jmattic1/Documents/RibosomalProteins_HVolc/"
	i=RPL8.HVolc.faa
	outdir="/Users/jmattic1/Documents/BactRPL_Alignment/"
	DB=/Users/jmattic1/Documents/BactBLAST_DB/Bact.faa

	###Use Fasttree first

	cat /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.faa | sed "s/\'//g" | sed 's/(.*)//g' > /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.noquotes.faa
	FastTree /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.noquotes.faa > /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.fasttree.treefile
	#Condense with R
	#outdir<-"/Users/jmattic1/Documents/BactRPL_Alignment/"
	#filename<-"/Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.fasttree.treefile"
	#name<-"AllBact.fasttree"
	#edgeLimit<-10
	#Manchar
	#nohup /usr/local/bin/Rscript RibosomalTreeReduction_basalnode_migratory.R > RibosomalTreeReduction_basalnode_migratory.log

	##Fixing the table -- one off, database has been edited
	cat /Users/jmattic1/Documents/BactRPL_Alignment/AllBact.fasttree.basalnode.1.tsv | sed 's/Burkholderia_L27/Burkholderia_L27_WP_158898489.1/g' | sed 's/Synechococcus_JA-2-3Ba/Synechococcus_JA-2-3Ba_WP_011434155.1/g' > /Users/jmattic1/Documents/BactRPL_Alignment/AllBact.fasttree.basalnode.1.fixed.tsv 

	###Subset the DB based on the nodes -- run IQTree for a better tree for more condensation

	conda activate /Applications/miniconda3/envs/samtools

	###Loop:
	rm  /Users/jmattic1/Documents/BactRPL_Alignment/Representative.fasta.list
	tail -n+2 /Users/jmattic1/Documents/BactRPL_Alignment/AllBact.fasttree.basalnode.1.fixed.tsv | awk '{print $1}' > /Users/jmattic1/Documents/BactRPL_Alignment/Representative.list
	while read line
	do
		hitcount=$(cat /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.noquotes.faa | grep $line | wc -l)
		if [[ $hitcount -gt 1 ]]
		then
			echo $line" is messed up -- multi-hit"
		elif [[ $hitcount -lt 1 ]]
		then
			echo $line" is messed up -- no hit"
		else
			cat /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.noquotes.faa | grep $line | sed 's/>//g' >> /Users/jmattic1/Documents/BactRPL_Alignment/Representative.fasta.list
		fi
	done < /Users/jmattic1/Documents/BactRPL_Alignment/Representative.list

	Subsetvars=$(cat /Users/jmattic1/Documents/BactRPL_Alignment/Representative.fasta.list | tr '\n' ' ')
	samtools faidx /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.noquotes.faa

	samtools faidx /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.noquotes.faa $Subsetvars > /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.noquotes.subset.faa
	##Redo MAFFT and trimal

	/Applications/miniconda3/bin/mafft --anysymbol --reorder --auto /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.noquotes.subset.faa > /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.noquotes.subset.realign.faa
	#Trimal to trim the alignment: gappyout parameter errs on the side of keeping gaps
	/Applications/miniconda3/bin/trimal -in /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.align.trimal.highcore.constrained.noquotes.subset.realign.faa -out /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.partiallyreduced.faa -automated1
	/Applications/miniconda3/bin/iqtree -T 14 -s /Users/jmattic1/Documents/BactRPL_Alignment/RPL8.HVolc.blast.1e-5.partiallyreduced.faa --alrt 1000 -B 1000
fi

##Bact from Proteins -- doesn't work
if [[ 1 -eq 2 ]]
then
	##Run R script again -- take cutoff of 7
	##Manually edit in our 14 species where possible
	##Replace any candidatus or sp. or strains with full species where possible
	##Pull from database
	rm -rf /Applications/Databases/Bact_Condensed_v4/
	mkdir /Applications/Databases/Bact_Condensed_v4/
	#conda activate samtools
	#samtools faidx /Users/jmattic1/Documents/BactBLAST_DB/Bact.faa
	blastDB=/Users/jmattic1/Documents/BactBLAST_DB/Bact.noquotes.faa
	guidetable=$(echo "/Users/jmattic1/Documents/BactRPL_Alignment/AllBact.iqtree.basalnode.7.Redited.manual.tsv")
	greplist=$(tail -n+2 $guidetable | awk '{print $1}' | awk -F "_" '{print $1"_"$2"_"}' | tr '\n' '|' | sed 's/\|$//g')
	cat /Users/jmattic1/Documents/BactBLAST_DB/Bact.noquotes.faa | grep -E $greplist | sed 's/>//g' > /Applications/Databases/Bact_Condensed_v4/tempgenelist.txt
	#for i in $(tail -n+2 $guidetable | awk '{print $1}' | awk -F "_" '{print $1"_"$2}')
	#do
	#	cat /Users/jmattic1/Documents/BactBLAST_DB/Bact.noquotes.faa | grep $i"_" | sed 's/>//g' >> /Applications/Databases/Bact_Condensed_v4/tempgenelist.txt

	#done
	cat /Applications/Databases/Bact_Condensed_v4/tempgenelist.txt | awk -F "_" '{print $1"_"$2}' | sort | uniq -c > /Applications/Databases/Bact_Condensed_v4/ProteomeInventory.txt

	/Applications/miniconda3/bin/blastdbcmd -db $blastDB -dbtype prot -entry_batch /Applications/Databases/Bact_Condensed_v4/tempgenelist.txt -outfmt %f -out /Applications/Databases/Bact_Condensed_v4/Bact_Condensed_v4.faa

fi

###Bact from downloads -- unzip
if [[ 1 -eq 1 ]]
then
	##Fix directory structure
	rm -rf /Applications/Databases/Bact_ManualDownloadv2/
	mkdir /Applications/Databases/Bact_ManualDownloadv2/
	\ls /Applications/Databases/Bact_ManualDownload/ | grep "zip" > /Applications/Databases/Bactzip.txt
	while read line
	do
		newfile=$(echo $line | tr ' ' '_')
		cp /Applications/Databases/Bact_ManualDownload/$line /Applications/Databases/Bact_ManualDownloadv2/$newfile
		newdir=$(echo $newfile | sed 's/.zip//g')
		mkdir /Applications/Databases/Bact_ManualDownloadv2/$newdir
		cd /Applications/Databases/Bact_ManualDownloadv2/$newdir
		unzip /Applications/Databases/Bact_ManualDownloadv2/$newfile
	done < /Applications/Databases/Bactzip.txt

fi

###Bact from downloads -- create db
if [[ 1 -eq 2 ]]
then
	##Run R script again -- take cutoff of 7
	rm -rf /Applications/Databases/Bact_Condensed_v4/
	mkdir /Applications/Databases/Bact_Condensed_v4/
	rm -rf /Users/jmattic1/Documents/Bact_Condensed_v4/
	mkdir /Users/jmattic1/Documents/Bact_Condensed_v4/
	#conda activate samtools


	#samtools faidx /Users/jmattic1/Documents/BactBLAST_DB/Bact.faa
	cat /Users/jmattic1/Documents/EukCondArchFullBactProteomeRGBlabel.v3.csv | grep ",Bacteria," > /Users/jmattic1/Documents/EukCondArchFullBactProteomeRGBlabel.v3.bact.csv
	##Bact is all NCBI
	while read line
	do
		species=$(echo $line | awk -F ',' '{print $1}')

		echo $species
		GCF=$(echo $line | awk -F ',' '{print $3}')
		faafile=$(\ls /Applications/Databases/Bact_ManualDownloadv2/*/*/*/$GCF/protein.faa)
		if [[ $(echo $faafile | grep "GCA_" | wc -l) -gt 0 || $(echo $faafile | grep "GCF_" | wc -l) -gt 0 ]]
		then
			cat $faafile | sed -e "s/>/>"$species"_/g" | sed 's/ /_/g' | sed 's/\_\[.*//g' >> /Users/jmattic1/Documents/Bact_Condensed_v4/Bact_Condensed_v4.faa
		else
			echo "WTF"
			break
		fi

	done < /Users/jmattic1/Documents/EukCondArchFullBactProteomeRGBlabel.v3.bact.csv

	cat /Users/jmattic1/Documents/Bact_Condensed_v4/Bact_Condensed_v4.faa | sed '/^$/d' | awk '{print $1}' | sed 's/,.*//g' | sed 's/\|/_/g' | sed 's/>/>lcl\|/g' | awk '{if($0 ~ /^>/)  print substr($0,1,50); else print $0;}' > /Applications/Databases/Bact_Condensed_v4/Bact_Condensed_v4.faa
	
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/Bact_Condensed_v4/Bact_Condensed_v4.faa -dbtype 'prot' -parse_seqids
	



fi



###Final database assembly
if [[ 1 -eq 2 ]]
then
	rm -rf /Applications/Databases/ArchBact_Condensed_v3/
	mkdir /Applications/Databases/ArchBact_Condensed_v3/
	cd /Applications/Databases/ArchBact_Condensed_v3/
	guidetable=$(echo "/Users/jmattic1/Documents/SR4concatlist.v3.tsv")
	#conda activate samtools
	for i in $(cat $guidetable | awk '{print $1}')
	do
		cat /Applications/Databases/ArchBactBLAST_DB_v2_xAsgard/CombinedBactArchxAsgard.faa | grep $i"_" | sed 's/>//g' > /Applications/Databases/ArchBact_Condensed_v3/tempgenelist.txt
		cat /Applications/Databases/ArchBact_Condensed_v3/tempgenelist.txt | awk -F "_" '{print $1"_"$2}' | sort | uniq -c >> /Applications/Databases/ArchBact_Condensed_v3/ProteomeInventory.txt
		samtoolsvar=$(cat /Applications/Databases/ArchBact_Condensed_v3/tempgenelist.txt | tr '\n' ' ')
		samtools faidx /Applications/Databases/ArchBactBLAST_DB_v2_xAsgard/CombinedBactArchxAsgard.faa $samtoolsvar >> /Applications/Databases/ArchBact_Condensed_v3/ArchBact_Condensed_v3.faa
	done
	rm -rf /Applications/Databases/Euk_ArchBact_Condensed_v3.1/
	mkdir /Applications/Databases/Euk_ArchBact_Condensed_v3.1/
	cd /Applications/Databases/Euk_ArchBact_Condensed_v3.1/
	cat /Applications/Databases/ArchBact_Condensed_v3/ArchBact_Condensed_v3.faa /Applications/Databases/EukBLASTDB/Euk.protein.db.faa > /Applications/Databases/Euk_ArchBact_Condensed_v3.1/Euk_ArchBact_Condensed_v3.1.faa
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/Euk_ArchBact_Condensed_v3.1/Euk_ArchBact_Condensed_v3.1.faa -dbtype 'prot' -parse_seqids

	rm -rf /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/
	mkdir /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/
	cd /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/
	cat /Applications/Databases/Euk_ArchBact_Condensed_v3.1/Euk_ArchBact_Condensed_v3.1.faa /Applications/Databases/Virus_DB/Viral.faa > /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/Euk_Virus_ArchBact_Condensed_v3.1/Euk_Virus_ArchBact_Condensed_v3.1.faa -dbtype 'prot' -parse_seqids

	###Final Database with fixed Bact
	##Use condensed Arch
	rm -rf /Applications/Databases/Arch_Condensed_v3/
	mkdir /Applications/Databases/Arch_Condensed_v3/
	grepspecieslist=$(cat /Applications/Databases/ArchBact_Condensed_v3/ArchBact_Condensed_v3.faa | grep ">" | awk -F "_" '{print $1"_"$2}' | sort | uniq | grep -E "L|arch" | sed 's/>lcl\|//g' | tr '\n' '|' | sed 's/\|$//g')
	cat /Applications/Databases/ArchBact_Condensed_v3/ArchBact_Condensed_v3.faa | grep ">" | grep -E $grepspecieslist | sed 's/>//g' > /Applications/Databases/Arch_Condensed_v3/genelist.txt
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/ArchBact_Condensed_v3/ArchBact_Condensed_v3.faa -dbtype 'prot' -parse_seqids
	/Applications/miniconda3/bin/blastdbcmd -db /Applications/Databases/ArchBact_Condensed_v3/ArchBact_Condensed_v3.faa -dbtype prot -entry_batch /Applications/Databases/Arch_Condensed_v3/genelist.txt -outfmt %f -out /Applications/Databases/Arch_Condensed_v3/Arch_Condensed_v3.faa

	rm -rf /Applications/Databases/ArchBact_Condensed_v4/
	mkdir /Applications/Databases/ArchBact_Condensed_v4/
	cat /Applications/Databases/Arch_Condensed_v3/Arch_Condensed_v3.faa /Applications/Databases/Bact_Condensed_v4/Bact_Condensed_v4.faa > /Applications/Databases/ArchBact_Condensed_v4/ArchBact_Condensed_v4.faa
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/ArchBact_Condensed_v4/ArchBact_Condensed_v4.faa -dbtype 'prot' -parse_seqids
	
	rm -rf /Applications/Databases/ArchBact_Condensed_v4.1/
	mkdir /Applications/Databases/ArchBact_Condensed_v4.1/
	cat <(cat /Applications/Databases/Arch_Condensed_v3/Arch_Condensed_v3.faa | sed 's/>/>lcl\|/g') /Applications/Databases/Bact_Condensed_v4/Bact_Condensed_v4.faa > /Applications/Databases/ArchBact_Condensed_v4.1/ArchBact_Condensed_v4.1.faa
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/ArchBact_Condensed_v4.1/ArchBact_Condensed_v4.1.faa -dbtype 'prot' -parse_seqids


	rm -rf /Applications/Databases/Euk_ArchBact_Condensed_v4/
	mkdir /Applications/Databases/Euk_ArchBact_Condensed_v4/
	cat /Applications/Databases/ArchBact_Condensed_v4/ArchBact_Condensed_v4.faa /Applications/Databases/EukBLASTDB/Euk.protein.db.faa > /Applications/Databases/Euk_ArchBact_Condensed_v4/Euk_ArchBact_Condensed_v4.faa
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/Euk_ArchBact_Condensed_v4/Euk_ArchBact_Condensed_v4.faa -dbtype 'prot' -parse_seqids

	rm -rf /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4/
	mkdir /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4/
	cat /Applications/Databases/Euk_ArchBact_Condensed_v4/Euk_ArchBact_Condensed_v4.faa /Applications/Databases/Virus_DB/Viral.faa > /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4/Viral_Euk_ArchBact_Condensed_v4.faa
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4/Viral_Euk_ArchBact_Condensed_v4.faa -dbtype 'prot' -parse_seqids

	rm -rf /Applications/Databases/Euk_ArchBact_Condensed_v4.1/
	mkdir /Applications/Databases/Euk_ArchBact_Condensed_v4.1/
	cat /Applications/Databases/ArchBact_Condensed_v4.1/ArchBact_Condensed_v4.1.faa /Applications/Databases/EukBLASTDB/Euk.protein.db.faa > /Applications/Databases/Euk_ArchBact_Condensed_v4.1/Euk_ArchBact_Condensed_v4.1.faa
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/Euk_ArchBact_Condensed_v4.1/Euk_ArchBact_Condensed_v4.1.faa -dbtype 'prot' -parse_seqids

	rm -rf /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1/
	mkdir /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1/
	cat /Applications/Databases/Euk_ArchBact_Condensed_v4.1/Euk_ArchBact_Condensed_v4.1.faa /Applications/Databases/Virus_DB/Viral.faa > /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1/Viral_Euk_ArchBact_Condensed_v4.1.faa
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1/Viral_Euk_ArchBact_Condensed_v4.1.faa -dbtype 'prot' -parse_seqids

	rm -rf mkdir /Applications/Databases/All_Arch_v4/
	mkdir /Applications/Databases/All_Arch_v4/
	cat /Users/jmattic1/Documents/ArchaealBLAST_FASTA/*.faa /Users/jmattic1/Documents/Laura_proteomes/*.faa > /Applications/Databases/All_Arch_v4/All_Arch_v4.faa
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Databases/All_Arch_v4/All_Arch_v4.faa -dbtype 'prot' -parse_seqids


fi



###Metazoa Assembly
if [[ 1 -eq 2 ]]
then
	rm -rf /Applications/Metazoa_DB/
	mkdir /Applications/Metazoa_DB/
	for i in $(\ls ~/Downloads/ProcB_metazoa_DataDryad/Original_Proteomes/*.fa | grep "fa")
	do
		orgname=$(echo $i | awk -F "/" '{print $NF}' | sed 's/.fa//g')
		#cat $i | grep ">" | awk '{print $1}' | awk -F "|" '{print $1}' | sed 's/>//g' | awk -v orgname="$orgname" '{print ">"orgname"_"$0}' | head -1
		cat $i | awk '{print $1}' | sed 's/*//g' | sed 's/\|/_/g' | awk -F "/" '{print $1}' | sed 's/Glottidia_pyramidata//g' | sed 's/Proneomenia_custodiens//g' | sed 's/Leptochiton_rugatus//g' | sed 's/Phoronis_psammophila//g' | sed 's/Hemithiris_psittacea//g' | sed 's/Macandrevia_cranium//g' | sed 's/Novocrania_anomala//g' | sed "s/>/>lcl\|$orgname\_/g" | awk '{if($0 ~ /^>/)  print substr($0"_"NR,1,75); else print $0;}' >> /Applications/Metazoa_DB/Metazoa.faa

	done
	/Applications/miniconda3/bin/makeblastdb -in /Applications/Metazoa_DB/Metazoa.faa -dbtype 'prot' -parse_seqids

fi
