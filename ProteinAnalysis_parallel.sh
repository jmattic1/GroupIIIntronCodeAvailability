for ((var=1; var<=$#; var++))
do
	if [[ "${!var}" = "-mod" ]] 
		then
		var=$(($var+1))
		mod=${!var}
		echo "Variable mod is "$mod
	fi
done


#for i in $(cat dbv2_alien_index.10min.faa | grep ">" | sed 's/>//g')
#do
#samtools faidx dbv2_alien_index.10min.faa $i > ./Wesley_Mesodinium_Proteome/$i.faa
#done

DB="/Applications/Databases/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard.faa"
#outdir="/Users/jmattic1/Documents/RibosomalProteinAlignmentsHVolc_reducedset_v3.1/"
#outdir="/Users/jmattic1/Documents/BanooUPs_blast_v3.1/"
outdir="/Users/jmattic1/Documents/Viral_Euk_ArchBact_Condensed_v4.1_AllAsgard_Proteome_Output/"

#faadir="/Users/jmattic1/Documents/RibosomalProteins_HVolc/"
faadir="/Users/jmattic1/Documents/Banoo_Candidate_UPs/"

count=1
for i in $(\ls $faadir | grep "faa")
do
	eval=$(expr $count % 5)
	if [[ $eval -eq $mod ]]
	then
		#/Applications/BLAST_wrapper_notree.sh -fasta $faadir/$i -blastDB $DB -eval 1e-5 -filter 60 -type psiblast -outdir $outdir -addorigin no -reciprocal no -trimal -automated1 -SR4 yes -num_hits 1
		/Applications/BLAST_wrapper_notree.sh -fasta $faadir/$i -blastDB $DB -eval 1e-5 -filter 60 -type psiblast -outdir $outdir -addorigin yes -reciprocal no -trimal -gappyout -SR4 no -num_hits 1

	fi
	count=$(( $count + 1 ))
done
