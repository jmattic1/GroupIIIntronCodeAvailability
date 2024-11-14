require(ape)
require(TreeTools)
require(ggplot2)
require(treeio)
require(castor)
require(phytools)

#filename=args[1]
filename<-"/Users/jmattic1/Documents/EukandRiboprotein_ConcatTree_60perc_v4//EukRibosomalProteins.concatenated.final.faa.treefile"
filename<-"/Users/jmattic1/Documents/Ploidy/FirstPass//Rad54_BLAST_16843_running/Rad54.blast.1e-10.align.trimal.faa.treefile"

filename<-"/Users/jmattic1/Documents/Wesley_Mesodinium_Proteome_Output/g4460_e1e-5f60nh1ni6_96x496_14797/g4460.blast.1e-5.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/GameteFusionPhylogenies/FirstPass/GCS1_e1e-5f40nh1ni6_76x599_19074/GCS1.blast.1e-5.align.trimal.faa.treefile"
#filename<-"/Users/jmattic1/Documents/Mei4/FirstPass//Mei4_all_e1e-5f30nh1ni6_65x207_14570/Mei4_all.blast.1e-5.align.trimal.faa.treefile"
#filename<-"/Users/jmattic1/Documents/GameteFusionPhylogenies/GCS_ConcatTree_60perc/EukRibosomalProteins.concatenated.final.SR4.fas.treefile"
#filename<-"/Users/jmattic1/Documents/GameteFusionPhylogenies/FirstPass/Fsx1_e1e-5f40nh1ni6_84x248_276/Fsx1.blast.1e-5.align.trimal.faa.treefile"
#filename<-"/Users/jmattic1/Documents/EukandRiboprotein_ConcatTree_60perc/EukRibosomalProteins.concatenated.final.faa.treefile"
#filename<-"/Users/jmattic1/Documents/TERT_reanalysis/new_db_rna_analysis/GroupII_Intron_single_e1e-5f60nh1ni6_35x437_14745/GroupII_Intron_single.blast.1e-5.align.trimal.faa.treefile"
#filename<-"/Users/jmattic1/Documents/Prp8AF/PRP8_e1e-5f40nh1ni6_182x2293_25360/PRP8.blast.1e-5.align.trimal.faa.treefile"

filename<-"/Users/jmattic1/Documents/GameteFusionPhylogenies/Metazoa/GCS1.blast.1e-10.origin.nostar_e1e-5f40nh1ni6_138x499_32591/GCS1.blast.1e-10.origin.nostar.blast.1e-5.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/GameteFusionPhylogenies/NumHit5/GCS1_e1e-5f40nh10ni6_88x660_21575/GCS1.blast.1e-5.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass/Cdc6_e1e-5f40nh1ni6_1034x330_32655/Cdc6.blast.1e-5.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass/Orc1_e1e-5f40nh1ni6_555x263_4120/Orc1.blast.1e-5.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass/Orc2_e1e-5f40nh1ni6_188x340_911/Orc2.blast.1e-5.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass/Orc3_e1e-5f40nh1ni6_118x509_14540/Orc3.blast.1e-5.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass/Orc4_e1e-5f40nh1ni6_878x231_27261/Orc4.blast.1e-5.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass/Orc5_e1e-5f40nh1ni6_1030x510_7214/Orc5.blast.1e-5.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass/Orc6_e1e-5f40nh1ni6_396x283_19934/Orc6.blast.1e-5.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass/Cdc6_e1e-20f60nh6ni6_686x328_7873/Cdc6.blast.1e-20.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass/OrcParalogues_e1e-20f80nh6ni6_941x362_1595/OrcParalogues.blast.1e-20.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass/OrcParalogues_e1e-20f80nh1ni6_649x337_18916/OrcParalogues.blast.1e-20.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass//OrcParalogues_e1e-20f80nh20ni6_941x362_7206/OrcParalogues.blast.1e-20.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass//OrcParalogueswithHaloferax_e1e-20f80nh20ni6_1060x354_13055/OrcParalogueswithHaloferax.blast.1e-20.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass//Orc6_e1e-20f40nh20ni6_77x269_27777/Orc6.blast.1e-20.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass//OrcParalogues_HVolc_Protist_e1e-20f40nh20ni6_1410x217_28816/OrcParalogues_HVolc_Protist.blast.1e-20.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass//OrcParalogues_HVolc_Protist_e1e-20f80nh20ni6_1086x283_12461/OrcParalogues_HVolc_Protist.blast.1e-20.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Orc_analysis/FirstPass//OrcParalogues_HVolc_Protist_e1e-20f60nh20ni6_1256x487_26515/OrcParalogues_HVolc_Protist.blast.1e-20.align.trimal.faa.treefile"
filename<-"/Users/jmattic1/Documents/Muts.species.uniq.align.trim.faa.treefile"

outfile<-gsub(".contree",".supergroups.contree",filename)
outfile<-gsub(".treefile",".supergroups.treefile",filename)
outfile2<-gsub(".contree",".supergroups.csv",filename)
outfile2<-gsub(".treefile",".supergroups.csv",filename)

myTree <- ape::read.tree(filename)

#translationTable<-read.csv("/Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.csv")
translationTable<-read.csv("/Users/jmattic1/Documents/EukCondArchFullBactProteomeRGBlabel.v3.csv")

translationTable[translationTable$Supergroup == "DIscoba",]$Supergroup<-"Discoba"
translationTable$actualTaxon<-translationTable$Taxon

euktable<-read.delim("/Users/jmattic1/Documents/EukSpeciesKey.tsv",sep="\t",header=F,stringsAsFactors = F)

for (x in 1:length(euktable$V1)) {
  translationTable[translationTable$Taxon == euktable[x,]$V1,]$actualTaxon<-gsub("\\.","",euktable[x,]$V2)
}

DF <- data.frame(do.call(rbind, strsplit(myTree$tip.label, "_", fixed=TRUE)),stringsAsFactors = F)
DF <- DF[,c(1:2)]
DF[,1]<-as.character(DF[,1])
DF[,2]<-as.character(DF[,2])

species<-apply(DF, 1, paste, collapse="_")
referenceframe<-data.frame()

for (x in 1:length(myTree$tip.label)) {
  if(length(translationTable[translationTable$actualTaxon == species[x],]$Supergroup) > 0){
    #subreferenceframe<-data.frame(myTree$tip.label[x],translationTable[translationTable$actualTaxon == species[x],]$Supergroup,translationTable[translationTable$actualTaxon == species[x],]$Phylum,translationTable[translationTable$actualTaxon == species[x],]$Domain)
    subreferenceframe<-data.frame(myTree$tip.label[x],translationTable[translationTable$actualTaxon == species[x],]$Supergroup)
    colnames(subreferenceframe)<-c("Species","Supergroup")
  myTree$tip.label[x]<-paste(myTree$tip.label[x],translationTable[translationTable$actualTaxon == species[x],]$Supergroup,translationTable[translationTable$actualTaxon == species[x],]$Phylum,sep="_")
  referenceframe<-rbind(referenceframe,subreferenceframe)
  }else{
    #subreferenceframe<-data.frame(myTree$tip.label[x],"Virus","Virus","Virus")
    subreferenceframe<-data.frame(myTree$tip.label[x],"Virus")
    colnames(subreferenceframe)<-c("Species","Supergroup")
    
    myTree$tip.label[x]<-paste(myTree$tip.label[x],"Virus",sep="_")
    #myTree$tip.label[x]<-paste(myTree$tip.label[x],"Metazoa",sep="_")
    referenceframe<-rbind(referenceframe,subreferenceframe)
    
  }
  
}
write.tree(myTree,file=outfile)
write.table(referenceframe,file=outfile2,quote=F,col.names = T, row.names = F,sep=",")

