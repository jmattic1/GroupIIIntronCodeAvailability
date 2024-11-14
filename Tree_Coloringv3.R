require(ape)
require(TreeTools)
require(ggplot2)
require(treeio)
require(castor)
require(phytools)
require(readr)
require(dendextend)
require(phylobase)
require(adephylo)

args = commandArgs(trailingOnly=TRUE)
#filename=args[1]
filename="/Users/jmattic1/Documents/TERT_reanalysis/Prp8concat/PFAM.1e-20Filter.allasgard.allbact.align.trim.faa.treefile"

outfile<-gsub(".contree",".pdf",filename)
outfile<-gsub(".treefile",".pdf",filename)

myTree <- ape::read.tree(filename)
myTree<- ladderize(myTree, right = FALSE)
#myTree_text <- as.character(readr::read_delim(filename, delim = "\t",col_names = F))

#translationTable<-read.csv("/Users/jmattic1/Documents/EukCondArchBactProteomeRGBlabel.csv")
translationTable<-read.csv("/Users/jmattic1/Documents/EukCondArchFullBactProteomeRGBlabel.v4.csv")

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

cols<-c()
###Tip colors
for (x in 1:length(myTree$tip.label)) {
  if(length(translationTable[translationTable$actualTaxon == species[x],]$Supergroup) > 0){
    tipcolor<-tolower(rgb(translationTable[translationTable$actualTaxon == species[x],]$R,translationTable[translationTable$actualTaxon == species[x],]$G,translationTable[translationTable$actualTaxon == species[x],]$B, maxColorValue=255))
    cols<-c(cols,tipcolor)
  }else{
    cols<-c(cols,"#000000")
  }
}
select.tip.or.node <- function(element, tree) {
  ifelse(element < Ntip(tree)+1,
         tree$tip.label[element],
         tree$node.label[element-Ntip(tree)])
}
###Branch colors
edge_table <- data.frame(
  "parent" = myTree$edge[,1],
  "par.name" = sapply(myTree$edge[,1],
                      select.tip.or.node,
                      tree = myTree),
  "child" = myTree$edge[,2],
  "chi.name" = sapply(myTree$edge[,2],
                      select.tip.or.node,
                      tree = myTree)
)

if(length(edge_table[edge_table$chi.name == "",]$chi.name) > 0){
  edge_table[edge_table$chi.name == "",]$chi.name<-"0/0"
}

edgecols<-c()
#Physcomitrella p XP 024388226.1 uncharacteriz
for (x in 1:length(edge_table$parent)) {
  subedges<-edge_table[x,]
  species_end<-paste(strsplit(subedges$chi.name,"_", fixed=TRUE)[[1]][1],strsplit(subedges$chi.name, "_", fixed=TRUE)[[1]][2],sep="_")

  if(length(translationTable[translationTable$actualTaxon == species_end,]$Supergroup) > 0){
    edgecolor<-tolower(rgb(translationTable[translationTable$actualTaxon == species_end,]$R,translationTable[translationTable$actualTaxon == species_end,]$G,translationTable[translationTable$actualTaxon == species_end,]$B, maxColorValue=255))
    edgecols<-c(edgecols,edgecolor)
    #edgecols[subedges$child] = edgecolor
  # }else{
  #    edgecolor<-"#000000"
  #    edgecols<-c(edgecols,edgecolor)
  # }
  }else{
    ##Track monophyly

    allparents<-edge_table[edge_table$parent == edge_table[x,]$parent,]
    while(length(allparents[grep("/",allparents$chi.name),]$parent) > 0){
      subparents<-allparents[grep("/",allparents$chi.name),]
      newchildren<-edge_table[edge_table$parent %in% subparents$child,]
      allparents<-allparents[grep("/",allparents$chi.name,invert=T),]
      allparents<-rbind(allparents,newchildren)
    }
    monophyly<-T
    count<-1
    monocolor<-c()
    while(monophyly & count < length(allparents$parent)){
      perspecies<-paste(strsplit(allparents[count,]$chi.name,"_", fixed=TRUE)[[1]][1],strsplit(allparents[count,]$chi.name, "_", fixed=TRUE)[[1]][2],sep="_")
      if(length(translationTable[translationTable$actualTaxon == perspecies,]$Supergroup) > 0){
        indcolor<-tolower(rgb(translationTable[translationTable$actualTaxon == perspecies,]$R,translationTable[translationTable$actualTaxon == perspecies,]$G,translationTable[translationTable$actualTaxon == perspecies,]$B, maxColorValue=255))
      }else{
        indcolor<-"#000000"
      }
      monocolor<-c(indcolor,monocolor)
      if(length(unique(monocolor)) > 1){
        monophyly<-F
      }
      count<-count+1
    }
    if(monophyly){
      edgecols<-c(edgecols,unique(monocolor))
    }else{
      edgecols<-c(edgecols,"#000000")
    }
  }
  
}

pdf(outfile,height = 50,width = 20)
print(plot(myTree,edge.color=edgecols,tip.color=cols,cex=0.3,show.node.label = TRUE,type="unrooted"))
#print(plot(myTree,edge.color=edgecols,tip.color=cols,cex=0.3,show.node.label = F,show.tip.label = F,type="unrooted"))
add.scale.bar(length = 0.5)
dev.off()

