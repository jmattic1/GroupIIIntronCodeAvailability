###Arg 1 is the table
###Arg 2 is the number of candidates
###Arg 3 is the threshold
###Arg 4 is the blastfile to be filtered
###Arg 5 is the output directory
library.path <- .libPaths()
args=commandArgs(trailingOnly=TRUE)

input=args[1]
#input="/Users/jmattic1/Documents/blastwrappr_test/XP_726415.1PlasmodiumyoeliiSpo112_Homo_s.faa.dist"
#input="/Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/GroupII_Introns_BLAST_10998_running/G2I_Necator_a_Tetraselmis_s.dist"
#print(args[2])
len=as.numeric(args[2])
#len=7
#len=1850
threshold=as.numeric(args[3])
#threshold=0.002
#threshold=10
blastfilter=args[4]
#blastfilter="/Users/jmattic1/Documents/blastwrappr_test/tempblast.eval.blast"
#blastfilter="/Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/GroupII_Introns_BLAST_10998_running/tempblast.eval.blast"
outdir=args[5]
#outdir="/Users/jmattic1/Documents/blastwrappr_test/"
#outdir="/Users/jmattic1/Documents/TERT_reanalysis/new_db_first_pass/GroupII_Introns_BLAST_10998_running/"

distmatrix<-read.table(input,sep = ",",col.names = paste0("V",seq_len(len)),fill = TRUE)
for (x in 2:len) {
  distmatrix[x,x:len]<-distmatrix[x,1:(len-x+1)]
  distmatrix[x,1:(x-1)]<-NA
}
outputframe<-data.frame()
family=1
currentmembers<-1
for (x in 1:len) {
  if(!(x %in% currentmembers)){
    suboutputframe<-data.frame(family,currentmembers)
    colnames(suboutputframe)<-c("Family","Members")
    outputframe<-unique(rbind(outputframe,suboutputframe))
    if(x %in% outputframe$Members){
      family<-unique(outputframe[outputframe$Members == x,]$Family)
      currentmembers<-outputframe[outputframe$Family == family,]$Members
    }else{
      family<-max(outputframe$Family)+1
      currentmembers<-x
    }
  }
  if(family == 5){
    break
  }
  firstseqcoord=x
  suboutputframe<-data.frame(family,firstseqcoord)
  colnames(suboutputframe)<-c("Family","Members")
  for (y in x:len) {
    secondseqcoord=y
    if(!(secondseqcoord %in% outputframe$Members)){
      if(distmatrix[x,y] <= threshold){
        currentmembers<-unique(c(currentmembers,secondseqcoord))
      }else if(x == (len-1)){
        suboutputframe<-data.frame(max(c(outputframe$Family,family))+1,secondseqcoord)
        colnames(suboutputframe)<-c("Family","Members")
        outputframe<-rbind(outputframe,suboutputframe)
        
      }
    }
  }
}
suboutputframe<-data.frame(family,currentmembers)
colnames(suboutputframe)<-c("Family","Members")
outputframe<-unique(rbind(outputframe,suboutputframe))

blasttable<-read.table(blastfilter)
finalblasttable<-data.frame()
totalfamilies<-unique(outputframe$Family)
for (x in 1:length(totalfamilies)) {
  currentseqs<-outputframe[outputframe$Family == totalfamilies[x],]$Members
  subblasttable<-blasttable[currentseqs,]
  subblasttable<-subblasttable[order(subblasttable$V12,subblasttable$V10,decreasing = T),]
  finalblasttable<-rbind(finalblasttable,subblasttable[1,])
}

write.table(finalblasttable,file=paste(outdir,"/tempblast.eval.post.blast",sep=""),quote=F,col.names = F, row.names = F,sep="\t")