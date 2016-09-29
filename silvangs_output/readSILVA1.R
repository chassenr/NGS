#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements: 
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# reads SILVA output and creates for each sample a table with reference OTU name, abundance, taxonomic path, status (classified or NoRelative)
# this function should work with most NGS data
# for an example of complete SILVA pipeline in R: see \trunk\clearingHouse\SILVApipe.R
# input:
# uses the reference OTU names extracted from the .fna files in SILVA\results\ssu\exports\otu_references
# extraction of names best done on unix:
#  grep '>' Project---ssu---otu_references---Sample----Classified.fna|cat>Sample_classified.txt
#  grep '>' Project---ssu---otu_references---Sample----NoRelative.fna|cat>Sample_NoRelative.txt
# Usage of function:
#  e.g. Sample=readSILVA("Sample_classified.txt", "Sample_NoRelative.txt")
# output:
# table with reference OTU name, abundance, taxonomic path, status (classified or NoRelative)
# for further information on status: see SILVA documentation
# dependencies:
# none
#=============================================================================
#documentation end


readSILVA1=function(class,NoRel){
  a=read.csv(class,sep="\t",h=F)
  a1=sapply(as.character(a[,1]),function(x){strsplit(x,">")[[1]][2]})
  a2=as.factor(sapply(a1,function(x){strsplit(x,"\\.*\\.")[[1]][1]}))
  
  b=read.csv(NoRel,sep="\t",h=F)
  b1=sapply(as.character(b[,1]),function(x){strsplit(x,">")[[1]][2]})
  b2=as.factor(sapply(b1,function(x){strsplit(x,"\\.*\\.")[[1]][1]}))
  
  aa=data.frame(a2,a[,2:3],rep("class",length(a2)))
  bb=data.frame(b2,b[,2:3],rep("NoRel",length(b2)))
  colnames(aa)=colnames(bb)=c("ref.otu","abund","path","status")
  
  OTU=rbind(aa,bb)
  OTU=OTU[order(OTU[,1]),]
  OTU[,1]=as.character(OTU[,1])
  OTU[,3]=as.character(OTU[,3])
  rownames(OTU)=c(seq(1,nrow(OTU)))
  return(OTU)
}