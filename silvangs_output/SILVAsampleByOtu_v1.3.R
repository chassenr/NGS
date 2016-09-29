#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements: 
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# creates Sample by OTU table from SILVA output, also changes taxonomic path to NoRelative where applicable
# in case of inconsistencies in taxonomic path of a references OTU (classified differently in separate samples), 
#  resolves these by pruning the path to the best resolved common taxonomic level (excluding NoRelative)
# for an example of complete SILVA pipeline in R: see \trunk\clearingHouse\SILVApipe.R
# input: 
# *otu_mapping.stats and samples.stats files from SILVA output (under SILVA\results\ssu\stats\data)
# read with:
#  Map=read.csv("Project---ssu---otu_mapping.stats",sep="\t",h=T)
#  Sample=read.csv("Project---ssu---samples.stats",sep="\t",h=T)
# output of readSILVA from all sample in 1 list read in SAME ORDER as in *samples.stats:
#  e.g. x=list(Sample1,Sample2,Sample3,...,Samplen)
# output:
# list with sample by OTU table and taxonomic path by OTU table
# dependencies:
# uses output from readSILVA
#=============================================================================
#documentation end

sampleByOTU=function(x, Map, Sample) {
  found.in=strsplit(as.character(Map$Found.In),split=",", fixed=T)
  n.found.in=sapply(found.in, function(x) {length(x)})
  
  multiple.occurences=Map[n.found.in>1,]
  members=strsplit(as.character(multiple.occurences$Members),split=",")
  
  x1=x
  
  for (i in 1:length(x1)){
    x1[[i]][x1[[i]][,4]=="NoRel",3]=c("NoRelative")
    
    for (j in 1:nrow(multiple.occurences)) {
      x1[[i]][,1][x1[[i]][,1] %in% members[[j]]]=as.character(multiple.occurences[j,3])
    } 
  }
  
  names_all=data.frame(c(NA),c(NA))
  for(i in 1:length(x1)){
    if(i==1){names_all=x1[[1]][,c(1,3)]}
    if(i!=1){names_all=rbind(names_all, x1[[i]][,c(1,3)])}
  }
  
  NamesOTU=unique(names_all)
  NamesOTU=NamesOTU[order(NamesOTU$ref.otu),]
  
  if (!length(NamesOTU$ref.otu)==length(unique(NamesOTU$ref.otu))) {
    duplicate_tax=NamesOTU[NamesOTU[,1] %in% names(which(table(NamesOTU[,1])!=1)),]
    duplicate_tax1=duplicate_tax[!duplicate_tax$path=="NoRelative",]
    
    for (i in 1:length(unique(duplicate_tax1$ref.otu))) {
      temp=duplicate_tax1$path[duplicate_tax1$ref.otu==unique(duplicate_tax1$ref.otu)[i]]
      if (length(temp)==1) {
        NamesOTU$path[NamesOTU$ref.otu==unique(duplicate_tax1$ref.otu)[i]]=temp
      }
      if (!length(temp)==1) {
        common=c()
        for (j in 1:min(sapply(strsplit(temp,";"), function(x) length(x)))){
          common[j]=length(unique(sapply(temp, function(x) strsplit(x, ";")[[1]][j])))
        }
        NamesOTU$path[NamesOTU$ref.otu==unique(duplicate_tax1$ref.otu)[i]]=
          paste(c(strsplit(temp,";")[[1]][common==1],""), collapse=";")
      }
    }
    NamesOTU=unique(NamesOTU)
    NamesOTU=NamesOTU[order(NamesOTU$ref.otu),]
  } else {duplicate_tax=c()}
  
  SampleOTU=matrix(NA,nrow(NamesOTU),length(x1))
  colnames(SampleOTU)=as.character(Sample[,2])
  rownames(SampleOTU)=NamesOTU$ref.otu
  for (i in 1:length(x1)){
    temp <- aggregate(x1[[i]][,2], list(x1[[i]][,1]), sum)
    ref.otu <- temp$Group.1
    SampleOTU[ref.otu,i]=temp[,2]
  }
  
  SampleOTU=SampleOTU[order(rownames(SampleOTU)),]
  SampleOTU[is.na(SampleOTU)]<-0
  
  result=list(NamesOTU,SampleOTU,duplicate_tax)
  return(result)
}

