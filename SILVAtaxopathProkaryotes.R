#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements: 
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# Function to allocate correct taxonomic rank to levels of taxonomic path (assigned using SILVA)
# It might be necessary to further curate the taxonomy afterwards (e.g. to account for NoRelative OTUs), but this function should do the majority of the work
# to avoid confusion, rownames of output should be named after OTU names ASAP
# input: 
# tax - list with length == number of OTUs and in each element an unsorted vector with taxonomic path of that OTU
#  e.g. [[1]]
#       [1] "Bacteria"          "Chlorobi"          "Ignavibacteria"    "Ignavibacteriales"
#
#       [[2]]
#       [1] "Bacteria"       "Chloroflexi"    "Ardenticatenia" "uncultured"    
#
#       [[3]]
#       [1] "Bacteria"         "Actinobacteria"   "Acidimicrobiia"   "Acidimicrobiales" "OM1 clade"       
#
# SILVA - reference table with columns path, node, rank
#  download: http://www.arb-silva.de/fileadmin/silva_databases/release_119/Exports/taxonomy/tax_slv_ssu_nr_119.txt (SILVA release 119)
#  Careful: this file has several errors!
#  or for prokaryotes (directly usable for below function): /trunk/clearingHouse/SILVA119taxmap_prok_curated.csv (read with: read.csv("SILVA119taxmap_prok_curated.csv",sep="\t",h=T))
#  e.g. path                                node                        rank
#       Archaea                             Archaea                     domain
#       Archaea;Ancient Archaeal Group(AAG) Ancient Archaeal Group(AAG) phylum
#       Archaea;Crenarchaeota               Crenarchaeota               phylum
# output:
# table with taxonomic path for each OTU, levels corresonding to correct taxonomic rank
#  e.g.    class           domain     family              genus        order               phylum          
#         "Ignavibacteria" "Bacteria" NA                  NA           "Ignavibacteriales" "Chlorobi"      
#         "Ardenticatenia" "Bacteria" NA                  NA           "uncultured"        "Chloroflexi"   
#         "Acidimicrobiia" "Bacteria" "OM1 clade"         NA           "Acidimicrobiales"  "Actinobacteria"# dependencies:
# to remove NA from the taxonomy table and replace them with unclassified (assuming there is no NA in domain column)
# manual curation still necessary afterwards
#Tax119b=Tax119a
#for(i in 2:5) {
#  if(sum(is.na(Tax119b[,i]))>1){
#    test=Tax119b[is.na(Tax119b[,i]),]
#    for(j in 1:nrow(test)){
#      if (sum(is.na(test[j,i:6]))==length(test[j,i:6])) {
#        test[j,i]=paste(test[j,(i-1)],"_unclassified",sep="")
#        test[j,(i+1):6]=test[j,i]
#      }
#    }
#    Tax119b[is.na(Tax119b[,i]),]=test
#  }
#  if(sum(is.na(Tax119b[,i]))==1){
#    test=Tax119b[is.na(Tax119b[,i]),]
#    if (sum(is.na(test[i:6]))==length(test[i:6])) {
#      test[i]=paste(test[(i-1)],"_unclassified",sep="")
#      test[(i+1):6]=test[i]
#    }
#    Tax119b[is.na(Tax119b[,i]),]=test
#  }
#}
#Tax119b[is.na(Tax119b[,6]),6]=paste(Tax119b[is.na(Tax119b[,6]),5],"_unclassified",sep="")
# dependencies:
# /trunk/clearingHouse/SILVA119taxmap_prok_curated.csv (or any other taxmap)
#=============================================================================
#documentation end

SILVAtaxopath <- function(tax, SILVA) {
  # create matrix with nrow = number of OTUs and ncol = number of possible taxonomic levels
  output <- matrix(NA, nrow = length(tax), ncol = length(levels(SILVA$rank))) 
  colnames(output) <- levels(SILVA$rank)
  
  for (i in 1:length(tax)) {
    if (length(tax[[i]]) == 6) {
      output[i, c("domain", "phylum", "class", "order", "family", "genus")] <- tax[[i]]
    } else {
      # this loop goes through the taxonomic path of an OTU one level at a time and assigns the name of each level to its correct rank according to SILVA
      for (j in 1:length(levels(SILVA$rank))) {
        if (paste(tax[[i]][1:j], collapse = ";") %in% SILVA$path) {
          output[i, as.character(SILVA[SILVA$path == paste(tax[[i]][1:j], collapse = ";"), "rank"])] <- as.character(SILVA[SILVA$path == paste(tax[[i]][1:j], collapse = ";"), "node"])
        }
      }
    }
  }
  return(output)
}

