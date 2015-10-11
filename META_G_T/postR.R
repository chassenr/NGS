if(!"plyr" %in% installed.packages()) {
 install.packages("plyr")
}
require(plyr)

###functional genes

#reading and formatting data
pfamacc2desc <- read.table(Sys.getenv("PFAMACC2DESC"), 
                           h = F,
                           sep = "\t",
                           quote = "",
                           comment.char = ""
                           )
colnames(pfamacc2desc) <- c("pfamacc", "pfamdesc")
rownames(pfamacc2desc) <- pfamacc2desc$pfamacc

pfam2go <- read.table("pfam2go.final",
                      h = F, 
                      sep = "\t",
                      quote = "",
                      comment.char = ""
                      )
colnames(pfam2go) <- c("pfamacc", "goterm")

goobo <- read.table("fullGO.NA", 
                    h = F,
                    sep = "\t", 
                    stringsAsFactors = F, 
                    fill = T, 
                    quote = "", 
                    comment.char = ""
                    )
colnames(goobo) <- c("goterm", "obsolete", "consider", "replaced", "alt_id", "go.name", "go.namespace")
rownames(goobo) <- goobo$goterm

SID <- strsplit(Sys.getenv("sample"), " ")[[1]]

pfam0 <- c()
fileName <- paste(SID, "_pfam_counts.txt", sep = "")
for (i in 1:length(SID)) {
  pfam0[[i]] <- read.table(fileName[i], sep = ",", h = F)
  colnames(pfam0[[i]]) <- c("pfamacc", fileName[i])
}
pfam <- join_all(pfam0, by = "pfamacc", type = "full", match = "all")
rownames(pfam) <- pfam$pfamacc
pfam_final <- pfam[, -which(colnames(pfam) == "pfamacc")]
pfam_final[is.na(pfam_final)] <- 0
pfam_final$pfamdesc <- pfamacc2desc[rownames(pfam_final), "pfamdesc"]

#creating GO tables
GO0 <- pfam2go[pfam2go$pfamacc %in% rownames(pfam_final), ]
GO_good <- GO0[GO0$goterm %in% goobo$goterm, ]
GO_good$goterm <- as.character(GO_good$goterm)
GO_good <- data.frame(GO_good, goobo[GO_good$goterm, 2:7])

#replacing obsolete terms
for (i in unique(GO_good$goterm)) {
  if(goobo[i, "replaced"] != "") {
    tmp <- strsplit(goobo[i, "replaced"], ",")[[1]]
    GO_good[GO_good$goterm == i, 2:8] <- goobo[tmp, 1:7]
  } else {
    if(goobo[i,"consider"]!=""){
      tmp <- strsplit(goobo[i, "consider"], ",")[[1]]
      GO_good[GO_good$goterm == i, 2:8] <- goobo[tmp, 1:7]
    }
  }
}


#adding GO terms to pfam table and checking alt_id if more than one goterm per namespace
pfam_go <- matrix(NA, nrow(pfam_final), 3)
rownames(pfam_go) <- rownames(pfam)
colnames(pfam_go) <- unique(GO_good$go.namespace)
for (i in 1:nrow(pfam_go)) { #for each pfamacc
  if (rownames(pfam_go)[i] %in% GO_good$pfamacc) { #if pfam has goterm
    for (j in 1:ncol(pfam_go)) { #for each namespace
      if (colnames(pfam_go)[j] %in% GO_good$go.namespace[GO_good$pfamacc == rownames(pfam_go)[i]]) { #if pfamacc has goterm in that namespace
        tmp <- GO_good$goterm[GO_good$pfamacc == rownames(pfam_go)[i] & GO_good$go.namespace == unique(GO_good$go.namespace)[j]] #vector with all goterms
        if (length(tmp) > 1) { #if more than 1 goterm
          for (k in 1:(length(tmp) - 1)) { #for each goterm
            if (any(grepl(tmp[k + 1], goobo[tmp[1:k], "alt_id"])) == TRUE) { #check if any goterms alt_ids of other goterms
              tmp[k + 1] <- NA #remove alt_ids
            }
          }
          pfam_go[i, j] <- paste(tmp[!is.na(tmp)], collapse = ",") #list goterms comma-separated
        } else {
          pfam_go[i, j] <- tmp #only one goterm for pfamacc
        }
      }
    }
  }
}

pfam_all <- data.frame(pfam_final, pfam_go)
write.table(pfam_all, "PfamCountsFinal.txt", sep = "\t", quote = F)
write.table(data.frame(GO_good), "GOdetail.txt", sep = "\t", quote = F)

#writing GO counts per namespace
GO_counts <- vector("list", length = 3)
names(GO_counts) <- colnames(pfam_go)
for (i in 1:3) {
  goterm <- unlist(strsplit(as.character(pfam_go[, i]), ","))
  tmp <- matrix(NA, length(unique(goterm[!is.na(goterm)])), length(SID))
  rownames(tmp) <- unique(goterm[!is.na(goterm)])
  colnames(tmp) <- paste(SID, "_go_counts.txt", sep = "")
  for (j in 1:nrow(tmp)) {
    tmp[j, ] <- colSums(pfam_all[grep(rownames(tmp)[j], pfam_go[, i]), 1:length(SID)])
  }
  GO_counts[[i]] <- data.frame(tmp, go.name = goobo[rownames(tmp), "go.name"])
}

write.table(GO_counts$molecular_function, "MF_counts.txt", sep = "\t", quote = F)
write.table(GO_counts$biological_process, "BP_counts.txt", sep = "\t", quote = F)
write.table(GO_counts$cellular_component, "CC_counts.txt", sep = "\t", quote = F)


#taxonomic composition
tax0 <- c()
fileName <- paste(SID, ".tax_slv", sep = "")
for (i in 1:length(SID)) {
  tmp <- read.table(fileName[i], sep = "\t", h = F, quote = "", comment.char = "")
  tax0[[i]] <- data.frame(table(droplevels(tmp[tmp != "Unclassified;", ])))
  colnames(tax0[[i]]) <- c("path", fileName[i])
}
tax <- join_all(tax0, by = "path", type = "full", match = "all")
rownames(tax) <- tax$path
tax_counts <- tax[, -which(colnames(tax) == "path")]
tax_counts[is.na(tax_counts)] <- 0

SILVA <- read.table(Sys.getenv("SILVA"), h = T, sep = "\t", quote = "", comment.char = "")
SILVAtaxopath <- function(tax, SILVA) {
  output <- matrix(NA, nrow = length(tax), ncol = length(levels(SILVA$rank)))
  colnames(output) <- levels(SILVA$rank) 
  
  for (i in 1:length(tax)) {
    for (j in 1:length(levels(SILVA$rank))) {
      if (paste(tax[[i]][1:j], collapse = ";") %in% SILVA$path) {
        output[i, as.character(SILVA[SILVA$path == paste(tax[[i]][1:j], collapse = ";"), "rank"])] <- as.character(SILVA[SILVA$path == paste(tax[[i]][1:j], collapse = ";"), "node"])
      }
    }
  }
  return(output)
}

tax_split <- strsplit(rownames(tax_counts), ";")
Tax <- SILVAtaxopath(tax_split, SILVA)
rownames(Tax) <- rownames(tax_counts)
Taxa <- Tax[,c("domain","major_clade","kingdom","phylum","class","order","family","genus")]
Taxa[Taxa == "uncultured"] <- "NA"
Taxpro <- Taxa[Taxa[, 1] != "Eukaryota", ]
for (i in 4:7) {
  if(sum(is.na(Taxpro[,i])) > 1) {
    test <- Taxpro[is.na(Taxpro[, i]), ]
    for (j in 1:nrow(test)) {
      if (sum(is.na(test[j, i:8])) == length(test[j, i:8])) {
        test[j, i] <- paste(test[j, (i - 1)], "_unclassified", sep = "")
        test[j, (i + 1):8] <- test[j, i]
      }
    }
    Taxpro[is.na(Taxpro[, i]), ] <- test
  }
  if(sum(is.na(Taxpro[, i])) == 1) {
    test <- Taxpro[is.na(Taxpro[, i]), ]
    if (sum(is.na(test[i:8])) == length(test[i:8])) {
      test[i] <- paste(test[(i - 1)], "_unclassified", sep = "")
      test[(i + 1):8] <- test[i]
    }
    Taxpro[is.na(Taxpro[, i]), ] <- test
  }
}
Taxpro[is.na(Taxpro[, 8]), 8] <- paste(Taxpro[is.na(Taxpro[, 8]), 7], "_unclassified", sep = "")

Taxeuk <- Taxa[Taxa[, 1] == "Eukaryota", ]
for (i in 2:7) {
  if (sum(is.na(Taxeuk[, i])) > 1) {
    test <- Taxeuk[is.na(Taxeuk[, i]), ]
    for (j in 1:nrow(test)) {
      if (sum(is.na(test[j, i:8])) == length(test[j, i:8])) {
        test[j, i] <- paste(test[j, (i - 1)], "_unclassified", sep = "")
        test[j, (i + 1):8] <- test[j, i]
      }
    }
    Taxeuk[is.na(Taxeuk[, i]), ] <- test
  }
  if (sum(is.na(Taxeuk[, i])) == 1) {
    test <- Taxeuk[is.na(Taxeuk[, i]), ]
    if (sum(is.na(test[i:8])) == length(test[i:8])) {
      test[i] <- paste(test[(i - 1)], "_unclassified", sep = "")
      test[(i + 1):8] <- test[i]
    }
    Taxeuk[is.na(Taxeuk[, i]), ] <- test
  }
}
Taxeuk[is.na(Taxeuk[, 8]), 8] <- paste(Taxeuk[is.na(Taxeuk[, 8]), 7], "_unclassified", sep = "")

tax_path <- rbind(Taxpro, Taxeuk)
tax_path <- tax_path[match(rownames(tax_counts), rownames(tax_path)), ]

write.table(data.frame(tax_counts, tax_path), "TaxonomyCountsFinal.txt", sep = "\t", quote = F)