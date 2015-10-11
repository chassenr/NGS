uid <- read.table("uid.txt", sep = "\t", h = F)
colnames(uid) <- c("uid")
path <- read.table("blastpPath.txt", 
                   sep = "\t", 
                   h = F, 
                   stringsAsFactors = F, 
                   quote = "",
                   comment.char = ""
                   )
colnames(path) <- c("path", "uid")
ipr <- read.table("ipr.txt", sep = "\t", h = F)

uid$path <- c(NA)
for (i in unique(uid$uid)) {
  uid$path[uid$uid == i] <- path$path[path$uid == i]
}

uid$IPR <- ipr[1, 1]

write.table(uid, "uid_new.txt", row.names = F, col.names = F, quote = F, sep = "\t")
