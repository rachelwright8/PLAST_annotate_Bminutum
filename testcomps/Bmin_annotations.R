setwd("~/Dropbox/genomes/Bminutum_JParkinson_Plast/testcomps/")
library(tidyverse)

blast_gene <- read.delim("blast_iso2gene.tab", header=F)
head(blast_gene)

# change the names to match the new annotations (isogroup instead of comp)
blast_gene$V1 <- gsub("comp", "isogroup", blast_gene$V1)
head(blast_gene)

# get rid of columns that don't have annotations
blast_gene <- blast_gene[!blast_gene$V2=="-",]
head(blast_gene)

plast_gene <- read.delim("plast_iso2gene.tab", header=F)
head(plast_gene)
# get rid of second column (not in BLAST)
plast_gene <- plast_gene %>% select(V1,V3)
head(plast_gene)

nrow(blast_gene)-nrow(plast_gene) # 5759 more genes annotated using BLAST than PLAST

table(blast_gene$V1 %in% plast_gene$V1) # 6921 in blast that aren't in plast
table(plast_gene$V1 %in% blast_gene$V1) # 1162 in plast that aren't in blast

matches <- merge(plast_gene,blast_gene,by=1)
names(matches) <- c("isogroup", "plastGene", "blastGene")
head(matches)

# see if the strings match
# get rid of characters besides the gene description
matches$plastClean <- gsub("OS.*","",matches$plastGene)
matches$plastClean <- gsub("^ *", "", matches$plastClean)
matches$plastClean <- gsub("* $", "", matches$plastClean)
matches$blastClean <- gsub(".*;", "", matches$blastGene)
head(matches)

# how many match?
table(matches$plastClean==matches$blastClean)
# FALSE  TRUE 
# 6573 10739

# which ones don't match?
nomatch <- matches[!matches$plastClean==matches$blastClean,][c(4,5)]
tail(nomatch)
# Checking by eye, it looks lie the genes are usually only slightly different
# Example:
# 60S ribosomal protein L3 vs 60S ribosomal protein L3-2

# summarize genes and plot
plast_results <- data.frame(matching=table(plast_gene$V1 %in% blast_gene$V1)[2],
                            unique=table(plast_gene$V1 %in% blast_gene$V1)[1])
blast_results <- data.frame(matching=table(blast_gene$V1 %in% plast_gene$V1)[2],
                            unique=table(blast_gene$V1 %in% plast_gene$V1)[1])

gene_results_numbers <- as.data.frame(rbind(blast_results,plast_results)) %>% 
  mutate(program=c("blast","plast")) %>%
  gather(type, number, 1:2) %>%
  ggplot() + geom_bar(aes(y = number, x = program, fill = factor(type,levels=c("unique", "matching"))), stat="identity") + 
  guides(fill=guide_legend(title="hit type")) +
  theme_bw()
gene_results_numbers


# In other cases, it looks like the plast has a description where the blast has "predicted protein" or similar. 
# This makes sense because the BLAST annotations are older. The uniprot database has been updated since these annotation files were made.
nomatch %>% filter(grepl("Predicted",blastClean)) %>% nrow() # 109 instances of "predicted protein" for blast
nomatch %>% filter(grepl("Predicted",plastClean)) %>% nrow() # 1 instance of "predicted protein" for plast
nomatch %>% filter(grepl("Predicted",blastClean) & grepl("Predicted",plastClean)) %>% nrow() # every "Predicted protein" in blast has an annotation in plast

nomatch %>% filter(grepl("Uncharacterized",blastClean)) %>% nrow() # 626 instances of "predicted protein" for blast
nomatch %>% filter(grepl("Uncharacterized",plastClean)) %>% nrow() # 160 instance of "predicted protein" for plast
nomatch %>% filter(grepl("Uncharacterized",blastClean) & grepl("Uncharacterized",plastClean)) %>% nrow() # 65 isogroups are uncharacterized in both

# summarize and plot

plast_results_annot <- data.frame( predicted=nomatch %>% filter(grepl("Predicted",plastClean)) %>% nrow(),
                             uncharacterized=nomatch %>% filter(grepl("Uncharacterized",plastClean)) %>% nrow())

blast_results_annot <- data.frame(predicted=nomatch %>% filter(grepl("Predicted",blastClean)) %>% nrow(),
                            uncharacterized=nomatch %>% filter(grepl("Uncharacterized",blastClean)) %>% nrow())
gene_results_annot <- as.data.frame(rbind(blast_results_annot,plast_results_annot)) %>% 
  mutate(program=c("blast","plast")) %>%
  gather(type, number, 1:2) %>%
  ggplot() + geom_bar(aes(y = number, x = program, fill = factor(type,levels=c("predicted", "uncharacterized"))), stat="identity") + 
  guides(fill=guide_legend(title="hit type")) +
  theme_bw()
gene_results_annot

##################################################################################################################################
##################################################################################################################################
# GO annotations -----
##################################################################################################################################
##################################################################################################################################

blast_go <- read.delim("blast_iso2go.tab", header=F)
head(blast_go)

# change the names to match the new annotations (isogroup instead of comp)
blast_go$V1 <- gsub("comp", "isogroup", blast_go$V1)
head(blast_go )

# get rid of isogroups that don't have annotations
blast_go  <- blast_go [!blast_go $V2=="-",]
head(blast_go )

# load plast annotations for GO
plast_go <- read.delim("../testgo.tab", header=F)
head(plast_go)

# get rid of second column (short gene ID) that isn't in the blast file
plast_go <- plast_go %>% select(V1,V3) %>% rename(V2="V3")
head(plast_go)

# get rid of isogroups that don't have annotations
plast_go  <- plast_go[!plast_go $V2=="noMatch",]
head(plast_go)

nrow(blast_go)-nrow(plast_go) # 1024 more GO annotations in BLAST

table(blast_go$V1 %in% plast_go$V1) # 3081 in blast that aren't in plast
table(plast_go$V1 %in% blast_go$V1) # 2057 in plast that aren't in blast

gomatches <- merge(plast_go,blast_go,by=1)
names(gomatches) <- c("isogroup", "plastGO", "blastGO")
head(gomatches)

# summarize and plot
blast_results_go <- data.frame(matching=table(blast_go$V1 %in% plast_go$V1)[2],
                            unique=table(blast_go$V1 %in% plast_go$V1)[1])
plast_results_go <- data.frame(matching=table(plast_go$V1 %in% blast_go$V1)[2],
                            unique=table(plast_go$V1 %in% blast_go$V1)[1])
go_results_numbers <- as.data.frame(rbind(blast_results_go,plast_results_go)) %>% 
  mutate(program=c("blast","plast")) %>%
  gather(type, number, 1:2) %>%
  ggplot() + geom_bar(aes(y = number, x = program, fill = factor(type,levels=c("unique", "matching"))), stat="identity") + 
  guides(fill=guide_legend(title="hit type")) +
  theme_bw()
go_results_numbers

# Compare GO annotations
plast_go_long <- strsplit(as.character(plast_go$V2), split = ";")
plast_go_long <- data.frame(V1 = rep(plast_go$V1, sapply(plast_go_long, length)), 
                            V2 = unlist(plast_go_long))
head(plast_go_long)
length(unique(plast_go_long$V1)) # 17934 isogroups with annotations

blast_go_long <- strsplit(as.character(blast_go$V2), split = ";")
blast_go_long <- data.frame(V1 = rep(blast_go$V1, sapply(blast_go_long, length)), 
                            V2 = unlist(blast_go_long))
head(blast_go_long)
length(unique(blast_go_long$V1)) # 18958 isogroups with annotations

# matches?
gomatchlong <- merge(plast_go_long, blast_go_long, by=c("V1","V2"))
nrow(gomatchlong) # 20347 total instances where the isogroup AND GO term matched
length(unique(gomatchlong$V1)) # 9964 isogroups have at least one matching GO term

