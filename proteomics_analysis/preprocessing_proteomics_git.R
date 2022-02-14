#Data manipulation packs
library(readxl)
library(reshape2)
library(robustbase)
library(stringr)
library(boot)
library(splitstackshape)
library(cometr)

#Bioconductor packs
library(Biobase)
library(convert)

#Set directory
input_path <- "/Users/Derek/Downloads/proteins_with_gene_names/S2BD1/"
output_path <- "/Users/Derek/Downloads/proteins_with_gene_names/"

#Recursively iterate over "protein.csv" files#
files <- list.files(input_path, pattern="protein[[:alnum:]]+.csv", recursive=TRUE, full.names=TRUE)

#Look for proteins that are concurrent across all subjects
#and save their accession, gene names, and synonyms.
codex <- data.frame()
counts <- c()
access <- NULL
les_noms <- colnames(df)

for (file in files) {
  #For keeping the segment and blood draw info
  df <- read.csv(file,header = T)
  df <- df[!is.na(df$Area_Sample), ]
  df <- df[df$Area_Sample > 0, ]
  colnames(df) <- les_noms
  counts <- c(counts, nrow(df))
  codex <- rbind(codex,df)
  if(is.null(access)) {
    access <- df$Accession
    gene <- df$Gene_Name
    syn <- df$Gene_Synonyms
  }
  else {
    access <- as.character(df[df$Accession%in%access,"Accession"])
    gene <- as.character(df[df$Accession%in%access, "Gene_Name"])
    syn <- as.character(df[df$Accession%in%access, "Gene_Synonyms"])
  }
}

#Code snippet for checking the protein coverage per blood draw across subjects;
#it specifically looks for proteins present across at least 50% of the subjects.
num_subjects1 <- length(files)
proteins1 <- unique(codex$Accession)
length(proteins1)
tbd1 <- sort(table(codex$Accession),decreasing=T)
sum(tbd1>=ceiling(0.5*num_subjects))
counts1 <- counts

(num_subjects2 <- length(files))
proteins2 <- unique(codex$Accession)
(length(proteins2))
tbd2 <- sort(table(codex$Accession),decreasing=T)
sum(tbd2>=ceiling(0.5*num_subjects2))
counts2 <- counts

(num_subjects3 <- length(files))
proteins3 <- unique(codex$Accession)
(length(proteins3))
tbd3 <- sort(table(codex$Accession),decreasing=T)
sum(tbd3>=ceiling(0.5*num_subjects3))
counts3 <- counts

sum(names(tbd1[tbd1>=ceiling(0.5*num_subjects1)])[names(tbd1[tbd1>=ceiling(0.5*num_subjects1)])%in%names(tbd2[tbd2>=ceiling(0.5*num_subjects2)])]%in%
  names(tbd3[tbd3>=ceiling(0.5*num_subjects3)]))

prot_list <- names(tbd1[tbd1>=ceiling(0.5*num_subjects1)])[names(tbd1[tbd1>=ceiling(0.5*num_subjects1)])%in%names(tbd2[tbd2>=ceiling(0.5*num_subjects2)])][
names(tbd1[tbd1>=ceiling(0.5*num_subjects1)])[names(tbd1[tbd1>=ceiling(0.5*num_subjects1)])%in%names(tbd2[tbd2>=ceiling(0.5*num_subjects2)])]%in%
  names(tbd3[tbd3>=ceiling(0.5*num_subjects3)])]

write.csv(prot_list, "access.csv")

#Make a quick visualization of the protein coverage levels.
barplot(tbd3, col="white")
abline(h=ceiling(0.5*num_subjects3), lty = 2, col = 'red')
title(main="Distribution of All Proteins in Segment 2 BD3")
legend(x="topright",legend = "Half of subjects",lty=2,col='red')

#Create a matrix of protein abundance figures with rows corresponding to subjects
prot <- c()
subjects <- c()
for (file in files) {
  loc <- gregexpr("proteins.csv",file)[[1]]
  file_name <- paste(substr(file, loc-18, loc-2), "proteins.csv", sep="_")
  loc_ <- gregexpr("_", file_name)[[1]][1]
  subjects <- c(subjects, paste("GEWP",substr(file_name,  loc_ + 1, loc_ + 3),sep=""))
  df <- read.csv(file, header = T)
  prot <- rbind(prot, as.numeric(df[df$Accession%in%access,'Area_Sample']))
}

#Creating a collapsed column name vector so that each feature name is either
#the gene name or, if the gene name is not available, the accession id.
gene[is.na(gene)] <- ""
syn[is.na(syn)] <- ""
ft_names <- gene
idx_empty <- sapply(ft_names, function(x){x == ""})
ft_names[idx_empty] <- access[idx_empty]
ft_names <- sapply(ft_names,
                   function(x){
                     idx <- gregexpr(";",x)[[1]][1]
                     if (idx != -1) return(substr(x,1,idx-1))
                     else return(x)
                     }
                   )

colnames(prot) <- ft_names
rownames(prot) <- subjects

#Try to batch normalize so that the protein data has a median of one per col
prot_med <- apply(prot, 2, median, na.rm = T)
prot_norm  <- sweep(prot, 2, prot_med, '/')

#Then take natural log of the protein data
log_prot <- log(prot_norm)

#Proteins df without NAs
log_prot_sans_na <- log_prot[,colSums(is.na(log_prot))==0]

#Replace -Inf with otherwise lowest real number as with standard batch normalization
log_prot_sans_na <- replace(log_prot_sans_na, log_prot_sans_na == -Inf, NA)
log_prot_mins <-  apply(log_prot_sans_na,2,min,na.rm=T)
idx_na <- is.na(log_prot_sans_na)
log_prot_sans_na[idx_na] <- log_prot_mins[col(log_prot_sans_na)][idx_na]

write.csv(log_prot_sans_na, file = paste(output_path,"log_prot_S2BD3.csv",sep=""))

log_prot_sans_na1 <- read.csv('log_prot_S2BD1.csv',header=T)[,-1]
log_prot_sans_na2 <- read.csv('log_prot_S2BD2.csv',header=T)[,-1]
log_prot_sans_na3 <- read.csv('log_prot_S2BD3.csv',header=T)[,-1]

#Check which proteins are common across the blood draws
c1 <- colnames(log_prot_sans_na1)
c2 <- colnames(log_prot_sans_na2)
c3 <- colnames(log_prot_sans_na)

c1_in_c2 <- c1[c1%in%c2]
in_all <-c1_in_c2[c1_in_c2%in%c3]
length(in_all)

log_prot1 <- log_prot
log_prot2 <- log_prot
log_prot3 <- log_prot

c1a <- colnames(log_prot1)
c2a <- colnames(log_prot2)
c3a <- colnames(log_prot3)

length(c1a[c1a[c1a%in%c2a]%in%c3a])

save.image(paste("/Volumes/GoogleDrive-115111199924997198421/My Drive/Objectives/Code/",
                 "processing_protein_v0.1",sep=""))
