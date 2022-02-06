#Data manipulation packs
library(readxl)
library(reshape2)
library(stringr)
library(splitstackshape)
library(gtools)

#Bioconductor pack
library(Biobase)
library(convert)
library(impute)

#Set directory.
input_path <- "/Users/Derek/Downloads/proteins_with_gene_names/S2BD1/"
output_path <- "/Volumes/GoogleDrive-115111199924997198421/My Drive/Objectives/
Code/Protein_Imputation/"

#Recursively iterate over "protein.csv" files#
files <- list.files(input_path, pattern="protein[[:alnum:]]+.csv", recursive=TRUE, full.names=TRUE)

#Look at the list of proteins that are common across at least 50% of the subjects.
access <- read.csv('/Users/Derek/Downloads/access.csv')

prot_list <- Reduce(intersect,list(names(tbd1[tbd1>=ceiling(max(tbd1)/2)]), 
                      names(tbd2[tbd2>=ceiling(max(tbd2)/2)]), 
                      names(tbd3[tbd3>=ceiling(max(tbd3)/2)])))

write.csv(prot_list,"check.csv")

names(tbd1[tbd1>=ceiling(max(tbd1)/2)])
names(tbd2[tbd2>=ceiling(max(tbd2)/2)])
names(tbd3[tbd3>=ceiling(max(tbd3)/2)])

#Create a matrix of protein abundance figures with rows corresponding to subjects.
prot <- data.frame(matrix(ncol=length(access[,1]),nrow=0))
colnames(prot) <- access[,1]
subjects <- c()
protl <- c()
for (file in files) {
  loc <- gregexpr("proteins.csv",file)[[1]]
  file_name <- paste(substr(file, loc-18, loc-2), "proteins.csv", sep="_")
  loc_ <- gregexpr("_", file_name)[[1]][1]
  subjects <- c(subjects, paste("GEWP",substr(file_name,  loc_ + 1, loc_ + 3),sep=""))
  df <- read.csv(file, header = T)
  df <- df[,c("Accession", "Area_Sample")]
  temp <- data.frame(matrix(df[,2],nrow=1))
  colnames(temp) <- df[,1]
  prot <- smartbind(prot, temp)
}

#KNN impute the missing values.
prot <- prot[,access[,1]]
rownames(prot) <- subjects
imputation <- impute.knn(t(as.matrix(prot)),k=5,rng.seed=23)
prot_imputed <- data.frame(t(imputation$data))
rownames(prot_imputed) <- rownames(prot)
colnames(prot_imputed) <- access[,2]

#First, save the knn imputed version as csv.
setwd(output_path)
write.csv(prot_imputed, "s2bd1_prot_imputed_raw.csv")

#Then divide each column by its median.
prot_imputed <- apply(prot_imputed,2,
                      function(x){
                        x/median(x)
                      })
#Set 0 abundance values to the lowest non-zero values.
prot_imputed <- apply(prot_imputed,2,
      function(x) {
        if (min(x)==0){
          x[x==0] = sort(unique(x),partial=2)[2]
          x
        } else x
      })

prot_imputed <- log(prot_imputed)
write.csv(prot_imputed, "s2bd1_log_prot_imputed.csv")
