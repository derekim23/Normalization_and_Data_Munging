#Data manipulation packs
library(readxl)
library(reshape2)
library(robustbase)
library(stringr)
library(plyr)

#Bioconductor packs
library(Biobase)
library(convert)

uniprot_mapping <- function(ids) {
  #Map UniProt IDs to gene names
  tryCatch(
  { 
    uri <- 'http://www.uniprot.org/uniprot/?query='
    idStr <- paste(ids, collapse="+or+")
    format <- '&format=tab'
    fullUri <- paste0(uri,idStr,format)
    dat <- read.delim(fullUri)
    dat <- dat[dat$Entry%in%ids,]
    if (length(dat$Gene.names)>1){
      gene_ids <- sapply(dat$Gene.names, extract_part)
    } else{
      gene_ids <- extract_part(dat$Gene.names)
    }
    if (length(dat$Gene.names)>1){
      remaining_gene_ids <- sapply(dat$Gene.names, extract_part, first = F)
    } else{
      remaining_gene_ids <- extract_part(dat$Gene.names, first = F)
    }
    dat$Gene.names <- gene_ids
    colnames(dat)[grep("Gene.names", colnames(dat))] <- "Gene_Name"
    dat <- cbind(dat,remaining_gene_ids)
    colnames(dat)[length(dat)] <- "Gene_Synonyms"
    dat
  },
  error = {
    function(e) NULL
  })
}

extract_part <- function(x, first = T){
  split_point <- gregexpr(" ",x)[[1]][1]
  if (split_point == -1){
    if(first)
      return(substr(x,1,nchar(x)))
    else
      return("")
  } else{
    if (first)
      return(substr(x,1,split_point - 1))
    else{
      return(substr(x,split_point + 1,nchar(x)))  
    }
  }
}

extract_id <- function(x, pattern = "\\|"){
#Extract the substring that corresponds to UniProt ID
  idx_matches <-gregexpr(x, pattern = pattern)[[1]]
  if (length(idx_matches)>1)
    return(substr(x,idx_matches[1]+1,idx_matches[2]-1))
  else
    return(substr(x,1,idx_matches[1]-1))
}

match_genes <- function(file_name){
  #Extract UniProt IDs from the "Accession" column and then merge them into the main DF ("file" DF)
  file <- read.csv(file_name, header=TRUE)
  colnames(file)[7] <- "Area_Sample"
  uni_prot <- sapply(file$Accession,extract_id,pattern="\\|")
  uni_prot <- cbind(names(uni_prot),uni_prot)
  colnames(uni_prot) <- c("Accession", "UniProt_ID")
  file <- merge(file, uni_prot, by = "Accession", all.x = T)

  #Create a mapping between UniProt IDs and gene names
  #Break the process down into chunks to prevent errors
  ids <- na.omit(file$UniProt_ID)
  first <- floor(length(ids)/3)
  second <- floor(length(ids)/3*2)
  third <- length(ids)
  ids2genes1 <- uniprot_mapping(ids[1:first]) #as.data.frame(t(sapply(file$UniProt_ID,uniprot_mapping)))
  ids2genes2 <- uniprot_mapping(ids[(first+1):second])
  ids2genes3 <- uniprot_mapping(ids[(second+1):third])
  ids2genes <- rbind(ids2genes1,ids2genes2,ids2genes3)

  if (is.null(ids2genes1)||is.null(ids2genes2)||is.null(ids2genes3)){
    print(paste("Retrieval error with file: ", file_name, sep=""))
    return(file_name)
  }
  
  #Merge gene names into the main DF
  file <- merge(file, ids2genes[,c("Entry", "Gene_Name", "Gene_Synonyms")], by.x="UniProt_ID", by.y = "Entry", all.x = T)
  #colnames(file)[dim(file)[2]] <- "Gene_Name"
  file <- file %>% relocate(c(UniProt_ID, Gene_Name, Gene_Synonyms, Area_Sample), .after = Accession)
  return(file)
}

#Set directory
input_path <- "/Users/Derek/Downloads/Results Analyzed/S2BD1 Analyzed/"
output_path <- "/Users/Derek/Downloads/proteomics/S2BD1/"
err_path <- "/Users/Derek/Downloads/didnt_work/"

#Recursively iterate over "protein.csv" files#
files <- list.files(input_path, pattern="protein[[:alnum:]]+.csv", recursive=TRUE, full.names=TRUE)
#files <- read.csv(paste(err_path,'err_files.csv',sep=''))[,2]
#files <- as.character(files)
#Create an empty vector for tracking error instances wherein the program fails to download from UniProt
err_vec <- c()
for (file in files) {
  #For keeping the segment and blood draw info
  loc <- gregexpr("proteins.csv",file)[[1]]
  file_name <- paste(substr(file, loc-18, loc-2), "proteins.csv", sep="_")
  #Read in file
  df <- match_genes(file)
  
  #Track error instances
  if(is.character(df))
    err_vec <- append(err_vec,df)
  else
    #Otherwise, write the output
    write.csv(df,paste(output_path,file_name,sep=""))
}

#Write the error instances to csv
if(length(err_vec) > 0){
  err_df <- as.data.frame(err_vec)
  colnames(err_df) <- "files_that_failed"
  write.csv(err_df, paste(err_path,"err_files.csv",sep=""))
}

#Now, check if the original and processed dataframes have the same number of rows
output <- list.files(output_path, pattern="protein[[:alnum:]]+.csv", recursive=TRUE, full.names=TRUE)
input <-  list.files(input_path, pattern="protein[[:alnum:]]+.csv", recursive=TRUE, full.names=TRUE)

check <- data.frame(cbind(input,output))

d <- data.frame(sort(input), sort(output))
mismatch_vec <- c()
for (i in 1:nrow(d)){
  df1 <- read.csv(d[i,1])
  df2 <- read.csv(d[i,2])
  if (nrow(df1)!=nrow(df2)){
    mismatch_vec <- c(mismatch_vec,d[i,1])
    print(paste(nrow(df1)-nrow(df2)))
  }
}
length(mismatch_vec)

#Q short script for correcting a column name
for (f in output){
  o <- read.csv(f,header=T)[,-c(1)]
  colnames(o)[4] <- "Gene_Synonyms"
  #colnames(o)[3] <- "Gene_Name"
  write.csv(o,f)
}
