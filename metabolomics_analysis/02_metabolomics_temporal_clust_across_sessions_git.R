#This version has no imputations - throws out data to keep things logitudinal
#Data manipulation packs
library(readxl)
library(reshape2)
library(robustbase)
library(stringr)
library(plyr)

#Bioconductor packs
library(Mfuzz)
library(Biobase)
library(convert)
library(tscR)

#SignedTest packs
library(psych)
library(BSDA)
library(DescTools)

library(Cluster)

getmode <- function(x) {
  #For finding the mode (works even when multimodal)
  unique_x <- unique(x)
  tabulate_x <- tabulate(match(x, unique_x))
  unique_x[tabulate_x == max(tabulate_x)]
}

unstack <- function(x, num_analytes=1269){
  #For unstacking stacked columns by blood draws.
  num_draws <- length(x)/num_analytes
  agg <- NULL
  for(i in 1:num_draws) { 
    sub <- x[(1:num_analytes)+(i-1)*num_analytes]
    agg <- cbind(agg,sub)
  }
  return(agg)
}

is.nan.data.frame <- function(x){
  #Because "is.nan" doesn't work with dataframes like "is.na," have to write this function. 
  do.call(cbind, lapply(x, is.nan))
}

standardize <- function(df){
  #Scale dataframes by column
  df <- apply(df,2,scale)
  return(df)
}

bic <- function(fit){
  #For calculating the BIC criterion for Clustering
  #source: https://stats.stackexchange.com/questions/90769/using-bic-to-estimate-the-number-of-k-in-kmeans/251169#251169
  #number of Clusters
  m = nrow(fit$centers)
  # size of the Clusters
  n = fit$size
  #size of data set
  N = sum(fit$size)
  d = ncol(fit$centers)
  
  #compute variance for all Clusters beforehand
  cl_var = (1.0 / (N - m) / d) * fit$withinerror
  BIC = - 0.5 * m * log(N) * (d+1)

  for (i in 1:m){
    BIC = BIC + sum(n[i] * log(n[i]+1e-12) - n[i] * log(N) -
                       ((n[i] * d) / 2) * log(2*pi*cl_var) -
                       ((n[i] - 1) * d/ 2))
  }
  return(-BIC)
}

longitudify <- function(df){
  #In cases where the number of draws is inconsistent across blood draw sessions, adjust for attrition to keep
  #the data longitudinal across blood draws (or even segments)
  #Assumes log2fc or some abundance data is already merged with metadata (i.e. subject ids, subpathways, superpathways).
  
  segments <- unique(df$CUSTOM_ATTRIBUTE_1)
  len_segments <- length(segments)
  lim_idx <- grep(F, sapply(df[1,], is.numeric))[c(1,2)]
  lim_idx[1] <- lim_idx[1] + 1
  lim_idx[2] <- lim_idx[2] - 1
  
  #Take the intersection of subject ids present across all sessions.
  cap <- Reduce(intersect, list(df[df$CUSTOM_ATTRIBUTE_1==segments[1],]$CLIENT_SAMPLE_ID,
                                df[df$CUSTOM_ATTRIBUTE_1==segments[2],]$CLIENT_SAMPLE_ID,
                                df[df$CUSTOM_ATTRIBUTE_1==segments[3],]$CLIENT_SAMPLE_ID))
  
  #Then filter the dataframe.
  df <- df[sapply(df$CLIENT_SAMPLE_ID, is.element, set = cap),]
  return(df[!duplicated(df[,c("CLIENT_SAMPLE_ID","CUSTOM_ATTRIBUTE_1")]),])
}

take_substr <- function(x){
  #Take substrings of rownames to extract molecule names
  #when R changes duplicate row names by adding prefixes and suffixes
  loc <- str_locate(x,'\\.')[1]
  if (is.na(loc))
    return(substr(x,2,nchar(x)))
  else
    return(substr(x,2,loc-1))
}

unique_mode <- function(vec){
  #Tells whether a given mode is unique for a given row
  #Used for dealing with multiple Cluster assignments for a given analyte.
  sum(is.na(vec)) == (length(vec) - 1)
}

hypergeo_test <- function(universe, select_m, alpha = 0.05){
  #Hypergeometric hypothesis testing
  #Right tailed alternative for checking the abundance of 
  #certain subpath members within a Cluster vs. the rest of the dataset.
  select_sub <- unique(select_m$SUB_PATHWAY)
  p_vals <- matrix(NA,nrow(select_m[,c("SUB_PATHWAY", "Cluster")]),8) #dataframes are incredibly slow; use matrices
  k <- 1
  N <- nrow(universe)
  for (i in 1:length(select_sub)){
    sub_patt <- unique(select_m[select_m$SUB_PATHWAY==select_sub[i],]$Cluster)
    for (j in 1:length(sub_patt)){
      select_metabolytes <- 
        select_m[
          select_m[,'Cluster']==sub_patt[j] & select_m[,'SUB_PATHWAY'] == select_sub[i],
        ]
      x <- sum(universe[,'Cluster']==sub_patt[j] &
                 universe[,'SUB_PATHWAY'] == select_sub[i] & universe$CHEM_ID%in%select_metabolytes$CHEM_ID)
      m <- sum(universe[,"SUB_PATHWAY"]==select_sub[i])
      p <- phyper(q = x-1, #p(Q > q) when using upper tail
                  m = m, 
                  n = N - m, 
                  k = sum(universe[,'Cluster']==sub_patt[j]), 
                  lower.tail = FALSE, log.p = FALSE
      )
      subgroup <- unique(select_metabolytes[,c("CHEM_ID", "HMDB", "CHEMICAL_NAME", "SUPER_PATHWAY")])
      for (r in 1:nrow(subgroup)){
        p_vals[k,] <- c(subgroup[r,]$CHEM_ID, subgroup[r,]$HMDB, subgroup[r,]$CHEMICAL_NAME, 
                        subgroup[r,]$SUPER_PATHWAY, 
                        select_sub[i], sub_patt[j], p, NA)
        k = k + 1
      }
    }
  }
  p_vals <- unique(p_vals)
  p_vals[,8] <- p.adjust(as.numeric(p_vals[,7]),method="hochberg") #p-val adjustment for multiple testing.
  p_vals <- p_vals[order(as.numeric(p_vals[,8])),]
  Significant <- as.numeric(p_vals[,8]) < alpha
  p_vals <- cbind(p_vals,Significant)
  #p_vals <- p_vals[p_vals[,9] == T,]
  colnames(p_vals)[1:8] <- c("CHEM_ID", "HMDB", "CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY", "Cluster", 
                             "p_value", "Hochberg_p_value")
  return(p_vals)
}

slope <- function(x) {
  #Calculate the slopes for consecutive column pairs per row.
  time <- 1:ncol(x)
  N <- dim(x)[1]
  M <- dim(x)[2]
  m1 <- x[,-c(M)]
  m2 <- x[,-c(1)]
  mdif <- m2 - m1
  slopes <- mdif %*% diag(1 / (time[2:M] - time[seq_len(M - 1)]))
  return(slopes)
}

calc_dx_across <- function(df){
  #Calculate the delta between consecutive time steps for the seg(k)_time dfs
  #Impute missing time diffs with medians between two "same" time steps
  dx <- as.matrix(apply(df[,c(-1,-2)],2,as.numeric)) - 
    as.matrix(apply(df[,c(-1,-ncol(df))],2,as.numeric))
  med <- colMedians(dx, na.rm = T)
  idx_na <- is.na(dx)
  dx[idx_na] <- med[col(dx)][idx_na]
  dx <- as.data.frame(cbind(df$`subject ID`,as.data.frame(dx)))
  colnames(dx)[1] <- "CLIENT_SAMPLE_ID"
  return(dx)
}

#For merging HMDB id to prize df and thens aving as tab delimited files
prize_to_csv <- function(x, metab_info,shared_name,exclude_na=T){
  Clusters <- unique(x$Cluster)
  x <- merge(x,metab_info,all.x=T,by="CHEM_ID")
  x <- cbind(x,"terminal")
  colnames(x)[ncol(x)] <- "type"
  for (i in Clusters){
    y <- x[x$Cluster == i,c("HMDB", "prize", "type")]
    y <- cSplit(y, "HMDB", sep=",", type.convert=FALSE)
    y <- cbind(y[,3],y[,1:2])
    colnames(y)[1] <- "name"
    if(exclude_na){
      y <- y[as.logical(!is.na(y[,1])),]
    }
    write.table(y, file = paste(shared_name,"cl",formatC(i,width=2,format="d",flag="0"),'.tsv',
                                sep=""),
                sep = "\t", col.names = T, row.names = F)  
  }
}

#Load existing workspace; assuming all commands were executed previously, none of the steps hereafter are necessary
path <- "/Volumes/GoogleDrive-115111199924997198421/My Drive/Objectives/Code/"
setwd(path)
load(file=paste('ma_metabolites_eda_ts_across v0.7r.RData',sep=''))

#Load main data
file <- 'Preprocessed_UntargetedMetabolomics-ReNormalizedData-Spring2109272021.xlsx'
dat <- read_excel(file, sheet = 'Log-Transformed and Merged Data')
time_file <- 'Phlebotomy-Spring2108172021.xlsx'

#Load the timestamps
seg2_time <- read_excel(time_file, sheet = 'Segment 2')
seg4_time <- read_excel(time_file, sheet = 'Segment 4')
seg5_time <- read_excel(time_file, sheet = 'Segment 5')
seg2_time <- seg2_time[,c(1,2,4,11,14,16)]
seg4_time <- seg4_time[,c(1,2,4,11,14,16)]
seg5_time <- seg5_time[,c(1,2,4,11,14,16)]

seg2_time[seg2_time == 'NA'] <- NA
seg4_time[seg4_time == 'NA'] <- NA
seg5_time[seg5_time == 'NA'] <- NA

#seg2_time <- seg2_time[complete.cases(seg2_time),]
#seg4_time <- seg4_time[complete.cases(seg4_time),]
#seg5_time <- seg5_time[complete.cases(seg5_time),]

#Add the dates to the times so that when it's time for diff between segments, can get diff in days
seg2_time[,-c(1,2)] <- sweep(apply(as.matrix(seg2_time[,-c(1,2)]),2,as.numeric), 1, 
                             as.numeric(seg2_time$`Date of Segment 2 Completion`), "+")
seg4_time[,-c(1,2)] <- sweep(apply(as.matrix(seg4_time[,-c(1,2)]),2,as.numeric), 1, 
                             as.numeric(seg4_time$`Date of Segment 4 Completion`), "+")
seg5_time[,-c(1,2)] <- sweep(apply(as.matrix(seg5_time[,-c(1,2)]),2,as.numeric), 1, 
                             as.numeric(seg5_time$`Date of Segment 5 Completion`), "+")

#Made adjustment for a subject with a erroneous seg5 completion date
#his or her seg5 blood draw 1 time is before seg4 blood draw 4
#to adjust for this, push the seg5 completion time by one day
#this doesn't seem like a far-fetched adjustment, as GEma037 and GEma045's
#seg 4 and 5 completion days are only a day apart as well.
seg4_time[seg4_time$`subject ID` == "GEma023",]
seg5_time[seg5_time$`subject ID` == "GEma023",]
seg5_time[seg5_time$`subject ID` == "GEma023",]$`Date of Segment 5 Completion` <- 
  as.character(as.numeric(seg5_time[seg5_time$`subject ID` == "GEma023",]$`Date of Segment 5 Completion`) + 1)

#Put the blood draw (k)'s together across segments
draw1_time <- join_all(list(seg2_time[,c(1,2,3)], seg4_time[,c(1,2,3)], seg5_time[,c(1,2,3)]), by = 'subject ID')
draw1_time <- draw1_time[complete.cases(draw1_time[,c(2,4,6)]),-c(2,4,6)]
draw2_time <- join_all(list(seg2_time[,c(1,2,4)], seg4_time[,c(1,2,4)], seg5_time[,c(1,2,4)]), by = 'subject ID')
draw2_time <- draw2_time[complete.cases(draw2_time[,c(2,4,6)]),-c(2,4,6)]
draw3_time <- join_all(list(seg2_time[,c(1,2,5)], seg4_time[,c(1,2,5)], seg5_time[,c(1,2,5)]), by = 'subject ID')
draw3_time <- draw3_time[complete.cases(draw3_time[,c(2,4,6)]),-c(2,4,6)]
draw4_time <- join_all(list(seg2_time[,c(1,2,6)], seg4_time[,c(1,2,6)], seg5_time[,c(1,2,6)]), by = 'subject ID')
draw4_time <- draw4_time[complete.cases(draw4_time[,c(2,4,6)]),-c(2,4,6)]

#Calculate dx first (impute missing time diffs w/ medians)
draw1_dx <- calc_dx_across(draw1_time)
draw2_dx <- calc_dx_across(draw2_time)
draw3_dx <- calc_dx_across(draw3_time)
draw4_dx <- calc_dx_across(draw4_time)

#Sort by subject ID, blood draw, then segments
dat <- dat[order(dat$CLIENT_SAMPLE_ID),]
dat <- dat[order(dat$CUSTOM_ATTRIBUTE_2),]
dat <- dat[order(dat$CUSTOM_ATTRIBUTE_1),]

#mean(draw1[draw1$CUSTOM_ATTRIBUTE_1=="Segment 2" & draw1$CUSTOM_ATTRIBUTE_2 =="Blood Draw 1",]$'482')
#mean(draw1[draw1$CUSTOM_ATTRIBUTE_1=="Segment 4" & draw1$CUSTOM_ATTRIBUTE_2 =="Blood Draw 1",]$'482')
#mean(draw1[draw1$CUSTOM_ATTRIBUTE_1=="Segment 5" & draw1$CUSTOM_ATTRIBUTE_2 =="Blood Draw 1",]$'482')
#Z-score scale by metabolite type
lim_idx <- grep(F, sapply(dat[1,], is.numeric))[2] - 1
dat[,2:lim_idx] <- apply(dat[,2:lim_idx],2,scale)

#Replace div by zero errors w/ 0
dat[,2:lim_idx][is.nan(dat[,2:lim_idx])] <- 0
num_metabo <- lim_idx-1

#Cut up dat by blood draws
draws <- unique(dat[,'CUSTOM_ATTRIBUTE_2'])
segments <- unique(dat[,'CUSTOM_ATTRIBUTE_1'])
draw1 <- dat[dat[,'CUSTOM_ATTRIBUTE_2'] == draws[1,]$CUSTOM_ATTRIBUTE_2,]
draw2 <- dat[dat[,'CUSTOM_ATTRIBUTE_2'] == draws[2,]$CUSTOM_ATTRIBUTE_2,]
draw3 <- dat[dat[,'CUSTOM_ATTRIBUTE_2'] == draws[3,]$CUSTOM_ATTRIBUTE_2,]
draw4 <- dat[dat[,'CUSTOM_ATTRIBUTE_2'] == draws[4,]$CUSTOM_ATTRIBUTE_2,]
#rm('dat')

draw1 <- longitudify(draw1)
draw2 <- longitudify(draw2)
draw3 <- longitudify(draw3)
draw4 <- longitudify(draw4)

lim_idx1 <- grep(F, sapply(draw1[1,], is.numeric))[2] - 1
lim_idx2 <- grep(F, sapply(draw2[1,], is.numeric))[2] - 1
lim_idx3 <- grep(F, sapply(draw3[1,], is.numeric))[2] - 1
lim_idx4 <- grep(F, sapply(draw4[1,], is.numeric))[2] - 1

#Critical to have only as many subjects as there are in the seg(k) frames
draw1_dx <- merge(draw1[draw1$CUSTOM_ATTRIBUTE_1 == 'Segment 2','CLIENT_SAMPLE_ID'], 
                 draw1_dx, by = 'CLIENT_SAMPLE_ID', all.x = T)
draw2_dx <- merge(draw2[draw2$CUSTOM_ATTRIBUTE_1 == 'Segment 2','CLIENT_SAMPLE_ID'], 
                  draw2_dx, by = 'CLIENT_SAMPLE_ID', all.x = T)
draw3_dx <- merge(draw3[draw3$CUSTOM_ATTRIBUTE_1 == 'Segment 2','CLIENT_SAMPLE_ID'], 
                  draw3_dx, by = 'CLIENT_SAMPLE_ID', all.x = T)
draw4_dx <- merge(draw4[draw4$CUSTOM_ATTRIBUTE_1 == 'Segment 2','CLIENT_SAMPLE_ID'], 
                  draw4_dx, by = 'CLIENT_SAMPLE_ID', all.x = T)

#Sort data in the order of subject ids.
draw1_dx <- draw1_dx[order(draw1_dx$CLIENT_SAMPLE_ID),]
draw2_dx <- draw2_dx[order(draw2_dx$CLIENT_SAMPLE_ID),]
draw3_dx <- draw3_dx[order(draw3_dx$CLIENT_SAMPLE_ID),]
draw4_dx <- draw4_dx[order(draw4_dx$CLIENT_SAMPLE_ID),]

excl <- "GEma097" #This one doesn't have any timestamps in seg 5 and therefore should be eliminated

draw1 <- draw1[draw1$CLIENT_SAMPLE_ID!= excl, ]
draw2 <- draw2[draw2$CLIENT_SAMPLE_ID!= excl, ]
draw3 <- draw3[draw3$CLIENT_SAMPLE_ID!= excl, ]
draw4 <- draw4[draw4$CLIENT_SAMPLE_ID!= excl, ]

draw1_dx <- draw1_dx[draw1_dx$CLIENT_SAMPLE_ID!= excl, ]
draw2_dx <- draw2_dx[draw2_dx$CLIENT_SAMPLE_ID!= excl, ]
draw3_dx <- draw3_dx[draw3_dx$CLIENT_SAMPLE_ID!= excl, ]
draw4_dx <- draw4_dx[draw4_dx$CLIENT_SAMPLE_ID!= excl, ]

#For assigning column names later.
headers1 <- draw1$CUSTOM_ATTRIBUTE_1
headers2 <- draw2$CUSTOM_ATTRIBUTE_1
headers3 <- draw3$CUSTOM_ATTRIBUTE_1
headers4 <- draw4$CUSTOM_ATTRIBUTE_1

#Take the log changes (not log2fc) only and transpose them.
draw1 <- t(draw1[,2:lim_idx1])
draw2 <- t(draw2[,2:lim_idx2])
draw3 <- t(draw3[,2:lim_idx3])
draw4 <- t(draw4[,2:lim_idx4])

#Make the segment numbers the column headers.
colnames(draw1) <- headers1
colnames(draw2) <- headers2
colnames(draw3) <- headers3
colnames(draw4) <- headers4

#Want to stack all the instances of time (segment) == 2 on top of each other and so on.
#This way, we get a dataframe of analytes vs. time steps.
where4_1 <- grep("Segment 4", ignore.case = T, colnames(draw1))[1]
where5_1 <- grep("Segment 5", ignore.case = T, colnames(draw1))[1]

where4_2 <- grep("Segment 4", ignore.case = T, colnames(draw2))[1]
where5_2 <- grep("Segment 5", ignore.case = T, colnames(draw2))[1]

where4_3 <- grep("Segment 4", ignore.case = T, colnames(draw3))[1]
where5_3 <- grep("Segment 5", ignore.case = T, colnames(draw3))[1]

where4_4 <- grep("Segment 4", ignore.case = T, colnames(draw4))[1]
where5_4 <- grep("Segment 5", ignore.case = T, colnames(draw4))[1]

seg2_1 <- melt(draw1[,1:(where4_1-1)],id.vars=1)
seg4_1 <- melt(draw1[,where4_1:(where5_1-1)],id.vars=1)
seg5_1 <- melt(draw1[,where5_1:ncol(draw1)],id.vars=1)

seg2_2 <- melt(draw2[,1:(where4_2-1)],id.vars=1)
seg4_2 <- melt(draw2[,where4_2:(where5_2-1)],id.vars=1)
seg5_2 <- melt(draw2[,where5_2:ncol(draw2)],id.vars=1)

seg2_3 <- melt(draw3[,1:(where4_3-1)],id.vars=1)
seg4_3 <- melt(draw3[,where4_3:(where5_3-1)],id.vars=1)
seg5_3 <- melt(draw3[,where5_3:ncol(draw3)],id.vars=1)

seg2_4 <- melt(draw4[,1:(where4_4-1)],id.vars=1)
seg4_4 <- melt(draw4[,where4_4:(where5_4-1)],id.vars=1)
seg5_4 <- melt(draw4[,where5_4:ncol(draw4)],id.vars=1)

draw1 <- cbind(seg2_1[,3],seg4_1[,3],seg5_1[,3])
draw2 <- cbind(seg2_2[,3],seg4_2[,3],seg5_2[,3])
draw3 <- cbind(seg2_3[,3],seg4_3[,3],seg5_3[,3])
draw4 <- cbind(seg2_4[,3],seg4_4[,3],seg5_4[,3])

#Match the time stamps for each subject id x metabolite combination
draw1_dx <- draw1_dx[rep(seq_len(nrow(draw1_dx)), num_metabo),]
draw2_dx <- draw2_dx[rep(seq_len(nrow(draw2_dx)), num_metabo),]
draw3_dx <- draw3_dx[rep(seq_len(nrow(draw3_dx)), num_metabo),]
draw4_dx <- draw4_dx[rep(seq_len(nrow(draw4_dx)), num_metabo),]

#Sanity check
nrow(draw1_dx) == nrow(draw1)
nrow(draw2_dx) == nrow(draw2)
nrow(draw3_dx) == nrow(draw3)
nrow(draw4_dx) == nrow(draw4)

#The row names must be made unique because of the peculiar integrity constraints for R dataframes.
rownames(draw1) <- make.names(seg2_1[,1],unique=T)
rownames(draw2) <- make.names(seg2_2[,1],unique=T)
rownames(draw3) <- make.names(seg2_3[,1],unique=T)
rownames(draw4) <- make.names(seg2_4[,1],unique=T)

seg_labels <- c("Segment 2", "Segment 4", "Segment 5")
colnames(draw1) <- seg_labels
colnames(draw2) <- seg_labels
colnames(draw3) <- seg_labels

#Stack everything together for a unified Clustering analysis
universe <- rbind(draw1,draw2,draw3,draw4)
rownames(universe) <- make.names(sapply(rownames(universe),take_substr),unique=T)
colnames(draw1_dx)[2:3] <- c('Draw Segment 4', 'Draw Segment 5')
colnames(draw2_dx)[2:3] <- c('Draw Segment 4', 'Draw Segment 5')
colnames(draw3_dx)[2:3] <- c('Draw Segment 4', 'Draw Segment 5')
colnames(draw4_dx)[2:3] <- c('Draw Segment 4', 'Draw Segment 5')
universe_dx <- rbind(draw1_dx,draw2_dx,draw3_dx,draw4_dx)
temp_names <- universe_dx$CLIENT_SAMPLE_ID
universe_dx <- as.matrix(universe_dx[,-1])
rownames(universe_dx) <- temp_names

temp1_names <- draw1_dx$CLIENT_SAMPLE_ID
temp2_names <- draw2_dx$CLIENT_SAMPLE_ID
draw1_dx <- as.matrix(draw1_dx[,-1])
draw2_dx <- as.matrix(draw2_dx[,-1])
rownames(draw1_dx) <- temp1_names
rownames(draw2_dx) <- temp2_names

#The log scales are for plotting later on
log_draw1_dx <- log(draw1_dx * 24 * 60) #log on a minute scale
log_draw2_dx <- log(draw2_dx * 24 * 60)
log_universe_dx <- log(universe_dx * 24 * 60)

#don't want negative time diffs
sum(log_draw1_dx < 0)
sum(log_draw2_dx < 0)
sum(log_universe_dx < 0)

#Calculate the slopes
draw1_slopes <- (draw1[,-1] - draw1[,-ncol(draw1)])/draw1_dx
draw2_slopes <- (draw2[,-1] - draw2[,-ncol(draw2)])/as.matrix(draw2_dx)
universe_slopes <- (universe[,-1] - universe[,-ncol(universe)])/as.matrix(universe_dx)
log_draw1_slopes <- (draw1[,-1] - draw1[,-ncol(draw1)])/log_draw1_dx
log_draw2_slopes <- (draw2[,-1] - draw2[,-ncol(draw2)])/log_draw2_dx
log_universe_slopes <- (universe[,-1] - universe[,-ncol(draw2)])/log_universe_dx

#Quick Fix: want to look at average slope per molecule rather than individuals to deal with noisiness of the measurements
log_draw1_slopes <- cbind(log_draw1_slopes,as.numeric(sapply(rownames(log_draw1_slopes),take_substr)))
colnames(log_draw1_slopes)[3] <- "CHEM_ID" 

log_draw2_slopes <- cbind(log_draw2_slopes,as.numeric(sapply(rownames(log_draw2_slopes),take_substr)))
colnames(log_draw2_slopes)[3] <- "CHEM_ID" 

log_draw1_slopes <- aggregate(log_draw1_slopes[,1:2], list(log_draw1_slopes[,3]), mean)
log_draw2_slopes <- aggregate(log_draw2_slopes[,1:2], list(log_draw2_slopes[,3]), mean)

rownames(log_draw1_slopes) <- log_draw1_slopes[,1]
rownames(log_draw2_slopes) <- log_draw2_slopes[,1]

log_draw1_slopes <- as.matrix(log_draw1_slopes[,2:3])
log_draw2_slopes <- as.matrix(log_draw2_slopes[,2:3])

#Create timestamps
universe_time <- cbind(0,universe_dx)
colnames(universe_time) <- seg_labels
universe_time[,3] <- universe_time[,3] +  universe_time[,2]

draw1_time <- cbind(0,draw1_dx)
colnames(draw1_time) <- seg_labels
draw1_time[,3] <- draw1_time[,3] +  draw1_time[,2]

draw2_time <- cbind(0,draw2_dx)
colnames(draw2_time) <- seg_labels
draw2_time[,3] <- draw2_time[,3] +  draw2_time[,2]

log_universe_time <- cbind(0,log_universe_dx)
colnames(log_universe_time) <- seg_labels
log_universe_time[,3] <- log_universe_time[,3] + log_universe_time[,2]

log_draw1_time <- cbind(0,log_draw1_dx)
colnames(log_draw1_time) <- seg_labels
log_draw1_time[,3] <- log_draw1_time[,3] + log_draw1_time[,2]

log_draw2_time <- cbind(0,log_draw2_dx)
colnames(log_draw2_time) <- seg_labels
log_draw2_time[,3] <- log_draw2_time[,3] + log_draw2_time[,2]

#For reproducibility
set.seed(23)

#Create expression sets for Clustering
draw1_slopes_exp <- ExpressionSet(draw1_slopes)
draw2_slopes_exp <- ExpressionSet(draw2_slopes)
uni_slopes_exp <- ExpressionSet(universe_slopes)

log_draw1_slopes_exp <- ExpressionSet(log_draw1_slopes)
log_draw2_slopes_exp <- ExpressionSet(log_draw2_slopes)
log_uni_slopes_exp <- ExpressionSet(log_universe_slopes)

#Below is for dmin calculation and the associated elbow pick. - should be done before Clustering.
#More Info: The minimum centroid distance can be used as Cluster validity index. 
#For an optimal Cluster number, we may see a 'drop' of minimum centroid distance when plotted versus 
#a range of Cluster number and a slower decrease of the minimum centroid distance for higher 
#Cluster number. More information and some examples can be found in the study of 
#Schwaemmle and Jensen (2010). However, it should be used with care, as the determination remains 
#difficult especially for short time series and overlapping Clusters.
#https://doi.org/10.1093/bioinformatics/btq534
dmin_slopes_draw1 <- Dmin(eset = draw1_slopes_exp, m = 1.2, crange = 3:20, repeats = 3, visu=T) #13
dmin_slopes_draw2 <- Dmin(eset = draw2_slopes_exp, m = 1.2, crange = 3:20, repeats = 3, visu=T) #4

dmin_slopes_log_draw1 <- Dmin(eset = log_draw1_slopes_exp, m = 1.2, crange = 3:25, repeats = 3, visu=T)
dmin_slopes_log_draw2 <- Dmin(eset = log_draw2_slopes_exp, m = 1.2, crange = 3:25, repeats = 3, visu=T)

#To plot the dmin series in case visu = F.
plot(dmin_slopes_draw1, x=3:20, main = "Minimum Distance between Cluster Centers vs. Number of Clusters", ylab = "Min. Euclidean Centroid Distance",
     xlab = "Number of Clusters") #13
plot(dmin_slopes_draw2[-1], x=4:20, main = "Minimum Distance between Cluster Centers vs. Number of Clusters", ylab = "Min. Euclidean Centroid Distance",
     xlab = "Number of Clusters") #15

plot(dmin_slopes_log_draw1, x=3:25, main = "Minimum Distance between Cluster Centers vs. Number of Clusters", ylab = "Min. Euclidean Centroid Distance",
     xlab = "Number of Clusters") #11
plot(dmin_slopes_log_draw2, x=3:25, main = "Minimum Distance between Cluster Centers vs. Number of Clusters", ylab = "Min. Euclidean Centroid Distance",
     xlab = "Number of Clusters") #18

#Clustering based on the elbow "inflection" point
draw1_slopes_cl <- mfuzz(draw1_slopes_exp, c=13, m=1.1)
draw1_membership_cl <- apply(draw1_slopes_cl$membership,1,which.max)
draw1 <- cbind(draw1[,1:3],draw1_membership_cl)

draw2_slopes_cl <- mfuzz(draw2_slopes_exp, c=4, m=1.1)
draw2_membership_cl <- apply(draw2_slopes_cl$membership,1,which.max)
draw2 <- cbind(draw2[,1:3],draw2_membership_cl)

log_draw1_slopes_cl <- mfuzz(log_draw1_slopes_exp, c=11, m=1.2)
log_draw1_membership_cl <- apply(log_draw1_slopes_cl$membership,1,which.max)
log_draw1 <- cbind(draw1[,1:3],log_draw1_membership_cl)

prize_log_draw1 <- log_draw1_slopes_cl$membership
prize_log_draw1 <- cbind(as.numeric(rownames(prize_log_draw1)), prize_log_draw1)
prize_log_draw1 <- as.data.frame(prize_log_draw1)
colnames(prize_log_draw1)[1] <- "CHEM_ID"
prize_log_draw1 <- prize_log_draw1 %>% group_by(CHEM_ID) %>% summarize_each(mean)

log_draw2_slopes_cl <- mfuzz(log_draw2_slopes_exp, c=18, m=1.2)
log_draw2_membership_cl <- apply(log_draw2_slopes_cl$membership,1,which.max)
log_draw2 <- cbind(draw2[,1:3],log_draw2_membership_cl)

prize_log_draw2 <- log_draw2_slopes_cl$membership
prize_log_draw2 <- cbind(as.numeric(rownames(prize_log_draw2)), prize_log_draw2)
prize_log_draw2 <- as.data.frame(prize_log_draw2)
colnames(prize_log_draw2)[1] <- "CHEM_ID"
prize_log_draw2 <- prize_log_draw2 %>% group_by(CHEM_ID) %>% summarize_each(mean)

#standardize the time points across the board for an "easier to understand" set of plots
draw1_int <- draw1
draw1_med_time <- apply(draw1_time,2,median, na.rm = T)
draw1_med_dx <- apply(draw1_dx,2,median, na.rm = T)

draw1_int[,2] <- draw1_med_dx[1] * draw1_slopes[,1] + draw1_int[,1]
draw1_int[,3] <- draw1_med_dx[2] * draw1_slopes[,2] + draw1_int[,2]

draw2_int <- draw2
draw2_med_time <- apply(draw2_time,2,median, na.rm = T)
draw2_med_dx <- apply(draw2_dx,2,median, na.rm = T)

draw2_int[,2] <- draw2_med_dx[1] * draw2_slopes[,1] + draw2_int[,1]
draw2_int[,3] <- draw2_med_dx[2] * draw2_slopes[,2] + draw2_int[,2]


log_seg2_int <- aggregate(log_seg2,list(log_seg2$CHEM_ID),mean)[,-1]
rownames(log_seg2_int) <- log_seg2_int$CHEM_ID
log_seg2_int <- log_seg2_int[,-1]
log_seg2_med_time <- apply(log_seg2_time,2,median, na.rm = T)
log_seg2_med_dx <- apply(log_seg2_dx,2,median, na.rm = T)
log_seg2_int[,2] <- log_seg2_med_dx[1] * log_seg2_slopes[,1] + log_seg2_int[,1]
log_seg2_int[,3] <- log_seg2_med_dx[2] * log_seg2_slopes[,2] + log_seg2_int[,2]
log_seg2_int[,4] <- log_seg2_med_dx[3] * log_seg2_slopes[,3] + log_seg2_int[,3]

log_draw1 <- cbind(as.numeric(sapply(rownames(log_draw1),take_substr)),log_draw1)
log_draw2 <- cbind(as.numeric(sapply(rownames(log_draw2),take_substr)),log_draw2)
colnames(log_draw1)[c(1,5)] <- c("CHEM_ID", "Cluster")
colnames(log_draw2)[c(1,5)] <- c("CHEM_ID", "Cluster")

log_draw1_int <- aggregate(log_draw1,list(log_draw1[,1]),mean)[,-1]
rownames(log_draw1_int) <- log_draw1_int$CHEM_ID
log_draw1_int <- log_draw1_int[,-1]
log_draw1_med_time <- apply(log_draw1_time,2,median, na.rm = T)
log_draw1_med_dx <- apply(log_draw1_dx,2,median, na.rm = T)

log_draw1_int[,2] <- log_draw1_med_dx[1] * log_draw1_slopes[,1] + log_draw1_int[,1]
log_draw1_int[,3] <- log_draw1_med_dx[2] * log_draw1_slopes[,2] + log_draw1_int[,2]

log_draw2_int <- aggregate(log_draw2,list(log_draw2[,1]),mean)[,-1]
rownames(log_draw2_int) <- log_draw2_int$CHEM_ID
log_draw2_int <- log_draw2_int[,-1]
log_draw2_med_time <- apply(log_draw2_time,2,median, na.rm = T)
log_draw2_med_dx <- apply(log_draw2_dx,2,median, na.rm = T)

log_draw2_int[,2] <- log_draw2_med_dx[1] * log_draw2_slopes[,1] + log_draw2_int[,1]
log_draw2_int[,3] <- log_draw2_med_dx[2] * log_draw2_slopes[,2] + log_draw2_int[,2]

#Plot an example cluster (change params as needed)
matplot(y = t(log_draw1_int[
  sample(
    seq(nrow(log_draw1_int))[
      log_draw1_int[,4]==1],
    min(100000, sum(log_draw1_int[,4]==1))
  ),1:3]), 
  x = log_draw1_med_time,
  xaxt="n",
  type = "l", col = "red", lty = 1, 
  xlab='Time (min, log-scaled)', ylab='Z-Score', 
  main = "Blood Draw 1 - Cluster 1", cex.lab = 2, cex.main = 1.9)
axis(1,at=log_draw1_med_time,labels=round(log_draw2_med_time,1),las=1)

#Clean up the Cluster assignment dataframe
# log_draw1 <- data.frame(log_draw1)
# log_draw1 <- cbind(rownames(log_draw1),log_draw1)
# colnames(log_draw1) <- c("CHEM_ID", "Segment 2", "Segment 4", "Segment 5", "Cluster")
# log_draw1$CHEM_ID <- as.numeric(sapply(log_draw1$CHEM_ID,take_substr))

#Retrieve pathway data for the above dataframe
metab_info <- read_excel(paste(getwd(),'/',file, sep = ''), sheet = 'Chemical Annotation')
metab_info <- data.frame(metab_info[order(metab_info$CHEM_ID), c("CHEM_ID", "HMDB", "CHEMICAL_NAME", 
                                                                 "SUPER_PATHWAY", "SUB_PATHWAY")])
rownames(metab_info) <- metab_info$CHEM_ID

#Merge pathway info into the Clustering dataframe
log_draw1 <- merge(log_draw1, metab_info[,c("HMDB", "CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY")],
                  by.x = "CHEM_ID", by.y = "row.names", all.x = T)
log_draw2 <- merge(log_draw2, metab_info[,c("HMDB", "CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY")],
                   by.x = "CHEM_ID", by.y = "row.names", all.x = T)

#Set all na values to 0
log_draw1[is.na(log_draw1)] <- 0
log_draw2[is.na(log_draw2)] <- 0

p_vals_log_draw1 <- hypergeo_test(log_draw1,log_draw1, alpha = 1)
p_vals_log_draw2 <- hypergeo_test(log_draw2,log_draw2, alpha = 1)

p_vals_log_draw1 <- p_vals_log_draw1[-1270,]
p_vals_log_draw2 <- p_vals_log_draw2[-1270,]

#Eliminate chemicals that got assigned to multiple Clusters
table_chem_id <- table(p_vals_log_draw1[,1])
p_vals_log_draw1 <- p_vals_log_draw1[p_vals_log_draw1[,1]%in%names(table_chem_id[table_chem_id==1]),]

table_chem_id <- table(p_vals_log_draw2[,1])
p_vals_log_draw2 <- p_vals_log_draw2[p_vals_log_draw2[,1]%in%names(table_chem_id[table_chem_id==1]),]

#Get the prize values
prize_log_draw1$CHEM_ID <- as.character(prize_log_draw1$CHEM_ID)
prize_log_draw2$CHEM_ID <- as.character(prize_log_draw2$CHEM_ID)

#Get max as the final prize value
prize_log_draw1_max <- apply(prize_log_draw1[,-1],1,max)
prize_log_draw1_max <- as.data.frame(cbind(prize_log_draw1$CHEM_ID, prize_log_draw1_max))
colnames(prize_log_draw1_max) <- c("CHEM_ID", "Prize")

prize_log_draw2_max <- apply(prize_log_draw2[,-1],1,max)
prize_log_draw2_max <- as.data.frame(cbind(prize_log_draw2$CHEM_ID, prize_log_draw2_max))
colnames(prize_log_draw2_max) <- c("CHEM_ID", "Prize")

#Merge prize values
p_vals_log_draw1 <- merge(p_vals_log_draw1, prize_log_draw1_max, by = "CHEM_ID", all.x = T)
p_vals_log_draw2 <- merge(p_vals_log_draw2, prize_log_draw2_max, by = "CHEM_ID", all.x = T)

#Reorder p_val datasets
p_vals_log_draw1 <- p_vals_log_draw1[order(as.numeric(p_vals_log_draw1$Hochberg_p_value)),]
p_vals_log_draw1 <- p_vals_log_draw1[order(p_vals_log_draw1$SUB_PATHWAY),]

#Calculate prize values per member per cluster for Omics Integration
log_draw1_cl <- cbind(as.numeric(names(log_draw1_membership_cl)),
                     log_draw1_membership_cl)
colnames(log_draw1_cl) <- c("CHEM_ID", "Cluster")
log_draw1_cl <- merge(unique(log_draw1_cl),prize_log_draw1, by="CHEM_ID", all.x=T)
log_draw1_cl <- cbind(log_draw1_cl[,1:2],
                     as.matrix(log_draw1_cl[,c(-1,-2)])[
                       as.matrix(cbind(seq_along(log_draw1_cl$Cluster),
                                       log_draw1_cl$Cluster))])
colnames(log_draw1_cl)[3] <- 'prize'


log_draw2_cl <- cbind(as.numeric(names(log_draw2_membership_cl)),
                      log_draw2_membership_cl)
colnames(log_draw2_cl) <- c("CHEM_ID", "Cluster")
log_draw2_cl <- merge(unique(log_draw2_cl),prize_log_draw2, by="CHEM_ID", all.x=T)
log_draw2_cl <- cbind(log_draw2_cl[,1:2],
                      as.matrix(log_draw2_cl[,c(-1,-2)])[
                        as.matrix(cbind(seq_along(log_draw2_cl$Cluster),
                                        log_draw2_cl$Cluster))])
colnames(log_draw2_cl)[3] <- 'prize'

setwd(paste(getwd(),'/Prizes/Across/Draw1',sep=''))

prize_to_csv(log_draw1_cl,metab_info,"log_draw1_avg_metab_prize_")
setwd('..'); setwd(paste(getwd(),'/Draw2',sep=''))
prize_to_csv(log_draw2_cl,metab_info,"log_draw2_avg_metab_prize_")
setwd('..'); setwd('..'); setwd('..');
write.csv(p_vals_log_draw1,paste("p_vals_slope_log_draw1_avg.csv"))
write.csv(p_vals_log_draw2,paste("p_vals_slope_log_draw2_avg.csv"))

#Save workspace
save.image(file=paste('ma_metabolites_eda_ts_across v0.7r.RData',sep=''))