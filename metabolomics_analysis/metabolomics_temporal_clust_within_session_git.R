#This version has no imputations - throws out data to keep things logitudinal
#Data manipulation packs
library(readxl)
library(reshape2)
library(robustbase)
library(stringr)
library(boot)
library(splitstackshape)
library(purrr)

#Bioconductor packs
library(Mfuzz)
library(Biobase)
library(convert)
library(tscR)

#SignedTest packs
library(psych)
library(BSDA)
library(DescTools)

library(cluster)

#Mode function
getmode <- function(x) {
  unique_x <- unique(x)
  tabulate_x <- tabulate(match(x, unique_x))
  unique_x[tabulate_x == max(tabulate_x)]
}

#For unstacking the membership vectors by blood draws.
unstack <- function(x, num_analytes=1269){
  num_draws <- length(x)/num_analytes
  agg <- NULL
  for(i in 1:num_draws) { 
    sub <- x[(1:num_analytes)+(i-1)*num_analytes]
    agg <- cbind(agg,sub)
  }
  return(agg)
}

scale_by_draw <- function(dat, idx_lim, sep=''){
  #Function for scaling each molecule within each blood draw within each segment.
  for (i in 1:nrow(sep)){
    dat[dat[,'time']==sep[i,]$time,2:idx_lim] = 
      scale(dat[dat[,'time']==sep[i,]$time,2:idx_lim])
  }
  return(dat)
}

is.nan.data.frame <- function(x){
  #Because is.nan doesn't work with dataframes like is.na.
  do.call(cbind, lapply(x, is.nan))
}

standardize <- function(df){
  #Scale dataframes by column
  df <- apply(df,2,scale)
  return(df)
}

longitudify <- function(df){
  #In cases where the number of blood draws is fewer in the current session than that of 
  #another draw session, fill the gap with the average for each compound in the current draw.  
  draws <- unique(df$CUSTOM_ATTRIBUTE_2)
  len_draws <- length(draws)
  lim_idx <- grep(F, sapply(df[1,], is.numeric))[c(1,2)]
  lim_idx[1] <- lim_idx[1] + 1
  lim_idx[2] <- lim_idx[2] - 1
  

  #Take the intersection of subject ids present across all sessions.
  cap <- Reduce(intersect, list(df[df$CUSTOM_ATTRIBUTE_2==draws[1],]$CLIENT_SAMPLE_ID,
                         df[df$CUSTOM_ATTRIBUTE_2==draws[2],]$CLIENT_SAMPLE_ID,
                         df[df$CUSTOM_ATTRIBUTE_2==draws[3],]$CLIENT_SAMPLE_ID,
                         df[df$CUSTOM_ATTRIBUTE_2==draws[4],]$CLIENT_SAMPLE_ID))
  
  #Then filter the dataframe.
  df <- df[sapply(df$CLIENT_SAMPLE_ID, is.element, set = cap),]
  return(df[!duplicated(df[,c("CLIENT_SAMPLE_ID","CUSTOM_ATTRIBUTE_2")]),])
}

#Must make the individual molecule names no longer unique
take_substr <- function(x){
  last <- str_locate(x,'\\.')[1]
  first <- str_locate(x,'X')[1]
  if (is.na(last) && is.na(first))
    return(x)
  else if (is.na(last))
    return(substr(x,first+1,nchar(x)))
  else if (is.na(first))
    return(substr(x,1,last-1))
  else
    return(substr(x,first+1,last-1))
}

unique_mode <- function(vec){
  sum(is.na(vec)) == (length(vec) - 1)
}

hypergeo_test <- function(universe, select_m, alpha = 0.05){
  #Hypergeometric hypothesis testing
  #Right tailed alternative for checking the abundance of 
  #certain subpath members within a cluster vs. the rest of the dataset.
  select_sub <- unique(select_m$SUB_PATHWAY)
  p_vals <- matrix(NA,nrow(select_m[,c("SUB_PATHWAY", "cluster")]),8) #dataframes are incredibly slow; use matrices
  a <- 1
  N <- nrow(universe)
  for (i in 1:length(select_sub)){
    sub_patt <- unique(select_m[select_m$SUB_PATHWAY==select_sub[i],]$cluster)
    for (j in 1:length(sub_patt)){
      select_metabolytes <- 
        select_m[
          select_m[,'cluster']==sub_patt[j] & select_m[,'SUB_PATHWAY'] == select_sub[i],
        ]
      x <- sum(universe[,'cluster']==sub_patt[j] &
                 universe[,'SUB_PATHWAY'] == select_sub[i] & universe$CHEM_ID%in%select_metabolytes$CHEM_ID)
      m <- sum(universe[,"SUB_PATHWAY"]==select_sub[i])
      p <- phyper(q = x-1, #p(Q > q) when using upper tail
                  m = m, 
                  n = N - m, 
                  k = sum(universe[,'cluster']==sub_patt[j]), 
                  lower.tail = FALSE, log.p = FALSE
      )
      subgroup <- unique(select_metabolytes[,c("CHEM_ID", "HMDB", "CHEMICAL_NAME", "SUPER_PATHWAY")])
      for (r in 1:nrow(subgroup)){
        p_vals[a,] <- c(subgroup[r,]$CHEM_ID, subgroup[r,]$HMDB, subgroup[r,]$CHEMICAL_NAME, subgroup[r,]$SUPER_PATHWAY, 
                        select_sub[i], sub_patt[j], p, NA)
        a = a + 1
      }
    }
  }
  p_vals <- unique(p_vals)
  p_vals[,8] <- p.adjust(as.numeric(p_vals[,7]),method="hochberg") #p-val adjustment for multiple testing.
  p_vals <- p_vals[order(as.numeric(p_vals[,8])),]
  Significant <- as.numeric(p_vals[,8]) < alpha
  p_vals <- cbind(p_vals,Significant)
  #p_vals <- p_vals[p_vals[,9] == T,]
  colnames(p_vals)[1:8] <- c("CHEM_ID", "HMDB", "CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY", "cluster", 
                             "p_value", "Hochberg_p_value")
  return(p_vals)
}

#For calculating the BIC criterion
bic <- function(fit){
  #source: https://stats.stackexchange.com/questions/90769/using-bic-to-estimate-the-number-of-k-in-kmeans/251169#251169
  #number of clusters
  m = nrow(fit$centers)
  # size of the clusters
  n = fit$size
  #size of data set
  N = sum(fit$size)
  d = ncol(fit$centers)
  
  #compute variance for all clusters beforehand
  cl_var = (1.0 / (N - m) / d) * fit$withinerror
  BIC = - 0.5 * m * log(N) * (d+1)
  
  for (i in 1:m){
    BIC = BIC + sum(n[i] * log(n[i]+1e-12) - n[i] * log(N) -
                      ((n[i] * d) / 2) * log(2*pi*cl_var) -
                      ((n[i] - 1) * d/ 2))
  }
  return(-BIC)
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

sig_test <- function(seg, membership, alpha = 0.01, thresh = 0.05, test = SignTest){
  #Matches every cluster of each metabolite to its time series pattern.
  #Keys: 
  #/ means up regulation from one time point to the next
  #\ means down regulation
  #- means no meaningful change based on a paired two-sample t test.
  u <- unique(membership)
  out <- rep(NA,length(membership))
  
  for (i in 1:length(u)){
    symbol = ""
    for (j in 1:(ncol(seg)-1)){
      test_l <- test(seg[membership==i,j],
                     seg[membership==i,j+1],paired=T, 
                     alternative="less", conf.level = 1-alpha)
      test_g <- test(seg[membership==i,j],
                     seg[membership==i,j+1],paired=T, 
                     alternative="greater", conf.level = 1-alpha)
      mean_latter <- mean(seg[membership==i,j+1])
      mean_former <- mean(seg[membership==i,j])
      change <- abs((mean_latter - mean_former)/mean_former)
      if(test_l$p.value < alpha && change > thresh){
        symbol <- paste(symbol,"u",sep="")
      } 
      else if (test_g$p.value < alpha && change > thresh){
        symbol <- paste(symbol,"d",sep="")
      }
      else{
        symbol <- paste(symbol,"s",sep="")
      }
    }
    out[membership==i] <- symbol
  }
  names(out) <- names(membership)
  return(out)
}

calc_dx <- function(df){
#Calculate the delta between consecutive time steps for the seg(k)_time dfs
#Impute missing time diffs with medians between two "same" time steps
    dx <- as.matrix(apply(df[,4:ncol(df)],2,as.numeric)) - 
      as.matrix(apply(df[,3:(ncol(df)-1)],2,as.numeric))
    med <- colMedians(dx, na.rm = T)
    idx_na <- is.na(dx)
    dx[idx_na] <- med[col(dx)][idx_na]
    dx[is.na(df[,2]),] <- NA
    dx <- as.data.frame(cbind(df[,1]$`subject ID`,as.data.frame(dx)))
    colnames(dx)[1] <- "CLIENT_SAMPLE_ID"
    return(dx)
}

#For merging HMDB id to prize df and thens aving as tab delimited files
prize_to_csv <- function(x, metab_info,shared_name,exclude_na=T){
  clusters <- unique(x$cluster)
  x <- merge(x,metab_info,all.x=T,by="CHEM_ID")
  x <- cbind(x,"terminal")
  colnames(x)[ncol(x)] <- "type"
  for (i in clusters){
    y <- x[x$cluster == i,c("HMDB", "prize", "type")]
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
load(file=paste('ma_metabolites_eda_ts_within v0.7r.RData',sep=''))

#Load main dataf
path <- "/Volumes/GoogleDrive-115111199924997198421/My Drive/Objectives/Code/"
setwd(path)
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

seg2_time <- seg2_time[complete.cases(seg2_time),]
seg4_time <- seg4_time[complete.cases(seg4_time),]
seg5_time <- seg5_time[complete.cases(seg5_time),]

#How to convert existing serialized Excel dates into R dates
as.POSIXct((
  as.numeric(seg2_time[1,'Date of Segment 2 Completion']) +
    as.numeric(seg2_time[1,'Time of Blood Draw 1'])
  + 5/24) 
  *3600*24, 
  origin = "1899-12-30",
  tz = "EST")

#Calculate dx first (impute missing time diffs w/ medians)
seg2_dx <- calc_dx(seg2_time)
seg4_dx <- calc_dx(seg4_time)
seg5_dx <- calc_dx(seg5_time)

#Sort by subject ID, blood draw, then segments
dat <- dat[order(dat$CLIENT_SAMPLE_ID),]
dat <- dat[order(dat$CUSTOM_ATTRIBUTE_2),]
dat <- dat[order(dat$CUSTOM_ATTRIBUTE_1),]

dat[dat$CUSTOM_ATTRIBUTE_1 == "Segment 2" & 
      dat$CUSTOM_ATTRIBUTE_2 == "Blood Draw 4","CLIENT_SAMPLE_ID"]

#Z-score scale by metabolite type
lim_idx <- grep(F, sapply(dat[1,], is.numeric))[2] - 1
dat[,2:lim_idx] <- apply(dat[,2:lim_idx],2,scale)

#Replace div by zero errors w/ 0
dat[,2:lim_idx][is.nan(dat[,2:lim_idx])] <- 0
num_metabo <- lim_idx-1

#Cut up dat into segments
segments <- unique(dat[,'CUSTOM_ATTRIBUTE_1'])
seg2 <- dat[dat[,'CUSTOM_ATTRIBUTE_1'] == segments[1,]$CUSTOM_ATTRIBUTE_1,]
seg4 <- dat[dat[,'CUSTOM_ATTRIBUTE_1'] == segments[2,]$CUSTOM_ATTRIBUTE_1,]
seg5 <- dat[dat[,'CUSTOM_ATTRIBUTE_1'] == segments[3,]$CUSTOM_ATTRIBUTE_1,]

seg2 <- longitudify(seg2)
seg4 <- longitudify(seg4)
seg5 <- longitudify(seg5)

lim_idx2 <- grep(F, sapply(seg2[1,], is.numeric))[2] - 1
lim_idx4 <- grep(F, sapply(seg4[1,], is.numeric))[2] - 1
lim_idx5 <- grep(F, sapply(seg5[1,], is.numeric))[2] - 1

#Critical to have only as many subjects as there are in the seg(k) frames
seg2_dx <- merge(seg2[seg2$CUSTOM_ATTRIBUTE_2 == 'Blood Draw 1','CLIENT_SAMPLE_ID'], 
                 seg2_dx, by = 'CLIENT_SAMPLE_ID', all.x = T)
seg4_dx <- merge(seg4[seg4$CUSTOM_ATTRIBUTE_2 == 'Blood Draw 1','CLIENT_SAMPLE_ID'], 
                 seg4_dx, by = 'CLIENT_SAMPLE_ID', all.x = T)
seg5_dx <- merge(seg5[seg5$CUSTOM_ATTRIBUTE_2 == 'Blood Draw 1','CLIENT_SAMPLE_ID'], 
                 seg5_dx, by = 'CLIENT_SAMPLE_ID', all.x = T)

#Sort data in the order of subject ids.
seg2_dx <- seg2_dx[order(seg2_dx$CLIENT_SAMPLE_ID),]
seg4_dx <- seg4_dx[order(seg4_dx$CLIENT_SAMPLE_ID),]
seg5_dx <- seg5_dx[order(seg5_dx$CLIENT_SAMPLE_ID),]

#Because the blood draw numbers are in character format, need to handle separately.
headers2 <- seg2$CUSTOM_ATTRIBUTE_2
headers4 <- seg4$CUSTOM_ATTRIBUTE_2
headers5 <- seg5$CUSTOM_ATTRIBUTE_2

#Take the log changes (not change) only and transpose them.
seg2 <- t(seg2[,2:lim_idx2])
seg4 <- t(seg4[,2:lim_idx4])
seg5 <- t(seg5[,2:lim_idx5])

#Make the blood draw numbers the column headers.
colnames(seg2) <- headers2
colnames(seg4) <- headers4
colnames(seg5) <- headers5

seg2 <- seg2[,order(colnames(seg2))]
seg4 <- seg4[,order(colnames(seg4))]
seg5 <- seg5[,order(colnames(seg5))]

#Want to stack all the instances of the first time_step on top of each other and so on.
where2_2 <- grep("Blood Draw 2",colnames(seg2))[1]
where3_2 <- grep("Blood Draw 3",colnames(seg2))[1]
where4_2 <- grep("Blood Draw 4",colnames(seg2))[1]

where2_4 <- grep("Blood Draw 2",colnames(seg4))[1]
where3_4 <- grep("Blood Draw 3",colnames(seg4))[1]
where4_4 <- grep("Blood Draw 4",colnames(seg4))[1]

where2_5 <- grep("Blood Draw 2",colnames(seg5))[1]
where3_5 <- grep("Blood Draw 3",colnames(seg5))[1]
where4_5 <- grep("Blood Draw 4",colnames(seg5))[1]

draw1_2 <- melt(seg2[,1:(where2_2-1)],id.vars=1)
draw2_2 <- melt(seg2[,where2_2:(where3_2-1)],id.vars=1)
draw3_2 <- melt(seg2[,where3_2:(where4_2-1)],id.vars=1)
draw4_2 <- melt(seg2[,where4_2:ncol(seg2)],id.vars=1)

draw1_4 <- melt(seg4[,1:(where2_4-1)],id.vars=1)
draw2_4 <- melt(seg4[,where2_4:(where3_4-1)],id.vars=1)
draw3_4 <- melt(seg4[,where3_4:(where4_4-1)],id.vars=1)
draw4_4 <- melt(seg4[,where4_4:ncol(seg4)],id.vars=1)

draw1_5 <- melt(seg5[,1:(where2_5-1)],id.vars=1)
draw2_5 <- melt(seg5[,where2_5:(where3_5-1)],id.vars=1)
draw3_5 <- melt(seg5[,where3_5:(where4_5-1)],id.vars=1)
draw4_5 <- melt(seg5[,where4_5:ncol(seg5)],id.vars=1)

seg2 <- cbind(draw1_2[,3],draw2_2[,3],draw3_2[,3],draw4_2[,3])
seg4 <- cbind(draw1_4[,3],draw2_4[,3],draw3_4[,3],draw4_4[,3])
seg5 <- cbind(draw1_5[,3],draw2_5[,3],draw3_5[,3],draw4_5[,3])

#Match the time stamps for each subject id x metabolite combination
seg2_dx <- seg2_dx[rep(seq_len(nrow(seg2_dx)), num_metabo),]
seg4_dx <- seg2_dx[rep(seq_len(nrow(seg4_dx)), num_metabo),]
seg5_dx <- seg2_dx[rep(seq_len(nrow(seg5_dx)), num_metabo),]

#The row names must be made unique because of the peculiar integrity constraints in R.
rownames(seg2) <- make.names(draw1_2[,1],unique=T)
rownames(seg4) <- make.names(draw1_4[,1],unique=T)
rownames(seg5) <- make.names(draw1_5[,1],unique=T)

bd_labels <- c("Blood Draw 1", "Blood Draw2", "Blood Draw 3", "Blood Draw 4")

colnames(seg2) <- bd_labels
colnames(seg4) <- bd_labels
colnames(seg5) <- bd_labels

#Stack everything together for a unified clustering analysis
universe <- rbind(seg2,seg4,seg5)
rownames(universe) <- make.names(sapply(rownames(universe),take_substr),unique=T)
universe_dx <- rbind(seg2_dx,seg4_dx,seg5_dx)
temp_names <- universe_dx$CLIENT_SAMPLE_ID
universe_dx <- as.matrix(universe_dx[,-1])
rownames(universe_dx) <- temp_names

temp_names <- seg2_dx[,1]
seg2_dx <- as.matrix(apply(seg2_dx[,-1],2,as.numeric))
rownames(seg2_dx) <- temp_names

temp_names <- seg4_dx[,1]
seg4_dx <- as.matrix(apply(seg4_dx[,-1],2,as.numeric))
rownames(seg4_dx) <- temp_names

temp_names <- seg5_dx[,1]
seg5_dx <- as.matrix(apply(seg5_dx[,-1],2,as.numeric))
rownames(seg5_dx) <- temp_names

#The log scales are for plotting later on
log_universe_dx <- log(universe_dx * 24 * 60) #log on a minute scale
log_seg2_dx <- log(seg2_dx * 24 * 60)
log_seg4_dx <- log(seg4_dx * 24 * 60)
log_seg5_dx <- log(seg5_dx * 24 * 60)

#Calculate the slopes
universe_slopes <- (universe[,-1] - universe[,-ncol(universe)])/as.matrix(universe_dx)
log_universe_slopes <- (universe[,-1] - universe[,-ncol(universe)])/as.matrix(log_universe_dx)

seg2_slopes <- (seg2[,-1] - seg2[,-ncol(seg2)])/as.matrix(seg2_dx)
log_seg2_slopes <- (seg2[,-1] - seg2[,-ncol(seg2)])/as.matrix(log_seg2_dx)

seg4_slopes <- (seg4[,-1] - seg4[,-ncol(seg4)])/as.matrix(seg4_dx)
log_seg4_slopes <- (seg4[,-1] - seg4[,-ncol(seg4)])/as.matrix(log_seg4_dx)

seg5_slopes <- (seg5[,-1] - seg5[,-ncol(seg5)])/as.matrix(seg5_dx)
log_seg5_slopes <- (seg5[,-1] - seg5[,-ncol(seg5)])/as.matrix(log_seg5_dx)

#Quick Fix: want to look at average slope per molecule rather than individuals to deal with noisiness of the measurements
log_seg2_slopes <- cbind(log_seg2_slopes,as.numeric(sapply(rownames(log_seg2_slopes),take_substr)))
colnames(log_seg2_slopes)[4] <- "CHEM_ID" 

log_seg4_slopes <- cbind(log_seg4_slopes,as.numeric(sapply(rownames(log_seg4_slopes),take_substr)))
colnames(log_seg4_slopes)[4] <- "CHEM_ID" 

log_seg5_slopes <- cbind(log_seg5_slopes,as.numeric(sapply(rownames(log_seg5_slopes),take_substr)))
colnames(log_seg5_slopes)[4] <- "CHEM_ID" 

log_seg2_slopes <- aggregate(log_seg2_slopes[,1:3], list(log_seg2_slopes[,4]), mean)
log_seg4_slopes <- aggregate(log_seg4_slopes[,1:3], list(log_seg4_slopes[,4]), mean)
log_seg5_slopes <- aggregate(log_seg5_slopes[,1:3], list(log_seg5_slopes[,4]), mean)

rownames(log_seg2_slopes) <- log_seg2_slopes[,1]
rownames(log_seg4_slopes) <- log_seg4_slopes[,1]
rownames(log_seg5_slopes) <- log_seg5_slopes[,1]

log_seg2_slopes <- as.matrix(log_seg2_slopes[,2:4])
log_seg4_slopes <- as.matrix(log_seg4_slopes[,2:4])
log_seg5_slopes <- as.matrix(log_seg5_slopes[,2:4])

#Create timestamps
universe_time <- cbind(rep(0,nrow(universe_dx)),universe_dx) * 24 * 60
colnames(universe_time) <- bd_labels
universe_time[,3] <- universe_time[,3] +  universe_time[,2]
universe_time[,4] <- universe_time[,4] +  universe_time[,3]

log_universe_time <- cbind(rep(0,nrow(log_universe_dx)),log_universe_dx)
colnames(log_universe_time) <- bd_labels
log_universe_time[,3] <- log_universe_time[,3] +  log_universe_time[,2]
log_universe_time[,4] <- log_universe_time[,4] +  log_universe_time[,3]


seg2_time <- cbind(rep(0,nrow(seg2_dx)),seg2_dx) * 24 * 60
colnames(seg2_time) <- bd_labels
seg2_time[,3] <- seg2_time[,3] +  seg2_time[,2]
seg2_time[,4] <- seg2_time[,4] +  seg2_time[,3]

log_seg2_time <- cbind(rep(0,nrow(log_seg2_dx)),log_seg2_dx)
colnames(log_seg2_time) <- bd_labels
log_seg2_time[,3] <- log_seg2_time[,3] +  log_seg2_time[,2]
log_seg2_time[,4] <- log_seg2_time[,4] +  log_seg2_time[,3]


seg4_time <- cbind(rep(0,nrow(seg4_dx)),seg4_dx) * 24 * 60
colnames(seg4_time) <- bd_labels
seg4_time[,3] <- seg4_time[,3] +  seg4_time[,2]
seg4_time[,4] <- seg4_time[,4] +  seg4_time[,3]

log_seg4_time <- cbind(rep(0,nrow(log_seg4_dx)),log_seg4_dx)
colnames(log_seg4_time) <- bd_labels
log_seg4_time[,3] <- log_seg4_time[,3] +  log_seg4_time[,2]
log_seg4_time[,4] <- log_seg4_time[,4] +  log_seg4_time[,3]


seg5_time <- cbind(rep(0,nrow(seg5_dx)),seg5_dx) * 24 * 60
colnames(seg5_time) <- bd_labels
seg5_time[,3] <- seg5_time[,3] +  seg5_time[,2]
seg5_time[,4] <- seg5_time[,4] +  seg5_time[,3]

log_seg5_time <- cbind(rep(0,nrow(log_seg5_dx)),log_seg5_dx)
colnames(log_seg5_time) <- bd_labels
log_seg5_time[,3] <- log_seg5_time[,3] +  log_seg5_time[,2]
log_seg5_time[,4] <- log_seg5_time[,4] +  log_seg5_time[,3]

#For reproducibility
set.seed(23)

#######Running the same on the "universe"#########
uni_slopes_exp <- ExpressionSet(universe_slopes)
log_uni_slopes_exp <- ExpressionSet(log_universe_slopes)

seg2_slopes_exp <- ExpressionSet(seg2_slopes)
log_seg2_slopes_exp <- ExpressionSet(log_seg2_slopes)

seg4_slopes_exp <- ExpressionSet(seg4_slopes)
log_seg4_slopes_exp <- ExpressionSet(log_seg4_slopes)

seg5_slopes_exp <- ExpressionSet(seg5_slopes)
log_seg5_slopes_exp <- ExpressionSet(log_seg5_slopes)

#Below is for dmin calculation and the associated elbow pick. - should be done before clustering.
#More Info: The minimum centroid distance can be used as cluster validity index. 
#For an optimal cluster number, we may see a 'drop' of minimum centroid distance when plotted versus 
#a range of cluster number and a slower decrease of the minimum centroid distance for higher 
#cluster number. More information and some examples can be found in the study of 
#Schwaemmle and Jensen (2010). However, it should be used with care, as the determination remains 
#difficult especially for short time series and overlapping clusters.
#https://doi.org/10.1093/bioinformatics/btq534
dmin_slopes_uni <- Dmin(eset = uni_slopes_exp, m = 1.2, crange = 3:20, repeats = 3, visu=T)
#dmin_slopes_uni <- as.numeric(read.csv('dmin_uni.csv',header = T)[,-1])
dmin_slopes_log_uni <- Dmin(eset = log_uni_slopes_exp, m = 1.2, crange = 3:25, repeats = 3, visu = T)

dmin_slopes_log_seg2 <- Dmin(eset = log_seg2_slopes_exp, m = 1.2, crange = 3:25, repeats = 3, visu = T)
dmin_slopes_log_seg4 <- Dmin(eset = log_seg4_slopes_exp, m = 1.2, crange = 3:25, repeats = 3, visu = T)
dmin_slopes_log_seg5 <- Dmin(eset = log_seg5_slopes_exp, m = 1.2, crange = 3:25, repeats = 3, visu = T)
#12

#To plot the dmin series in case visu = F.
plot(dmin_slopes_uni, x=3:20, main = "Minimum Distance between Cluster Centers vs. Number of Clusters", 
     ylab = "Min. Euclidean Pairwise Centroid Distance",
     xlab = "Number of Clusters")

plot(dmin_slopes_log_uni, x=3:25, main = "Minimum Pairwise Centorid Distance vs. Number of Clusters", 
     ylab = "Min. Pairwise Centroid Distance",
     xlab = "Number of Clusters", cex.lab = 2, cex.main = 1.9)

plot(dmin_slopes_log_seg2, x=3:25, main = "Minimum Pairwise Centorid Distance vs. Number of Clusters", 
     ylab = "Min. Pairwise Centroid Distance",
     xlab = "Number of Clusters", cex.lab = 2, cex.main = 1.9)
#8 for seg2 #9 for 1.2
#9 for seg4
#13 for seg5

#Clustering based on the elbow "inflection" point
uni_slopes_cl10 <- mfuzz(uni_slopes_exp, c=10, m=1.1)
membership_cl10 <- apply(uni_slopes_cl10$membership,1,which.max)
universe <- cbind(universe[,1:4],membership_cl10)

log_uni_slopes_cl12 <- mfuzz(log_uni_slopes_exp, c=12, m=1.1)
log_membership_cl12 <- apply(log_uni_slopes_cl12$membership,1,which.max)
log_universe <- cbind(universe[,1:4],log_membership_cl12)

log_seg2_slopes_cl12 <- mfuzz(log_seg2_slopes_exp, c=12, m=1.1)
log_seg2_membership_cl12 <- apply(log_seg2_slopes_cl12$membership,1,which.max)
log_seg2 <- cbind(seg2[,1:4],log_seg2_membership_cl12)

log_seg4_slopes_cl12 <- mfuzz(log_seg4_slopes_exp, c=12, m=1.1)
log_seg4_membership_cl12 <- apply(log_seg4_slopes_cl12$membership,1,which.max)
log_seg4 <- cbind(seg4[,1:4],log_seg4_membership_cl12)

log_seg5_slopes_cl12 <- mfuzz(log_seg5_slopes_exp, c=12, m=1.1)
log_seg5_membership_cl12 <- apply(log_seg5_slopes_cl12$membership,1,which.max)
log_seg5 <- cbind(seg5[,1:4],log_seg5_membership_cl12)

log_seg2_slopes_cl9 <- mfuzz(log_seg2_slopes_exp, c=9, m=1.2)
log_seg2_membership_cl9 <- apply(log_seg2_slopes_cl9$membership,1,which.max)
log_seg2_membership_cl9 <- cbind(as.numeric(names(log_seg2_membership_cl9)),log_seg2_membership_cl9)
colnames(log_seg2_membership_cl9) <- c("CHEM_ID","cluster")
seg2 <- cbind(seg2,as.numeric(sapply(rownames(seg2),take_substr)))
colnames(seg2)[5] <- "CHEM_ID"
log_seg2 <- merge(seg2,log_seg2_membership_cl9, by = "CHEM_ID")

log_seg4_slopes_cl9 <- mfuzz(log_seg4_slopes_exp, c=9, m=1.2)
log_seg4_membership_cl9 <- apply(log_seg4_slopes_cl9$membership,1,which.max)
log_seg4_membership_cl9 <- cbind(as.numeric(names(log_seg4_membership_cl9)),log_seg4_membership_cl9)
colnames(log_seg4_membership_cl9) <- c("CHEM_ID","cluster")
seg4 <- cbind(seg4,as.numeric(sapply(rownames(seg4),take_substr)))
colnames(seg4)[5] <- "CHEM_ID"
log_seg4 <- merge(seg4,log_seg4_membership_cl9, by = "CHEM_ID")

log_seg5_slopes_cl13 <- mfuzz(log_seg5_slopes_exp, c=13, m=1.2)
log_seg5_membership_cl13 <- apply(log_seg5_slopes_cl13$membership,1,which.max)
log_seg5_membership_cl13 <- cbind(as.numeric(names(log_seg5_membership_cl13)),log_seg5_membership_cl13)
colnames(log_seg5_membership_cl13) <- c("CHEM_ID","cluster")
seg5 <- cbind(seg5,as.numeric(sapply(rownames(seg5),take_substr)))
colnames(seg5)[5] <- "CHEM_ID"
log_seg5 <- merge(seg5,log_seg5_membership_cl13, by = "CHEM_ID")

prize_log_seg2 <- log_seg2_slopes_cl9$membership
prize_log_seg2 <- cbind(as.numeric(sapply(rownames(prize_log_seg2),take_substr)), prize_log_seg2)
prize_log_seg2 <- as.data.frame(prize_log_seg2)
colnames(prize_log_seg2)[1] <- "CHEM_ID"
prize_log_seg2 <- prize_log_seg2 %>% group_by(CHEM_ID) %>% summarize_each(mean)

prize_log_seg4 <- log_seg4_slopes_cl9$membership
prize_log_seg4 <- cbind(as.numeric(sapply(rownames(prize_log_seg4),take_substr)), prize_log_seg4)
prize_log_seg4 <- as.data.frame(prize_log_seg4)
colnames(prize_log_seg4)[1] <- "CHEM_ID"
prize_log_seg4 <- prize_log_seg4 %>% group_by(CHEM_ID) %>% summarize_each(mean)

prize_log_seg5 <- log_seg5_slopes_cl13$membership
prize_log_seg5 <- cbind(as.numeric(sapply(rownames(prize_log_seg5),take_substr)), prize_log_seg5)
prize_log_seg5 <- as.data.frame(prize_log_seg5)
colnames(prize_log_seg5)[1] <- "CHEM_ID"
prize_log_seg5 <- prize_log_seg5 %>% group_by(CHEM_ID) %>% summarize_each(mean)

#standardize the time points across the board for an "easier to understand" set of plots
universe_int <- universe
med_time <- apply(universe_time,2,median, na.rm = T)
med_dx <- apply(universe_dx,2,median, na.rm = T)

universe_int[,2] <- med_dx[1] * universe_slopes[,1] + universe_int[,1]
universe_int[,3] <- med_dx[2] * universe_slopes[,2] + universe_int[,2]
universe_int[,4] <- med_dx[3] * universe_slopes[,3] + universe_int[,3]


seg2_int <- seg2
seg2_med_time <- apply(seg2_time,2,median, na.rm = T)
seg2_med_dx <- apply(seg2_dx,2,median, na.rm = T)

seg2_int[,2] <- seg2_med_dx[1] * seg2_slopes[,1] + seg2_int[,1]
seg2_int[,3] <- seg2_med_dx[2] * seg2_slopes[,2] + seg2_int[,2]
seg2_int[,4] <- seg2_med_dx[3] * seg2_slopes[,3] + seg2_int[,3]


seg4_int <- seg4
seg4_med_time <- apply(seg4_time,2,median, na.rm = T)
seg4_med_dx <- apply(seg4_dx,2,median, na.rm = T)

seg4_int[,2] <- seg4_med_dx[1] * seg4_slopes[,1] + seg4_int[,1]
seg4_int[,3] <- seg4_med_dx[2] * seg4_slopes[,2] + seg4_int[,2]
seg4_int[,4] <- seg4_med_dx[3] * seg4_slopes[,3] + seg4_int[,3]


seg5_int <- seg5
seg5_med_time <- apply(seg5_time,2,median, na.rm = T)
seg5_med_dx <- apply(seg5_dx,2,median, na.rm = T)

seg5_int[,2] <- seg5_med_dx[1] * seg5_slopes[,1] + seg5_int[,1]
seg5_int[,3] <- seg5_med_dx[2] * seg5_slopes[,2] + seg5_int[,2]
seg5_int[,4] <- seg5_med_dx[3] * seg5_slopes[,3] + seg5_int[,3]

#standardize the time points across the board for an "easier to understand" set of plots
log_universe_int <- log_universe
log_med_time <- apply(log_universe_time,2,median, na.rm = T)
log_med_dx <- apply(log_universe_dx,2,median, na.rm = T)

log_universe_int[,2] <- log_med_dx[1] * log_universe_slopes[,1] + log_universe_int[,1]
log_universe_int[,3] <- log_med_dx[2] * log_universe_slopes[,2] + log_universe_int[,2]
log_universe_int[,4] <- log_med_dx[3] * log_universe_slopes[,3] + log_universe_int[,3]


log_seg2_int <- aggregate(log_seg2,list(log_seg2$CHEM_ID),mean)[,-1]
rownames(log_seg2_int) <- log_seg2_int$CHEM_ID
log_seg2_int <- log_seg2_int[,-1]
log_seg2_med_time <- apply(log_seg2_time,2,median, na.rm = T)
log_seg2_med_dx <- apply(log_seg2_dx,2,median, na.rm = T)
log_seg2_int[,2] <- log_seg2_med_dx[1] * log_seg2_slopes[,1] + log_seg2_int[,1]
log_seg2_int[,3] <- log_seg2_med_dx[2] * log_seg2_slopes[,2] + log_seg2_int[,2]
log_seg2_int[,4] <- log_seg2_med_dx[3] * log_seg2_slopes[,3] + log_seg2_int[,3]


log_seg4_int <- aggregate(log_seg4,list(log_seg4$CHEM_ID),mean)[,-1]
rownames(log_seg4_int) <- log_seg4_int$CHEM_ID
log_seg4_int <- log_seg4_int[,-1]
log_seg4_med_time <- apply(log_seg4_time,2,median, na.rm = T)
log_seg4_med_dx <- apply(log_seg4_dx,2,median, na.rm = T)
log_seg4_int[,2] <- log_seg4_med_dx[1] * log_seg4_slopes[,1] + log_seg4_int[,1]
log_seg4_int[,3] <- log_seg4_med_dx[2] * log_seg4_slopes[,2] + log_seg4_int[,2]
log_seg4_int[,4] <- log_seg4_med_dx[3] * log_seg4_slopes[,3] + log_seg4_int[,3]


log_seg5_int <- aggregate(log_seg5,list(log_seg5$CHEM_ID),mean)[,-1]
rownames(log_seg5_int) <- log_seg5_int$CHEM_ID
log_seg5_int <- log_seg5_int[,-1]
log_seg5_med_time <- apply(log_seg5_time,2,median, na.rm = T)
log_seg5_med_dx <- apply(log_seg5_dx,2,median, na.rm = T)
log_seg5_int[,2] <- log_seg5_med_dx[1] * log_seg5_slopes[,1] + log_seg5_int[,1]
log_seg5_int[,3] <- log_seg5_med_dx[2] * log_seg5_slopes[,2] + log_seg5_int[,2]
log_seg5_int[,4] <- log_seg5_med_dx[3] * log_seg5_slopes[,3] + log_seg5_int[,3]

par(mar = c(5,5,4,2) + 0.1)
#Plot clusters for uni
matplot(y = t(universe_int[
  sample(
    seq(nrow(universe_int))[
      universe_int[,5]==11],
        min(5000, sum(universe_int[,5]==11))
    ),1:4]), 
        x = med_time,
        xaxt="n",
        type = "l", col = "red", lty = 1, 
        xlab='Time (min)',ylab='Z-Score of Log Abundance', 
        main = "Cluster 1")
axis(1,at=med_time,labels=med_time,las=1)


#Ad-hoc plot: want Seg2 cholesterol levels stratified by VO2Max
#First find the subjects that belong to segment 2 across all BDs
seg2subjects <- rownames(log_seg2_dx[seg2$CHEM_ID==266,1:3])

#Then isolate the part of Seg2 that corresponds to cholesterol
seg2chol <- seg2[seg2$CHEM_ID==266,1:4]
row.names(seg2chol) <- seg2subjects
log_seg2chol_dx <- log_seg2_dx[seg2$CHEM_ID==266,1:3]
log_seg2chol_slopes <- (seg2chol[,-1]-seg2chol[,-ncol(seg2chol)])/
  as.matrix(log_seg2chol_dx)

#Find the median time differences
log_seg2chol_med_dx <- apply(log_seg2chol_dx,2, median, na.rm = T)

#Create median timestamps
log_seg2chol_med_time <- c(0, log_seg2chol_med_dx)
names(log_seg2chol_med_time) <- bd_labels
log_seg2chol_med_time[3] <- log_seg2chol_med_time[3] +  log_seg2chol_med_time[2]
log_seg2chol_med_time[4] <- log_seg2chol_med_time[4] +  log_seg2chol_med_time[3]

#Interpolate the y values for the standardized times
log_seg2chol_int <- seg2chol
log_seg2chol_int[,2] <- log_seg2chol_med_dx[1] * log_seg2chol_slopes[,1] + log_seg2chol_int[,1]
log_seg2chol_int[,3] <- log_seg2chol_med_dx[2] * log_seg2chol_slopes[,2] + log_seg2chol_int[,2]
log_seg2chol_int[,4] <- log_seg2chol_med_dx[3] * log_seg2chol_slopes[,3] + log_seg2chol_int[,3]

rownames(log_seg2chol_int) <- seg2subjects

#Merge in the VO2 data while quatifying the ordinal bin data
log_seg2chol_int <- merge(log_seg2chol_int,
                          vo2[,c("subject ID", "Gender", "VO2 Rank (Bin)")],
                          by.x = 0, by.y = "subject ID")
row.names(log_seg2chol_int) <- log_seg2chol_int$Row.names
log_seg2chol_int <- log_seg2chol_int[,-1]
log_seg2chol_int[,6][log_seg2chol_int[,6]=="Very Poor"] <- 1
log_seg2chol_int[,6][log_seg2chol_int[,6]=="Poor"] <- 2
log_seg2chol_int[,6][log_seg2chol_int[,6]=="Fair"] <- 3
log_seg2chol_int[,6][log_seg2chol_int[,6]=="Good"] <- 4
log_seg2chol_int[,6][log_seg2chol_int[,6]=="Excellent"] <- 5
log_seg2chol_int[,6][log_seg2chol_int[,6]=="Superior"] <- 6
log_seg2chol_int[,6] <- as.numeric(log_seg2chol_int[,6])

#Plot
matplot(y=t(log_seg2chol_int[log_seg2chol_int[,4]>-3 &
                               log_seg2chol_int[,6]%in%c(4),1:4]),
        x = log_seg2chol_med_time,
        xaxt="n",
        type = "l", 
        col = "red", 
        lty = 1, 
        xlab='Time (min, standardized and log-scaled)',
        ylab='Z-Score of Log Abundance', 
        main = "Segment 2 Cholesterol (HMDB0000067) 
        Trends For \"Good\"",
        cex.lab = 1.5,
        cex.main = 1.5)
axis(1,at=log_seg2chol_med_time,labels=round(log_seg2chol_med_time,2),las=1)

#Create an additional dataframe for stratified wilcoxon tests
seg2chol <- merge(seg2chol, vo2[, c("subject ID", "Gender", "VO2 Rank (Bin)")], 
                  by.x = 0, by.y = "subject ID")
rownames(seg2chol) <- seg2chol$Row.names
seg2chol <- seg2chol[,-1]
seg2chol$`VO2 Rank (Bin)` <- log_seg2chol_int[,6]
seg2chol <- seg2chol[order(row.names(seg2chol)),]

wilcox.test(seg2chol[seg2chol$`VO2 Rank (Bin)`%in%c(1,2,3),3],
              seg2chol[seg2chol$`VO2 Rank (Bin)`%in%c(1,2,3),1],
              paired=T)

seg2chol_time <- seg2_time[seg2$CHEM_ID==266,]
seg2chol_time <- seg2chol_time[order(row.names(seg2chol_time)),]

mean(seg2chol_time[seg2chol$`VO2 Rank (Bin)`%in%c(5,6),3]
     -seg2chol_time[seg2chol$`VO2 Rank (Bin)`%in%c(5,6),1])

#Plot clusters for log
par(mfrow=c(1,1))
par(mar = c(5,5,4,2) + 0.1)

matplot(y = t(log_seg5_int[
  sample(
    seq(nrow(log_seg5_int))[
      log_seg5_int[,5]==13],
    min(10000, sum(log_seg5_int[,5]==13))
  ),1:4]), 
  x = log_seg5_med_time,
  xaxt="n",
  type = "l", col = "red", lty = 1, 
  xlab='Time (min, standardized and log-scaled)', ylab='Z-Score', 
  main = "Segment 5 - Cluster 13",
  cex.lab = 1.8,
  cex.main = 1.9)
axis(1,at=log_seg5_med_time,labels=round(log_med_time,1),las=1)

matplot(y = t(log_seg2_int[log_seg2_int[,5]==9,1:4]), 
  x = log_seg2_med_time,
  xaxt="n",
  type = "l", col = "red", lty = 1, 
  xlab='Time (min, standardized and log-scaled)', ylab='Z-Score of Log Abundance', 
  main = "Segment 2 - Cluster 9",
  cex.lab = 1.8,
  cex.main = 1.9)
axis(1,at=log_seg2_med_time,labels=round(log_med_time,1),las=1)

wilcox.test(x=log_seg2_int[log_seg2_int[,5]==9,2], 
            log_seg2_int[log_seg2_int[,5]==9,4],
            paired = TRUE, alternative = "two.sided")

top_decile <- log_seg2_int[log_seg2_int$log_seg2_membership_cl12==5 & 
                             log_seg2_int$`5`>=
                             quantile(log_seg2_int$`5`[log_seg2_int$log_seg2_membership_cl12==5], 
                                      0.99), 
                           2:5]

bottom_decile <- log_seg2_int[log_seg2_int$log_seg2_membership_cl12==5 & 
                                log_seg2_int$`5`<
                                quantile(log_seg2_int$`5`[log_seg2_int$log_seg2_membership_cl12==5],
                                         0.1),
                              2:5]

top_decile <- top_decile + rnorm(nrow(top_decile),0,5)
bottom_decile <- bottom_decile + rnorm(nrow(top_decile),0,5)

matplot(y = t(top_decile), 
        x = log_seg2_med_time,
        xaxt="n",
        type = "l", col = "red", lty = 1, 
        xlab='Time (min, standardized and log-scaled)', ylab=expression(paste('Z-Score of Log Abundance + ',epsilon,sep =" ")), 
        main = "Cluster 5 (Top 1% in Prize Value)",
        cex.lab = 1.8,
        cex.main = 1.9, ylim = c(-15,15))
axis(1,at=log_seg2_med_time,labels=round(log_med_time,1),las=1)

matplot(y = t(bottom_decile), 
  x = log_seg2_med_time,
  xaxt="n",
  type = "l", col = "dodgerblue", lty = 1, 
  xlab='Time (min, standardized and log-scaled)', ylab=expression(paste('Z-Score of Log Abundance + ',epsilon,sep =" ")), 
  main = "Cluster 5 (Lowest 1% in Prize Value)",
  cex.lab = 1.8,
  cex.main = 1.9)
axis(1,at=log_seg2_med_time,labels=round(log_med_time,1),las=1)

#Plot a cluster for a deck visual
matplot(y = t(universe[log_universe_int[,5]==1,1:4]), 
  x = t(universe_time[log_universe_int[,5]==1,1:4]),
  type = "l", col = "red", lty = 1, 
  xlab='Time (min)', ylab='Z-Score of Log Abundance', 
  main = "Cluster 5 (Loewest 1% in Prize Value)",
  cex.lab = 2,
  cex.main = 1.9)
axis(1,at=log_seg2_med_time,labels=round(log_med_time,1),las=1)

#axis(1,at=log_med_time,labels=round(log_med_time,1),las=1)
#xaxp = c(1, 4, 3)
#Clean up the cluster assignment dataframes
universe <- data.frame(universe)
universe <- cbind(rownames(universe),universe)
colnames(universe) <- c("CHEM_ID", "Blood Draw 1", "Blood Draw 2", "Blood Draw 3", "Blood Draw 4", "cluster")
universe$CHEM_ID <- as.numeric(sapply(universe$CHEM_ID,take_substr))

log_universe <- data.frame(log_universe)
log_universe <- cbind(rownames(log_universe),log_universe)
colnames(log_universe) <- c("CHEM_ID", "Blood Draw 1", "Blood Draw 2", "Blood Draw 3", "Blood Draw 4", "cluster")
log_universe$CHEM_ID <- as.numeric(sapply(log_universe$CHEM_ID,take_substr))


log_seg2 <- data.frame(log_seg2)
log_seg2 <- cbind(rownames(log_seg2),log_seg2)
colnames(log_seg2) <- c("CHEM_ID", "Blood Draw 1", "Blood Draw 2", "Blood Draw 3", "Blood Draw 4", "cluster")
log_seg2$CHEM_ID <- as.numeric(sapply(log_seg2$CHEM_ID,take_substr))


log_seg4 <- data.frame(log_seg4)
log_seg4 <- cbind(rownames(log_seg4),log_seg4)
colnames(log_seg4) <- c("CHEM_ID", "Blood Draw 1", "Blood Draw 2", "Blood Draw 3", "Blood Draw 4", "cluster")
log_seg4$CHEM_ID <- as.numeric(sapply(log_seg4$CHEM_ID,take_substr))


log_seg5 <- data.frame(log_seg5)
log_seg5 <- cbind(rownames(log_seg5),log_seg5)
colnames(log_seg5) <- c("CHEM_ID", "Blood Draw 1", "Blood Draw 2", "Blood Draw 3", "Blood Draw 4", "cluster")
log_seg5$CHEM_ID <- as.numeric(sapply(log_seg5$CHEM_ID,take_substr))

#Retrieve pathway data for the above dataframe
file <- '/Preprocessed_UntargetedMetabolomics-ReNormalizedData-Spring2109272021.xlsx'
metab_info <- read_excel(paste(getwd(),file, sep = ''), sheet = 'Chemical Annotation')
metab_info <- data.frame(metab_info[order(metab_info$CHEM_ID), c("CHEM_ID", "HMDB", "CHEMICAL_NAME", 
                                                                 "SUPER_PATHWAY", "SUB_PATHWAY")])
rownames(metab_info) <- metab_info$CHEM_ID

#Merge pathway info into the clustering dataframe
universe <- merge(universe, metab_info[,c("HMDB", "CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY")],
                  by.x = "CHEM_ID", by.y = "row.names", all.x = T)

log_universe <- merge(log_universe, metab_info[,c("HMDB", "CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY")],
                  by.x = "CHEM_ID", by.y = "row.names", all.x = T)

log_seg2 <- merge(log_seg2, metab_info[,c("HMDB", "CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY")],
                      by.x = "CHEM_ID", by.y = "row.names", all.x = T)

log_seg4 <- merge(log_seg4, metab_info[,c("HMDB", "CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY")],
                      by.x = "CHEM_ID", by.y = "row.names", all.x = T)

log_seg5 <- merge(log_seg5, metab_info[,c("HMDB", "CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY")],
                      by.x = "CHEM_ID", by.y = "row.names", all.x = T)

#Set all na values to 0
universe[is.na(universe)] <- 0
log_universe[is.na(log_universe)] <- 0
log_seg2[is.na(log_seg2)] <- 0
log_seg4[is.na(log_seg4)] <- 0
log_seg5[is.na(log_seg5)] <- 0

log_seg2 <- log_seg2[!duplicated(log_seg2[,c("CHEM_ID","Membership")]),]
log_seg4 <- log_seg4[!duplicated(log_seg4[,c("CHEM_ID","Membership")]),]
log_seg5 <- log_seg5[!duplicated(log_seg5[,c("CHEM_ID","Membership")]),]

colnames(log_seg2)[6] <- "cluster"
colnames(log_seg4)[6] <- "cluster"
colnames(log_seg5)[6] <- "cluster"

p_vals <- hypergeo_test(universe,universe)
p_vals_log <- hypergeo_test(log_universe, log_universe)

p_vals_log_seg2 <- hypergeo_test(log_seg2, log_seg2)
p_vals_log_seg4 <- hypergeo_test(log_seg4, log_seg4)
p_vals_log_seg5 <- hypergeo_test(log_seg5, log_seg5)

seg2_sample <- read.csv('seg2_sample.csv',header=F)
seg4_sample <- read.csv('seg4_sample.csv',header=F)

log_seg2_sample <- log_seg2[log_seg2$HMDB%in%seg2_sample$V1,]
log_seg4_sample <- log_seg4[log_seg4$HMDB%in%seg4_sample$V1,]

log_seg2_sample[is.na(log_seg2_sample)] <- 0
log_seg4_sample[is.na(log_seg4_sample)] <- 0

p_vals_log_seg2_sample <- hypergeo_test(log_seg2, log_seg2_sample, alpha = 1)
p_vals_log_seg4_sample <- hypergeo_test(log_seg2, log_seg4_sample, alpha = 1)

write.csv(p_vals_log_seg2_sample, "p_vals_log_seg2_sample.csv")
write.csv(p_vals_log_seg4_sample, "p_vals_log_seg4_sample.csv")