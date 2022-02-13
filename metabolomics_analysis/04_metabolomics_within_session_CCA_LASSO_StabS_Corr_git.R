#This version has no imputations - throws out data to keep things logitudinal
#Data manipulation packs
library(readxl)
library(reshape2)
library(robustbase)
library(stringr)
library(boot)
library(splitstackshape)

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
load(file=paste('ma_metabolites_eda_ts_within v0.7corr.RData',sep=''))

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

#correlation analysis
wearables_file <- 'WearablesDuringSegments-Spring2108162021.xlsx'
pwr <- read_excel(wearables_file, sheet = 'pwr')
hr <- read_excel(wearables_file, sheet = 'hr')

seg2 <- seg2[order(seg2$CLIENT_SAMPLE_ID),]
seg4 <- seg4[order(seg4$CLIENT_SAMPLE_ID),]
seg5 <- seg5[order(seg5$CLIENT_SAMPLE_ID),]

#exclude GEma048, who has janky metabolite measurements
seg2_subjects <- pwr[!is.na(pwr$seg2) & pwr$subjectID%in%seg2$CLIENT_SAMPLE_ID,]$subjectID
seg4_subjects <- pwr[!is.na(pwr$seg4) & pwr$subjectID%in%seg4$CLIENT_SAMPLE_ID,]$subjectID
seg5_subjects <- pwr[!is.na(pwr$seg5) & pwr$subjectID%in%seg5$CLIENT_SAMPLE_ID,]$subjectID

seg2_pwr <- pwr[pwr$subjectID%in%seg2_subjects,]$seg2
seg4_pwr <- pwr[pwr$subjectID%in%seg4_subjects,]$seg4
seg5_pwr <- pwr[pwr$subjectID%in%seg5_subjects,]$seg5

seg2_hr <- hr[hr$subjectID%in%seg2_subjects,]$seg2
seg4_hr <- hr[hr$subjectID%in%seg4_subjects,]$seg4
seg5_hr <- hr[hr$subjectID%in%seg5_subjects,]$seg5

seg2_bd1 <- seg2[seg2$CLIENT_SAMPLE_ID%in%seg2_subjects & seg2$CUSTOM_ATTRIBUTE_2 == "Blood Draw 1",]
seg2_bd2 <- seg2[seg2$CLIENT_SAMPLE_ID%in%seg2_subjects & seg2$CUSTOM_ATTRIBUTE_2 == "Blood Draw 2",]
seg2_bd3 <- seg2[seg2$CLIENT_SAMPLE_ID%in%seg2_subjects & seg2$CUSTOM_ATTRIBUTE_2 == "Blood Draw 3",]
seg2_bd4 <- seg2[seg2$CLIENT_SAMPLE_ID%in%seg2_subjects & seg2$CUSTOM_ATTRIBUTE_2 == "Blood Draw 4",]

seg2_bd1 <- seg2_bd1[!duplicated(seg2_bd1$CLIENT_SAMPLE_ID),]
seg2_bd2 <- seg2_bd2[!duplicated(seg2_bd2$CLIENT_SAMPLE_ID),]
seg2_bd3 <- seg2_bd3[!duplicated(seg2_bd3$CLIENT_SAMPLE_ID),]
seg2_bd4 <- seg2_bd4[!duplicated(seg2_bd4$CLIENT_SAMPLE_ID),]

seg4_bd1 <- seg4[seg4$CLIENT_SAMPLE_ID%in%seg4_subjects & seg4$CUSTOM_ATTRIBUTE_2 == "Blood Draw 1",]
seg4_bd2 <- seg4[seg4$CLIENT_SAMPLE_ID%in%seg4_subjects & seg4$CUSTOM_ATTRIBUTE_2 == "Blood Draw 2",]
seg4_bd3 <- seg4[seg4$CLIENT_SAMPLE_ID%in%seg4_subjects & seg4$CUSTOM_ATTRIBUTE_2 == "Blood Draw 3",]
seg4_bd4 <- seg4[seg4$CLIENT_SAMPLE_ID%in%seg4_subjects & seg4$CUSTOM_ATTRIBUTE_2 == "Blood Draw 4",]

seg4_bd1 <- seg4_bd1[!duplicated(seg4_bd1$CLIENT_SAMPLE_ID),]
seg4_bd2 <- seg4_bd2[!duplicated(seg4_bd2$CLIENT_SAMPLE_ID),]
seg4_bd3 <- seg4_bd3[!duplicated(seg4_bd3$CLIENT_SAMPLE_ID),]
seg4_bd4 <- seg4_bd4[!duplicated(seg4_bd4$CLIENT_SAMPLE_ID),]

seg5_bd1 <- seg5[seg5$CLIENT_SAMPLE_ID%in%seg5_subjects & seg5$CUSTOM_ATTRIBUTE_2 == "Blood Draw 1",]
seg5_bd2 <- seg5[seg5$CLIENT_SAMPLE_ID%in%seg5_subjects & seg5$CUSTOM_ATTRIBUTE_2 == "Blood Draw 2",]
seg5_bd3 <- seg5[seg5$CLIENT_SAMPLE_ID%in%seg5_subjects & seg5$CUSTOM_ATTRIBUTE_2 == "Blood Draw 3",]
seg5_bd4 <- seg5[seg5$CLIENT_SAMPLE_ID%in%seg5_subjects & seg5$CUSTOM_ATTRIBUTE_2 == "Blood Draw 4",]

seg5_bd1 <- seg5_bd1[!duplicated(seg5_bd1$CLIENT_SAMPLE_ID),]
seg5_bd2 <- seg5_bd2[!duplicated(seg5_bd2$CLIENT_SAMPLE_ID),]
seg5_bd3 <- seg5_bd3[!duplicated(seg5_bd3$CLIENT_SAMPLE_ID),]
seg5_bd4 <- seg5_bd4[!duplicated(seg5_bd4$CLIENT_SAMPLE_ID),]

#check for subject alignment
all.equal(seg2_subjects, seg2_bd1$CLIENT_SAMPLE_ID)
all.equal(seg2_subjects, seg2_bd2$CLIENT_SAMPLE_ID)
all.equal(seg2_subjects_2, seg2_bd3$CLIENT_SAMPLE_ID)
all.equal(seg2_subjects_2, seg2_bd4$CLIENT_SAMPLE_ID)

all.equal(seg4_subjects, seg4_bd1$CLIENT_SAMPLE_ID)
all.equal(seg4_subjects, seg4_bd2$CLIENT_SAMPLE_ID)
all.equal(seg4_subjects, seg4_bd3$CLIENT_SAMPLE_ID)
all.equal(seg4_subjects, seg4_bd4$CLIENT_SAMPLE_ID)

all.equal(seg5_subjects, seg5_bd1$CLIENT_SAMPLE_ID)
all.equal(seg5_subjects, seg5_bd2$CLIENT_SAMPLE_ID)
all.equal(seg5_subjects, seg5_bd3$CLIENT_SAMPLE_ID)
all.equal(seg5_subjects, seg5_bd4$CLIENT_SAMPLE_ID)

#log2fc
seg2_l21 <- log(exp(seg2_bd2[,2:1270] - seg2_bd1[,2:1270]),2)
seg2_l32 <- log(exp(seg2_bd3[,2:1270] - seg2_bd2[,2:1270]),2)
seg2_l42 <- log(exp(seg2_bd4[,2:1270] - seg2_bd2[,2:1270]),2)

seg4_l21 <- log(exp(seg4_bd2[,2:1270] - seg4_bd1[,2:1270]),2)
seg4_l32 <- log(exp(seg4_bd3[,2:1270] - seg4_bd2[,2:1270]),2)
seg4_l42 <- log(exp(seg4_bd4[,2:1270] - seg4_bd2[,2:1270]),2)

seg5_l21 <- log(exp(seg5_bd2[,2:1270] - seg5_bd1[,2:1270]),2)
seg5_l32 <- log(exp(seg5_bd3[,2:1270] - seg5_bd2[,2:1270]),2)
seg5_l42 <- log(exp(seg5_bd4[,2:1270] - seg5_bd2[,2:1270]),2)

seg2_bd1 <- seg2_bd1[,2:1270]
seg2_bd2 <- seg2_bd2[,2:1270]
seg2_bd3 <- seg2_bd3[,2:1270]
seg2_bd4 <- seg2_bd4[,2:1270]

seg4_bd1 <- seg4_bd1[,2:1270]
seg4_bd2 <- seg4_bd2[,2:1270]
seg4_bd3 <- seg4_bd3[,2:1270]
seg4_bd4 <- seg4_bd4[,2:1270]

seg5_bd1 <- seg5_bd1[,2:1270]
seg5_bd2 <- seg5_bd2[,2:1270]
seg5_bd3 <- seg5_bd3[,2:1270]
seg5_bd4 <- seg5_bd4[,2:1270]

####normalize but only for regression
seg2_bd1_n <- standardize(seg2_bd1)
seg2_bd2_n <- standardize(seg2_bd2)
seg2_bd3_n <- standardize(seg2_bd3)
seg2_bd4_n <- standardize(seg2_bd4)

seg4_bd1_n <- standardize(seg4_bd1)
seg4_bd2_n <- standardize(seg4_bd2)
seg4_bd3_n <- standardize(seg4_bd3)
seg4_bd4_n <- standardize(seg4_bd4)

seg5_bd1_n <- standardize(seg5_bd1)
seg5_bd2_n <- standardize(seg5_bd2)
seg5_bd3_n <- standardize(seg5_bd3)
seg5_bd4_n <- standardize(seg5_bd4)

seg2_bd1_n[is.na(seg2_bd1_n)] <- 0
seg2_bd2_n[is.na(seg2_bd2_n)] <- 0
seg2_bd3_n[is.na(seg2_bd3_n)] <- 0
seg2_bd4_n[is.na(seg2_bd4_n)] <- 0

seg4_bd1_n[is.na(seg4_bd1_n)] <- 0
seg4_bd2_n[is.na(seg4_bd2_n)] <- 0
seg4_bd3_n[is.na(seg4_bd3_n)] <- 0
seg4_bd4_n[is.na(seg4_bd4_n)] <- 0

seg5_bd1_n[is.na(seg5_bd1_n)] <- 0
seg5_bd2_n[is.na(seg5_bd2_n)] <- 0
seg5_bd3_n[is.na(seg5_bd3_n)] <- 0
seg5_bd4_n[is.na(seg5_bd4_n)] <- 0

seg2_l21_n <- standardize(seg2_l21)
seg2_l32_n <- standardize(seg2_l32)
seg2_l42_n <- standardize(seg2_l42)

seg4_l21_n <- standardize(seg4_l21)
seg4_l32_n <- standardize(seg4_l32)
seg4_l42_n <- standardize(seg4_l42)

seg5_l21_n <- standardize(seg5_l21)
seg5_l32_n <- standardize(seg5_l32)
seg5_l42_n <- standardize(seg5_l42)

seg2_l21_n[is.na(seg2_l21_n)] <- 0
seg2_l32_n[is.na(seg2_l32_n)] <- 0
seg2_l42_n[is.na(seg2_l42_n)] <- 0

seg4_l21_n[is.na(seg4_l21_n)] <- 0
seg4_l32_n[is.na(seg4_l32_n)] <- 0
seg4_l42_n[is.na(seg4_l42_n)] <- 0

seg5_l21_n[is.na(seg5_l21_n)] <- 0
seg5_l32_n[is.na(seg5_l32_n)] <- 0
seg5_l42_n[is.na(seg5_l42_n)] <- 0

log_prot_sans_na1 <- read.csv('log_prot_S2BD1.csv',header=T)[,-1]
log_prot_sans_na2 <- read.csv('log_prot_S2BD2.csv',header=T)[,-1]
log_prot_sans_na3 <- read.csv('log_prot_S2BD3.csv',header=T)[,-1]

seg2_subjects_prot_bd1 <- read.csv('log_prot_S2BD1.csv',header=T)[,1]
seg2_subjects_prot_bd2 <- read.csv('log_prot_S2BD2.csv',header=T)[,1]
seg2_subjects_prot_bd3 <- read.csv('log_prot_S2BD3.csv',header=T)[,1]

seg2_bd1_prot <- log_prot_sans_na1[seg2_subjects_prot_bd1%in%seg2_subjects,]
seg2_bd2_prot <- log_prot_sans_na2[seg2_subjects_prot_bd2%in%seg2_subjects,]
seg2_bd3_prot <- log_prot_sans_na3[seg2_subjects_prot_bd3%in%seg2_subjects,]

temp <- colnames(seg2_bd1_prot)[colnames(seg2_bd1_prot)%in%colnames(seg2_bd2_prot)]
shared_cols <- temp[temp%in%colnames(seg2_bd3_prot)]

seg2_bd1_prot_c <- seg2_bd1_prot[,shared_cols]
seg2_bd2_prot_c <- seg2_bd2_prot[,shared_cols]
seg2_bd3_prot_c <- seg2_bd3_prot[,shared_cols]

temp <- seg2_subjects_prot_bd1[seg2_subjects_prot_bd1%in%seg2_subjects]
temp2 <- seg2_subjects_prot_bd2[seg2_subjects_prot_bd2%in%seg2_subjects]
temp3 <- seg2_subjects_prot_bd3[seg2_subjects_prot_bd3%in%seg2_subjects]

seg2_bd1_prot_cc <- seg2_bd1_prot[temp%in%temp2,shared_cols]
seg2_bd2_prot_cc <- seg2_bd2_prot[temp2%in%temp2,shared_cols]
seg2_bd3_prot_cc <- seg2_bd3_prot[temp3%in%temp2,shared_cols]

seg2_bd1_np <- cbind(seg2_bd1_n,seg2_bd1_prot)
seg2_bd2_np <- cbind(seg2_bd2_n[seg2_subjects%in%seg2_subjects_prot_bd2,],seg2_bd2_prot)
seg2_bd3_np <- cbind(seg2_bd3_n,seg2_bd3_prot)

seg2_bd1_np <- seg2_bd1_np[seg2_subjects_prot_bd1[(seg2_subjects_prot_bd1%in%seg2_subjects)]%in%seg2_subjects_prot_bd2,]
seg2_bd3_np <- seg2_bd3_np[seg2_subjects_prot_bd3[(seg2_subjects_prot_bd3%in%seg2_subjects)]%in%seg2_subjects_prot_bd2,]

seg2_bd1_npc <- cbind(seg2_bd1_n,seg2_bd1_prot_c)
seg2_bd2_npc <- cbind(seg2_bd2_n[seg2_subjects%in%seg2_subjects_prot_bd2,],seg2_bd2_prot_c)
seg2_bd3_npc <- cbind(seg2_bd3_n,seg2_bd3_prot_c)

seg2_bd1_npc <- seg2_bd1_npc[seg2_subjects_prot_bd1[(seg2_subjects_prot_bd1%in%seg2_subjects)]%in%seg2_subjects_prot_bd2,]
seg2_bd3_npc <- seg2_bd3_npc[seg2_subjects_prot_bd3[(seg2_subjects_prot_bd3%in%seg2_subjects)]%in%seg2_subjects_prot_bd2,]

seg2_l21_npc <- log(exp(seg2_bd2_npc-seg2_bd1_npc),2)
seg2_l32_npc <- log(exp(seg2_bd3_npc-seg2_bd2_npc),2)

#A quick side excursion to random forest and how it performs on 
#regression tasks for the given dataset
library('randomForest')
library('caTools')

input_list <- list(seg2_bd2_n, seg4_bd2_n, seg5_bd2_n)
import_list <- matrix(NA,300,length(input_list))
pr2_list <- c()

k = 1
tgt <- as.numeric(as.numeric(vo2$`Percent Body Fat (%)`))
for (each in input_list){
  each_t<- each
  n <- (1:(nrow(each_t)-1))
  train <- sample(n,floor(0.8*length(n)))
  test <- (n)[!(n)%in%train]
  
  train_set <- cbind(tgt[train], each_t[train,])
  colnames(train_set)[1] <- "target"
  rf <- randomForest(target~.,data=data.frame(train_set))
  pred <- predict(rf, newdata=data.frame(each_t[test,]))
  pr2 <- 1 - mean((tgt[test] - pred)^2)/var(tgt[test])
  pr2_list <- c(pr2_list,pr2)
  import <- rf$importance[order(rf$importance[,1], decreasing = T),]
  import_list[,k] <- sapply(names(import[1:300]), function(x){substr(x,2,nchar(x))})
  k = k + 1
}

out <- Reduce(intersect, list(import_list[,1], import_list[,2], import_list[,3]))
(out <- c(out, pr2_list[1:3]))
write.csv(out, file = paste(getwd(),"/metab_vo2bf_bd2.csv",sep=""))

placeholder <- vo2[vo2$`subject ID`%in%seg2_subjects_prot_bd2,]$`Sandhurst Spring 2021 Participation on Record`
placeholder <- placeholder[-17]
placeholder[placeholder == "N"] <- 0
placeholder[placeholder!="0"] <- 1
placeholder <- as.factor(placeholder)

data_list <- list(
  seg2_pwr,
  seg2_hr,
  as.numeric(vo2$`VO2 Max (mL/Kg/min)`),
  as.numeric(vo2$`Power @ VO2max (From Fitmate)- W`),
  as.numeric(vo2$`Percent Body Fat (%)`),
  as.numeric(vo2$`Skeletal Muscle Mass (lbs)`),
  as.factor(vo2$`Sandhurst Spring 2021 Participation on Record`)
)

ol_subject <- seg2_subjects[seg2_subjects%in%seg2_subjects_prot_bd2][8] #Possible outlier
post_ol_idx <- seg2_subjects[seg2_subjects%in%seg2_subjects_prot_bd2]!=ol_subject

##CCA starts here
library(whitening)
library(PMA)

elim_no_var <- function(x){
#For eliminating featues without any variance
  y <- apply(x,2,var)
  y <- t(matrix(rep(as.logical(y != 0),nrow(x)),ncol=nrow(x)))
  x <- x*y
  x <- x[,colSums(x != 0) != 0]
  x
}


feats <- read.csv('temp_feat.csv')

#Preprocessing
x1 = seg2_bd2_np
x2 = cbind(seg2_pwr[seg2_subjects%in%seg2_subjects_prot_bd2],
           seg2_hr[seg2_subjects%in%seg2_subjects_prot_bd2],
  as.numeric(vo2$`VO2 Max (mL/Kg/min)`[seg2_subjects%in%seg2_subjects_prot_bd2]), #seg2_subjects%in%seg2_subjects_prot_bd2
  as.numeric(vo2$`Power @ VO2max (From Fitmate)- W`[seg2_subjects%in%seg2_subjects_prot_bd2]),
  as.numeric(vo2$`Percent Body Fat (%)`[seg2_subjects%in%seg2_subjects_prot_bd2]),
           as.numeric(vo2$`Skeletal Muscle Mass (lbs)`[seg2_subjects%in%seg2_subjects_prot_bd2])
)

check <- apply(x2,1,sum)
leave <- -(1:length(check))[is.na(check)]
if (length(leave)!=0){
  x2 <- x2[leave,]
  x1 <- x1[leave,]
}
x1 <- elim_no_var(x1)
x1 <- scale(x1)
x2 <- scale(x2)
colnames(x2) <- c("seg2pwr", "seg2hr", "vo2volume", 
                  "vo2power", "vo2bodyFat", "vo2muscleMass")

colsize <- ncol(x1)
set.seed(123)
x1m <- cbind(x1, matrix(rnorm(n=nrow(x1) * 50 ,sd=1), nrow = nrow(x1)))
x1m <- x1m[rownames(x1m)!='8',]
x2m <- x2[-7,]
x1m=scale(x1m)
x2m=scale(x2m)

var_select <- matrix(NA, nrow=750,ncol=100)
tar_select <- matrix(NA, nrow=6,ncol=100)
var_loadings <- matrix(NA, nrow=750,ncol=100)
tar_loadings <- matrix(NA, nrow=6,ncol=100)
corz <- c()

for(i in 1:100){
  x1m <- cbind(x1, matrix(rnorm(n=nrow(x1) * 50 ,sd=1), nrow = nrow(x1)))
  x1m <- x1m[rownames(x1m)!='8',]
  x1m=scale(x1m)
  #CCA, unsupervised, starts here
  perm.out1 <- CCA.permute(x=x1m,z=x2m,typex="standard",
                          typez="standard",nperms=250, 
                          penaltyxs=seq(0,0.475,len=10),
                          penaltyzs=seq(0,0.475,len=10))
  
  out1 <- CCA(x=x1m,z=x2m, typex="standard", typez="standard",
              v=perm.out1$v.init,
              penaltyx=perm.out1$bestpenaltyx,
              penaltyz=perm.out1$bestpenaltyz)
  
  # Identifying features that have non-zero weights
  to_csv1 <- cbind(colnames(x1m_sub)[out1$u!=0],out1$u[out1$u!=0],colnames(x2m_sub)[out1$v!=0],out1$v[out1$v!=0], out1$cors)
  to_csv1[-(1:length(colnames(x2m_sub)[out1$v!=0])),3] <- ""
  to_csv1[-1,5] <- ""
  to_csv1[-(1:sum(out1$u!=0)),2] <- ""
  to_csv1[-(1:sum(out1$v!=0)),4] <- ""
  var_select[1:length(colnames(x1m_sub)[out1$u!=0]),i] <- colnames(x1m_sub)[out1$u!=0]
  tar_select[1:length(colnames(x2m_sub)[out1$v!=0]),i] <- colnames(x2m_sub)[out1$v!=0]
  var_loadings[1:length(colnames(x1m_sub)[out1$u!=0]),i] <- out1$u[out1$u!=0]
  tar_loadings[1:length(colnames(x2m_sub)[out1$v!=0]),i] <- out1$v[out1$v!=0]
  corz <- c(corz,out1$cors)
}

colnames(to_csv1) <- c("feat", "feat_loading", "target", "target_loading", "corr")
save.image(file=paste('ma_metabolites_eda_ts_within v0.82corr_sCCA_testing.RData',sep=''))

plot_noise <- apply(var_select=="",2,sum)

names(plot_noise) <- 1:100

plot_selection <- apply(var_select!="" & var_select != "NA",2,sum)
names(plot_selection) <- 1:100

barplot(plot_noise, xlab="Run", ylab="# Noisy Variables",ylim=c(0,25),
        cex.lab=1.25, main = "# of Noisy Variables Selected Out of 50 vs. Selection Run", cex.main=1.5)

barplot(plot_selection, xlab="Run", ylab="# Features",ylim=c(0,550),
        cex.lab=1.25, main = "# of Features Selected Out of 1284 vs. Selection Run", cex.main=1.5)

#CCA, OG
set.seed(123) 
set.seed(23) 
perm.out <- CCA.permute(x=x1m,z=x2m,typex="standard",
                         typez="standard",nperms=250, 
                         penaltyxs=seq(0,0.475,len=10),
                         penaltyzs=seq(0,0.475,len=10))

out <- CCA(x=x1m,z=x2m, typex="standard", typez="standard",
            v=perm.out$v.init,
            penaltyx= perm.out$bestpenaltyx,
            penaltyz= perm.out$bestpenaltyz) 

print(perm.out$bestpenaltyx)
print(perm.out$bestpenaltyz)
# Identifying features that have non-zero weights
(to_csv <- cbind(colnames(x1)[out$u!=0],out$u[out$u!=0],colnames(x2)[out$v!=0],out$v[out$v!=0],
                 out$cors))
sum(is.na(to_csv[,1]))
to_csv[-(1:length(colnames(x2)[out$v!=0])),3] <- ""
to_csv[-1,5] <- ""
to_csv[-(1:sum(out$u!=0)),2] <- ""
to_csv[-(1:sum(out$v!=0)),4] <- ""
colnames(to_csv) <- c("feat", "feat_loading", "target", "target_loading", "corr")
write.csv(to_csv,file=paste(getwd(),"/cca_test_low_0.1_corr1.csv",sep=""))

###Try LASSO stability selection for dimensionality reduction
library('stabs')
library('lars')
library('glmnet')
data("bodyfat", package = "TH.data")

vo2 <- read_excel(paste(path,'VO2max-InBody-Spring2108142021.xlsx',sep=''))
vo2 <- vo2[vo2$`subject ID`%in%seg2_subjects,]

#Sandhurst participation, vo2max volume
vo2$`VO2 Max (mL/Kg/min)`
vo2$`Power @ VO2max (From Fitmate)- W`
vo2$`Percent Body Fat (%)`
vo2$`Skeletal Muscle Mass (lbs)`
vo2$`Sandhurst Spring 2021 Participation on Record`

placeholder <- vo2[vo2$`subject ID`%in%seg2_subjects_prot_bd2,]$`Sandhurst Spring 2021 Participation on Record`
placeholder <- placeholder[-17]
placeholder[placeholder == "N"] <- 0
placeholder[placeholder!="0"] <- 1
placeholder <- as.factor(placeholder)

set.seed(123)
start_time <- Sys.time()
seg5_pwr_bd4_stab <- stabsel(seg2_bd1_np[-38,], 
                             as.numeric(vo2$`Power @ VO2max (From Fitmate)- W`)[-38],
        cutoff = 0.6, PFER = 1, fitfun = glmnet.lasso)

seg5_pwr_bd4_stab <- stabsel(seg2_l21_npc, 
                             as.numeric(vo2$`Skeletal Muscle Mass (lbs)`
                                        [vo2$`subject ID`%in%seg2_subjects_prot_bd2]),
                             cutoff = 0.6, PFER = 1, fitfun = glmnet.lasso)

seg5_pwr_bd4_stab <- stabsel(seg2_bd3_npc[-3,], 
                             as.numeric(vo2$`HR (Fitmate) @ VO2 max- bpm`
                                        [seg2_subjects%in%seg2_subjects_prot_bd2])[-3],
                             cutoff = 0.6, PFER = 1, fitfun = glmnet.lasso)

seg5_pwr_bd4_stab <- stabsel(seg2_bd3_prot_cc, 
                             seg2_hr[seg2_subjects%in%temp2],
                             cutoff = 0.6, PFER = 1, fitfun = glmnet.lasso)


end_time <- Sys.time()
end_time - start_time

seg5_pwr_bd4_stab$selected

plot(seg5_pwr_bd4_stab, main = "LASSO Stability Selection for Stress Test Heart Rate vs. BD2 (Segment 2)")

mod <- lm(seg2_hr~seg2_bd4_n[,'272'])
summary(mod)

library(superpc)

set.seed(123)

# create train and test data objects. censoring.status=1 means the event occurred;
# censoring.status=0 means censored

is_train <- sample(69,floor(0.8*69))
is_not_train <- (1:69)[!(1:69)%in%is_train]
is_val <- sample(is_not_train, floor(length(is_not_train)/2))
is_test <- is_not_train[!is_not_train%in%is_val]
featurenames <- colnames(elim_no_var(seg2_bd2_npc))

data <- list(x=t(elim_no_var(seg2_bd3_prot_cc)),
             y=as.numeric(as.numeric(vo2$`Skeletal Muscle Mass (lbs)`
                                     [vo2$`subject ID`%in%seg2_subjects_prot_bd2])
                          ), featurenames=featurenames)
data_train<-list(x=t(elim_no_var(seg2_bd2_n[is_train,])),
           y=as.numeric(vo2$`VO2 Max (mL/Kg/min)`[is_train]), featurenames=featurenames)
data_val <- list(x=t(elim_no_var(seg2_bd2_n[is_val,])),
                 y=as.numeric(vo2$`VO2 Max (mL/Kg/min)`[is_val]), featurenames=featurenames)

data_test <- list(x=t(elim_no_var(seg2_bd2_n[is_test,])),
           y=as.numeric(vo2$`VO2 Max (mL/Kg/min)`[is_test]), featurenames=featurenames)

data <- list(x=t(elim_no_var(seg2_bd3_prot_cc)[-17,]),
             y=as.logical(as.numeric(placeholder)-1), featurenames=featurenames)
# train  the model. This step just computes the  scores for each feature
train.obj<- superpc.train(data, type="regression", s0.perc = 0.5)
cv.obj<-superpc.cv(train.obj, data)
superpc.plotcv(cv.obj)
# here we have the luxury of test data, so we can compute the  likelihood ratio statistic
# over the test data and plot them. We see that the threshold of 0.7
# works pretty well
lrtest.obj<-superpc.lrtest.curv(train.obj, data,data_val)
lrtest.obj$threshold
superpc.plot.lrtest(lrtest.obj)

# now we derive the predictor of survival  for the test data, 
# and then then use it
# as the predictor in a Cox model . We see that the 1st supervised PC is
# highly significant; the next two are not
fit.cts<- superpc.predict(train.obj, data, data_val, threshold=0.5, n.components=3, 
                          prediction.type="continuous")
superpc.fit.to.outcome(train.obj, data_val, fit.cts$v.pred)

fit.red<- superpc.predict.red(train.obj, data, data_val, threshold=0.5)

fit.redcv<- superpc.predict.red.cv(fit.red, cv.obj,  data,  threshold=0.5)

superpc.plotred.lrtest(fit.redcv)

##Try default dataset####
data("bodyfat", package = "TH.data")
## set seed
set.seed(1234)

## lasso
(stab.lasso <- stabsel(x = bodyfat[, -2], y = bodyfat[,2],
                       fitfun = glmnet.lasso, cutoff = 0.9,
                       PFER = 1))

par(mfrow = c(1, 1))
plot(stab.lasso, main = "Lasso")
#create synthetic dataset using iris to see if the pkg is working correctly#
data('iris')
test_in <- matrix(NA, nrow = 150, ncol = 1269)
test_in[,1:3] <- as.matrix(iris[,2:4])
test_in[,1:3] <- standardize(test_in[,1:3])
test_in[,5] <- as.factor(iris[,5])
test_out <-  as.matrix(iris[,1])
test_in[,6:1269] <- rnorm(1264*150)

(stab_iris <- stabsel(x = test_in, y = test_out,
                       fitfun = glmnet.lasso, cutoff = 0.8,
                       PFER = 1))
plot(stab_iris, main = "Lasso")

#Here, check several pearson correlations vs. targetes to calculate prize values
#for omics integrator.
seg2_x_list <- list(seg2_bd1, seg2_bd2, seg2_bd3, seg2_bd4, seg2_l21, seg2_l32, seg2_l42)
seg4_x_list <- list(seg4_bd1, seg4_bd2, seg4_bd3, seg4_bd4, seg4_l21, seg4_l32, seg4_l42)
seg5_x_list <- list(seg5_bd1, seg5_bd2, seg5_bd3, seg5_bd4, seg5_l21, seg5_l32, seg5_l42)

seg2_pwr_pcor_list <- list()
seg4_pwr_pcor_list <- list()
seg5_pwr_pcor_list <- list()
seg2_hr_pcor_list <- list()
seg4_hr_pcor_list <- list()
seg5_hr_pcor_list <- list()

seg2_pwr_scor_list <- list()
seg4_pwr_scor_list <- list()
seg5_pwr_scor_list <- list()
seg2_hr_scor_list <- list()
seg4_hr_scor_list <- list()
seg5_hr_scor_list <- list()

for (i in 1:length(seg2_x_list)){
  seg2_pwr_pcor_list <- 
    append(seg2_pwr_pcor_list, apply(seg2_x_list[[i]],2,cor,y=seg2_pwr,method="pearson"))
}
for (i in 1:length(seg4_x_list)){
  seg4_pwr_pcor_list <- 
    append(seg4_pwr_pcor_list, apply(seg4_x_list[[i]],2,cor,y=seg4_pwr,method="pearson"))
}
for (i in 1:length(seg5_x_list)){
  seg5_pwr_pcor_list <- 
    append(seg5_pwr_pcor_list, apply(seg5_x_list[[i]],2,cor,y=seg5_pwr,method="pearson"))
}

for (i in 1:length(seg2_x_list)){
  seg2_hr_pcor_list <- 
    append(seg2_hr_pcor_list, apply(seg2_x_list[[i]],2,cor,y=seg2_hr,method="pearson"))
}
for (i in 1:length(seg4_x_list)){
  seg4_hr_pcor_list <- 
    append(seg4_hr_pcor_list, apply(seg4_x_list[[i]],2,cor,y=seg4_hr,method="pearson"))
}
for (i in 1:length(seg5_x_list)){
  seg5_hr_pcor_list <- 
    append(seg5_hr_pcor_list, apply(seg5_x_list[[i]],2,cor,y=seg5_hr,method="pearson"))
}

seg2_pwr_scor_list <- list()
seg4_pwr_scor_list <- list()
seg5_pwr_scor_list <- list()
seg2_hr_scor_list <- list()
seg4_hr_scor_list <- list()
seg5_hr_scor_list <- list()

seg2_pwr_scor_list <- list()
seg4_pwr_scor_list <- list()
seg5_pwr_scor_list <- list()
seg2_hr_scor_list <- list()
seg4_hr_scor_list <- list()
seg5_hr_scor_list <- list()

for (i in 1:length(seg2_x_list)){
  seg2_pwr_scor_list <- 
    append(seg2_pwr_scor_list, apply(seg2_x_list[[i]],2,cor,y=seg2_pwr,method="spearman"))
}
for (i in 1:length(seg4_x_list)){
  seg4_pwr_scor_list <- 
    append(seg4_pwr_scor_list, apply(seg4_x_list[[i]],2,cor,y=seg4_pwr,method="spearman"))
}
for (i in 1:length(seg5_x_list)){
  seg5_pwr_scor_list <- 
    append(seg5_pwr_scor_list, apply(seg5_x_list[[i]],2,cor,y=seg5_pwr,method="spearman"))
}

for (i in 1:length(seg2_x_list)){
  seg2_hr_scor_list <- 
    append(seg2_hr_scor_list, apply(seg2_x_list[[i]],2,cor,y=seg2_hr,method="spearman"))
}
for (i in 1:length(seg4_x_list)){
  seg4_hr_scor_list <- 
    append(seg4_hr_scor_list, apply(seg4_x_list[[i]],2,cor,y=seg4_hr,method="spearman"))
}
for (i in 1:length(seg5_x_list)){
  seg5_hr_scor_list <- 
    append(seg5_hr_scor_list, apply(seg5_x_list[[i]],2,cor,y=seg5_hr,method="spearman"))
}

seg2_pwr_pcor <- matrix(seg2_pwr_pcor_list,nrow=1269)
seg4_pwr_pcor <- matrix(seg4_pwr_pcor_list,nrow=1269)
seg5_pwr_pcor <- matrix(seg5_pwr_pcor_list,nrow=1269)

seg2_hr_pcor <- matrix(seg2_hr_pcor_list,nrow=1269)
seg4_hr_pcor <- matrix(seg4_hr_pcor_list,nrow=1269)
seg5_hr_pcor <- matrix(seg5_hr_pcor_list,nrow=1269)

seg2_pwr_scor <- matrix(seg2_pwr_scor_list,nrow=1269)
seg4_pwr_scor <- matrix(seg4_pwr_scor_list,nrow=1269)
seg5_pwr_scor <- matrix(seg5_pwr_scor_list,nrow=1269)

seg2_hr_scor <- matrix(seg2_hr_scor_list,nrow=1269)
seg4_hr_scor <- matrix(seg4_hr_scor_list,nrow=1269)
seg5_hr_scor <- matrix(seg5_hr_scor_list,nrow=1269)
  
rownames(seg2_pwr_pcor) <- as.numeric(colnames(dat[,2:1270]))
rownames(seg4_pwr_pcor) <- as.numeric(colnames(dat[,2:1270]))
rownames(seg5_pwr_pcor) <- as.numeric(colnames(dat[,2:1270]))

rownames(seg2_hr_pcor) <- as.numeric(colnames(dat[,2:1270]))
rownames(seg4_hr_pcor) <- as.numeric(colnames(dat[,2:1270]))
rownames(seg5_hr_pcor) <- as.numeric(colnames(dat[,2:1270]))

rownames(seg2_pwr_scor) <- as.numeric(colnames(dat[,2:1270]))
rownames(seg4_pwr_scor) <- as.numeric(colnames(dat[,2:1270]))
rownames(seg5_pwr_scor) <- as.numeric(colnames(dat[,2:1270]))

rownames(seg2_hr_scor) <- as.numeric(colnames(dat[,2:1270]))
rownames(seg4_hr_scor) <- as.numeric(colnames(dat[,2:1270]))
rownames(seg5_hr_scor) <- as.numeric(colnames(dat[,2:1270]))

cols <- c("BD1", "BD2", "BD3", "BD4", "log2FC21", "log2FC32", "log2FC42")

colnames(seg2_pwr_pcor) <- cols
colnames(seg4_pwr_pcor) <- cols
colnames(seg5_pwr_pcor) <- cols

colnames(seg2_hr_pcor) <- cols
colnames(seg4_hr_pcor) <- cols
colnames(seg5_hr_pcor) <- cols

colnames(seg2_pwr_scor) <- cols
colnames(seg4_pwr_scor) <- cols
colnames(seg5_pwr_scor) <- cols

colnames(seg2_hr_scor) <- cols
colnames(seg4_hr_scor) <- cols
colnames(seg5_hr_scor) <- cols

seg2_pwr_pcor <- merge(seg2_pwr_pcor,metab_info, by.x =0, by.y = "CHEM_ID")[,1:9]
seg4_pwr_pcor <- merge(seg4_pwr_pcor,metab_info, by.x =0, by.y = "CHEM_ID")[,1:9]
seg5_pwr_pcor <- merge(seg5_pwr_pcor,metab_info, by.x =0, by.y = "CHEM_ID")[,1:9]
colnames(seg2_pwr_pcor)[1] <- "CHEM_ID"
colnames(seg4_pwr_pcor)[1] <- "CHEM_ID"
colnames(seg5_pwr_pcor)[1] <- "CHEM_ID"
seg2_pwr_pcor <- seg2_pwr_pcor[,c(1, 9, 2:8)]
seg4_pwr_pcor <- seg4_pwr_pcor[,c(1, 9, 2:8)]
seg5_pwr_pcor <- seg5_pwr_pcor[,c(1, 9, 2:8)]

seg2_hr_pcor <- merge(seg2_hr_pcor,metab_info, by.x =0, by.y = "CHEM_ID")[,1:9]
seg4_hr_pcor <- merge(seg4_hr_pcor,metab_info, by.x =0, by.y = "CHEM_ID")[,1:9]
seg5_hr_pcor <- merge(seg5_hr_pcor,metab_info, by.x =0, by.y = "CHEM_ID")[,1:9]
colnames(seg2_hr_pcor)[1] <- "CHEM_ID"
colnames(seg4_hr_pcor)[1] <- "CHEM_ID"
colnames(seg5_hr_pcor)[1] <- "CHEM_ID"
seg2_hr_pcor <- seg2_hr_pcor[,c(1, 9, 2:8)]
seg4_hr_pcor <- seg4_hr_pcor[,c(1, 9, 2:8)]
seg5_hr_pcor <- seg5_hr_pcor[,c(1, 9, 2:8)]

seg2_pwr_pcor <- apply(seg2_pwr_pcor,2,as.character)
seg4_pwr_pcor <- apply(seg4_pwr_pcor,2,as.character)
seg5_pwr_pcor <- apply(seg5_pwr_pcor,2,as.character)

seg2_hr_pcor <- apply(seg2_hr_pcor,2,as.character)
seg4_hr_pcor <- apply(seg4_hr_pcor,2,as.character)
seg5_hr_pcor <- apply(seg5_hr_pcor,2,as.character)

write.csv(seg2_pwr_pcor, 'seg2_pwr_pcor.csv')
write.csv(seg4_pwr_pcor, 'seg4_pwr_pcor.csv')
write.csv(seg5_pwr_pcor, 'seg5_pwr_pcor.csv')

write.csv(seg2_hr_pcor, 'seg2_hr_pcor.csv')
write.csv(seg4_hr_pcor, 'seg4_hr_pcor.csv')
write.csv(seg5_hr_pcor, 'seg5_hr_pcor.csv')

write.csv(seg2_pwr_scor)
write.csv(seg4_pwr_scor)
write.csv(seg5_pwr_scor)

write.csv(seg2_hr_scor)
write.csv(seg4_hr_scor)
write.csv(seg5_hr_scor)

#####end of corr######