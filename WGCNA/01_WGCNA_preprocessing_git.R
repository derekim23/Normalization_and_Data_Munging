#This version has no imputations - throws out data to keep things logitudinal
#Data manipulation packs
library(readxl)
library(reshape2)
library(robustbase)
library(stringr)
library(boot)
library(splitstackshape)

#Bioconductor packs
library(Biobase)
library(convert)

#WGCNA
library(WGCNA)

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
  
  #Take the intersection of cadet ids present across all sessions.
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

calc_dx <- function(df){
  #Calculate the delta between consecutive time steps for the seg(k)_time dfs
  #Impute missing time diffs with medians between two "same" time steps
  dx <- as.matrix(apply(df[,4:ncol(df)],2,as.numeric)) - 
    as.matrix(apply(df[,3:(ncol(df)-1)],2,as.numeric))
  med <- colMedians(dx, na.rm = T)
  idx_na <- is.na(dx)
  dx[idx_na] <- med[col(dx)][idx_na]
  dx[is.na(df[,2]),] <- NA
  dx <- as.data.frame(cbind(df[,1]$`Cadet ID`,as.data.frame(dx)))
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

#Load main dataframe
path <- "/Volumes/GoogleDrive-115111199924997198421/My Drive/Objectives/Code/"
setwd(path)
file <- 'Preprocessed_DARPA-MBA-UntargetedMetabolomics-ReNormalizedData-Spring21WestPoint-09272021-AA-DK.xlsx'
dat <- read_excel(file, sheet = 'Log-Trans. and Merged')
time_file <- 'DARPA-MBA-Phlebotomy-Spring21WestPoint-08172021-AA.xlsx'

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

#Sort by cadet ID, blood draw, then segments
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
wearables_file <- 'DARPA-MBA-WearablesDuringSegments-Spring21WestPoint-08162021-AA.xlsx'
pwr <- read_excel(wearables_file, sheet = 'pwr')
hr <- read_excel(wearables_file, sheet = 'hr')

seg2 <- seg2[order(seg2$CLIENT_SAMPLE_ID),]
seg4 <- seg4[order(seg4$CLIENT_SAMPLE_ID),]
seg5 <- seg5[order(seg5$CLIENT_SAMPLE_ID),]

#exclude GEWP048, who has janky metabolite measurements
seg2_cadets <- pwr[!is.na(pwr$seg2) & pwr$CadetID%in%seg2$CLIENT_SAMPLE_ID,]$CadetID
seg4_cadets <- pwr[!is.na(pwr$seg4) & pwr$CadetID%in%seg4$CLIENT_SAMPLE_ID,]$CadetID
seg5_cadets <- pwr[!is.na(pwr$seg5) & pwr$CadetID%in%seg5$CLIENT_SAMPLE_ID,]$CadetID

seg2_pwr <- pwr[pwr$CadetID%in%seg2_cadets,]$seg2
seg4_pwr <- pwr[pwr$CadetID%in%seg4_cadets,]$seg4
seg5_pwr <- pwr[pwr$CadetID%in%seg5_cadets,]$seg5

seg2_hr <- hr[hr$CadetID%in%seg2_cadets,]$seg2
seg4_hr <- hr[hr$CadetID%in%seg4_cadets,]$seg4
seg5_hr <- hr[hr$CadetID%in%seg5_cadets,]$seg5

seg2_bd1 <- seg2[seg2$CLIENT_SAMPLE_ID%in%seg2_cadets & seg2$CUSTOM_ATTRIBUTE_2 == "Blood Draw 1",]
seg2_bd2 <- seg2[seg2$CLIENT_SAMPLE_ID%in%seg2_cadets & seg2$CUSTOM_ATTRIBUTE_2 == "Blood Draw 2",]
seg2_bd3 <- seg2[seg2$CLIENT_SAMPLE_ID%in%seg2_cadets & seg2$CUSTOM_ATTRIBUTE_2 == "Blood Draw 3",]
seg2_bd4 <- seg2[seg2$CLIENT_SAMPLE_ID%in%seg2_cadets & seg2$CUSTOM_ATTRIBUTE_2 == "Blood Draw 4",]

seg2_bd1 <- seg2_bd1[!duplicated(seg2_bd1$CLIENT_SAMPLE_ID),]
seg2_bd2 <- seg2_bd2[!duplicated(seg2_bd2$CLIENT_SAMPLE_ID),]
seg2_bd3 <- seg2_bd3[!duplicated(seg2_bd3$CLIENT_SAMPLE_ID),]
seg2_bd4 <- seg2_bd4[!duplicated(seg2_bd4$CLIENT_SAMPLE_ID),]

seg4_bd1 <- seg4[seg4$CLIENT_SAMPLE_ID%in%seg4_cadets & seg4$CUSTOM_ATTRIBUTE_2 == "Blood Draw 1",]
seg4_bd2 <- seg4[seg4$CLIENT_SAMPLE_ID%in%seg4_cadets & seg4$CUSTOM_ATTRIBUTE_2 == "Blood Draw 2",]
seg4_bd3 <- seg4[seg4$CLIENT_SAMPLE_ID%in%seg4_cadets & seg4$CUSTOM_ATTRIBUTE_2 == "Blood Draw 3",]
seg4_bd4 <- seg4[seg4$CLIENT_SAMPLE_ID%in%seg4_cadets & seg4$CUSTOM_ATTRIBUTE_2 == "Blood Draw 4",]

seg4_bd1 <- seg4_bd1[!duplicated(seg4_bd1$CLIENT_SAMPLE_ID),]
seg4_bd2 <- seg4_bd2[!duplicated(seg4_bd2$CLIENT_SAMPLE_ID),]
seg4_bd3 <- seg4_bd3[!duplicated(seg4_bd3$CLIENT_SAMPLE_ID),]
seg4_bd4 <- seg4_bd4[!duplicated(seg4_bd4$CLIENT_SAMPLE_ID),]

seg5_bd1 <- seg5[seg5$CLIENT_SAMPLE_ID%in%seg5_cadets & seg5$CUSTOM_ATTRIBUTE_2 == "Blood Draw 1",]
seg5_bd2 <- seg5[seg5$CLIENT_SAMPLE_ID%in%seg5_cadets & seg5$CUSTOM_ATTRIBUTE_2 == "Blood Draw 2",]
seg5_bd3 <- seg5[seg5$CLIENT_SAMPLE_ID%in%seg5_cadets & seg5$CUSTOM_ATTRIBUTE_2 == "Blood Draw 3",]
seg5_bd4 <- seg5[seg5$CLIENT_SAMPLE_ID%in%seg5_cadets & seg5$CUSTOM_ATTRIBUTE_2 == "Blood Draw 4",]

seg5_bd1 <- seg5_bd1[!duplicated(seg5_bd1$CLIENT_SAMPLE_ID),]
seg5_bd2 <- seg5_bd2[!duplicated(seg5_bd2$CLIENT_SAMPLE_ID),]
seg5_bd3 <- seg5_bd3[!duplicated(seg5_bd3$CLIENT_SAMPLE_ID),]
seg5_bd4 <- seg5_bd4[!duplicated(seg5_bd4$CLIENT_SAMPLE_ID),]

#Check for cadet alignment
all.equal(seg2_cadets, seg2_bd1$CLIENT_SAMPLE_ID)
all.equal(seg2_cadets, seg2_bd2$CLIENT_SAMPLE_ID)
all.equal(seg2_cadets_2, seg2_bd3$CLIENT_SAMPLE_ID)
all.equal(seg2_cadets_2, seg2_bd4$CLIENT_SAMPLE_ID)

all.equal(seg4_cadets, seg4_bd1$CLIENT_SAMPLE_ID)
all.equal(seg4_cadets, seg4_bd2$CLIENT_SAMPLE_ID)
all.equal(seg4_cadets, seg4_bd3$CLIENT_SAMPLE_ID)
all.equal(seg4_cadets, seg4_bd4$CLIENT_SAMPLE_ID)

all.equal(seg5_cadets, seg5_bd1$CLIENT_SAMPLE_ID)
all.equal(seg5_cadets, seg5_bd2$CLIENT_SAMPLE_ID)
all.equal(seg5_cadets, seg5_bd3$CLIENT_SAMPLE_ID)
all.equal(seg5_cadets, seg5_bd4$CLIENT_SAMPLE_ID)

#Calculate log2-fold change.
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

seg2_cadets_prot_bd1 <- read.csv('log_prot_S2BD1.csv',header=T)[,1]
seg2_cadets_prot_bd2 <- read.csv('log_prot_S2BD2.csv',header=T)[,1]
seg2_cadets_prot_bd3 <- read.csv('log_prot_S2BD3.csv',header=T)[,1]

seg2_bd1_prot <- log_prot_sans_na1[seg2_cadets_prot_bd1%in%seg2_cadets,]
seg2_bd2_prot <- log_prot_sans_na2[seg2_cadets_prot_bd2%in%seg2_cadets,]
seg2_bd3_prot <- log_prot_sans_na3[seg2_cadets_prot_bd3%in%seg2_cadets,]

temp <- colnames(seg2_bd1_prot)[colnames(seg2_bd1_prot)%in%colnames(seg2_bd2_prot)]
shared_cols <- temp[temp%in%colnames(seg2_bd3_prot)]

seg2_bd1_prot_c <- seg2_bd1_prot[,shared_cols]
seg2_bd2_prot_c <- seg2_bd2_prot[,shared_cols]
seg2_bd3_prot_c <- seg2_bd3_prot[,shared_cols]

temp <- seg2_cadets_prot_bd1[seg2_cadets_prot_bd1%in%seg2_cadets]
temp2 <- seg2_cadets_prot_bd2[seg2_cadets_prot_bd2%in%seg2_cadets]
temp3 <- seg2_cadets_prot_bd3[seg2_cadets_prot_bd3%in%seg2_cadets]

seg2_bd1_prot_cc <- seg2_bd1_prot[temp%in%temp2,shared_cols]
seg2_bd2_prot_cc <- seg2_bd2_prot[temp2%in%temp2,shared_cols]
seg2_bd3_prot_cc <- seg2_bd3_prot[temp3%in%temp2,shared_cols]

seg2_bd1_np <- cbind(seg2_bd1_n,seg2_bd1_prot)
seg2_bd2_np <- cbind(seg2_bd2_n[seg2_cadets%in%seg2_cadets_prot_bd2,],seg2_bd2_prot)
seg2_bd3_np <- cbind(seg2_bd3_n,seg2_bd3_prot)

seg2_bd1_np <- seg2_bd1_np[seg2_cadets_prot_bd1[(seg2_cadets_prot_bd1%in%seg2_cadets)]%in%seg2_cadets_prot_bd2,]
seg2_bd3_np <- seg2_bd3_np[seg2_cadets_prot_bd3[(seg2_cadets_prot_bd3%in%seg2_cadets)]%in%seg2_cadets_prot_bd2,]

seg2_bd1_npc <- cbind(seg2_bd1_n,seg2_bd1_prot_c)
seg2_bd2_npc <- cbind(seg2_bd2_n[seg2_cadets%in%seg2_cadets_prot_bd2,],seg2_bd2_prot_c)
seg2_bd3_npc <- cbind(seg2_bd3_n,seg2_bd3_prot_c)

seg2_bd1_npc <- seg2_bd1_npc[seg2_cadets_prot_bd1[(seg2_cadets_prot_bd1%in%seg2_cadets)]%in%seg2_cadets_prot_bd2,]
seg2_bd3_npc <- seg2_bd3_npc[seg2_cadets_prot_bd3[(seg2_cadets_prot_bd3%in%seg2_cadets)]%in%seg2_cadets_prot_bd2,]

seg2_l21_npc <- log(exp(seg2_bd2_npc-seg2_bd1_npc),2)
seg2_l32_npc <- log(exp(seg2_bd3_npc-seg2_bd2_npc),2)

elim_no_var <- function(x){
  #For eliminating features without any variance
  y <- apply(x,2,var)
  y <- t(matrix(rep(as.logical(y != 0),nrow(x)),ncol=nrow(x)))
  x <- x*y
  x <- x[,colSums(x != 0) != 0]
  x
}

#Additional processing needed; adding a few target variables.
#First, set the blood draw 2 metabolites dataframe as the main dataset.
#Also create a separate data-frame tracking the desired phenotypes.
x1 = seg2_bd2_np
vo2 <- read_excel(paste(path,'DARPA-MBA-VO2max-InBody-Spring21WestPoint-08142021-AA.xlsx',sep=''))
vo2 <- vo2[vo2$`Cadet ID`%in%seg2_cadets,]
x2 = cbind(seg2_pwr[seg2_cadets%in%seg2_cadets_prot_bd2],
           seg2_hr[seg2_cadets%in%seg2_cadets_prot_bd2],
           as.numeric(vo2$`VO2 Max (mL/Kg/min)`[seg2_cadets%in%seg2_cadets_prot_bd2]),
           as.numeric(vo2$`Power @ VO2max (From Fitmate)- W`[seg2_cadets%in%seg2_cadets_prot_bd2]),
           as.numeric(vo2$`Percent Body Fat (%)`[seg2_cadets%in%seg2_cadets_prot_bd2]),
           as.numeric(vo2$`Skeletal Muscle Mass (lbs)`[seg2_cadets%in%seg2_cadets_prot_bd2])
)

#There should be a few NA's.
#Get rid of rows associated with those.
check <- apply(x2,1,sum)
leave <- -(1:length(check))[is.na(check)]
if (length(leave)!=0){
  x2 <- x2[leave,]
  x1 <- x1[leave,]
}
x1 <- elim_no_var(x1)
x1 <- scale(x1)

#Give proper column names for the taregts
x2 <- scale(x2)
colnames(x2) <- c("seg2pwr", "seg2hr", "vo2volume", 
                  "vo2power", "vo2bodyFat", "vo2muscleMass")

#Eliminate another row with a data integrity issue.
colsize <- ncol(x1)
set.seed(123)
x1m <- cbind(x1, matrix(rnorm(n=nrow(x1) * 50 ,sd=1), nrow = nrow(x1)))
x1m <- x1m[rownames(x1m)!='8',]
x2m <- x2[-7,]
x1m=scale(x1m)
x2m=scale(x2m)

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


sampleTree = hclust(dist(x1m), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# Plot a line to show the cut
abline(h = 62, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 62, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
x1mo = x1m[keepSamples, ]
x1mo = scale(x1mo)

#No major outliers...

#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Re-cluster samples
sampleTree2 = hclust(dist(x1mo), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")