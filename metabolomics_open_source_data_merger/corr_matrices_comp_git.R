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

#Additional packs
library(psych)
library(BSDA)
library(DescTools)

#Import two correlation matrices
setwd('/Users/Derek/Downloads/')
ma_pre <- read_excel('Correlation_Matrices_Updated.xlsx', sheet = "ma_PreExercise", col_names=F)
ma_post <- read_excel('Correlation_Matrices_Updated.xlsx', sheet = "ma_PostExercise", col_names = F)
ma_fc <- read_excel('Correlation_Matrices_Updated.xlsx', sheet = "ma_FC", col_names = F)
c1pre <- read_excel('Correlation_Matrices_Updated.xlsx', sheet = "Cohort1_PreExercise", col_names = F)
c1post <- read_excel('Correlation_Matrices_Updated.xlsx', sheet = "Cohort1_PostExercise", col_names = F)
c1fc <- read_excel('Correlation_Matrices_Updated.xlsx', sheet = "Cohort1_FC", col_names = F)

#Dedupe
ma_pre <- ma_pre[!duplicated(ma_pre$...1),colnames(ma_pre)[!duplicated(as.character(ma_pre[1,]))]]
ma_post <- ma_post[!duplicated(ma_post$...1),colnames(ma_post)[!duplicated(as.character(ma_post[1,]))]]
ma_fc <- ma_fc[!duplicated(ma_fc$...1),colnames(ma_fc)[!duplicated(as.character(ma_fc[1,]))]]
c1pre <- c1pre[!duplicated(c1pre$...1),colnames(c1pre)[!duplicated(as.character(c1pre[1,]))]]
c1post <- c1post[!duplicated(c1post$...1),colnames(c1post)[!duplicated(as.character(c1post[1,]))]]
c1fc <- c1fc[!duplicated(c1fc$...1),colnames(c1fc)[!duplicated(as.character(c1fc[1,]))]]

#Eliminate first col and row and set row and col names
ma_pre <- ma_pre[,-1]
ma_post <- ma_post[,-1]
ma_fc <- ma_fc[,-1]
c1pre <- c1pre[,-1]
c1post <- c1post[,-1]
c1fc <- c1fc[,-1]

colnames(ma_pre) <- as.character(ma_pre[1,])
colnames(ma_post) <- as.character(ma_post[1,])
colnames(ma_fc) <- as.character(ma_fc[1,])
colnames(c1pre) <- as.character(c1pre[1,])
colnames(c1post) <- as.character(c1post[1,])
colnames(c1fc) <- as.character(c1fc[1,])

ma_pre <- ma_pre[-1,]
ma_post <- ma_post[-1,]
ma_fc <- ma_fc[-1,]
c1pre <- c1pre[-1,]
c1post <- c1post[-1,]
c1fc <- c1fc[-1,]

#Turn the dataframes into matrices for ease of computation
ma_pre <- apply(ma_pre,2,as.numeric)
ma_post <- apply(ma_post,2,as.numeric)
ma_fc <- apply(ma_fc,2,as.numeric)
c1pre <- apply(c1pre,2,as.numeric)
c1post <- apply(c1post,2,as.numeric)
c1fc <- apply(c1fc,2,as.numeric)

#Set row names
row.names(ma_pre) <- colnames(ma_pre)
row.names(ma_post) <- colnames(ma_post)
row.names(ma_fc) <- colnames(ma_fc)
row.names(c1pre) <- colnames(c1pre)
row.names(c1post) <- colnames(c1post)
row.names(c1fc) <- colnames(c1fc)

#Check whether correlation values are greater a threshold
threshold = 0.1
cutoff_mask_pre <- (abs(ma_pre) > threshold) & (abs(c1pre) > threshold)
cutoff_mask_post <- (abs(ma_post) > threshold) & (abs(c1post) > threshold)
cutoff_mask_fc <- (abs(ma_fc) > threshold) & (abs(c1fc) > threshold)

#Check whether correlation signs align
sgn_mask_pre <- sign(ma_pre) * sign(c1pre) > 0
sgn_mask_post <- sign(ma_post) * sign(c1post) > 0
sgn_mask_fc <- sign(ma_fc) * sign(c1fc) > 0

#Put the masks together
mask_pre <- cutoff_mask_pre & sgn_mask_pre & ma_pre != 1
mask_post <- cutoff_mask_post & sgn_mask_post & ma_post != 1
mask_fc <- cutoff_mask_fc & sgn_mask_fc & ma_fc != 1

#Set all NAs due to non-numerical instances to F
mask_pre[is.na(mask_pre)] <- F
mask_post[is.na(mask_post)] <- F
mask_fc[is.na(mask_fc)] <- F

#Make copies of the original matrices for safety
ma_pre_c <- ma_pre
ma_post_c <- ma_post
ma_fc_c <- ma_fc
c1pre_c <- c1pre
c1post_c <- c1post
c1fc_c <- c1fc

#Set the non-masked values to NA so that we can perform numerical operations
#only on the masked values
ma_pre_c[!mask_pre] <- NA
c1pre_c[!mask_pre] <- NA
ma_post_c[!mask_post] <- NA
c1post_c[!mask_post] <- NA
ma_fc_c[!mask_fc] <- NA
c1fc_c[!mask_fc] <- NA

#Select those whose correlation absolute differences are less than 10%
mask_pre_fin <- abs((ma_pre_c-c1pre_c)/c1pre_c) < 0.1
mask_pre_fin[is.na(mask_pre_fin)] <- F

mask_post_fin <- abs((ma_post_c-c1post_c)/c1post_c) < 0.1
mask_post_fin[is.na(mask_post_fin)] <- F

mask_fc_fin <- abs((ma_fc_c-c1fc_c)/c1fc_c) < 0.1
mask_fc_fin[is.na(mask_fc_fin)] <- F

#Select the row and column combination that meet all criteria
pre_out <- outer(colnames(ma_pre_c),colnames(ma_pre_c),paste,sep=", ")[mask_pre_fin]
post_out <- outer(colnames(ma_post_c),colnames(ma_post_c),paste,sep=", ")[mask_post_fin]
fc_out <- outer(colnames(ma_fc_c),colnames(ma_fc_c),paste,sep=", ")[mask_fc_fin]

#Save
write.csv(pre_out,"pre-stress_pairs_v0.1.csv")
write.csv(post_out,"post-stress_pairs_v0.1.csv")
write.csv(fc_out,"fc_pairs_v0.1.csv")

#Look for intersections of metabolites that appear across modalities
all <- Reduce(intersect,list(pre_out, post_out, fc_out))
pre_post <- Reduce(intersect,list(pre_out, post_out))
post_fc <- Reduce(intersect,list(post_out, fc_out))
pre_fc <- Reduce(intersect,list(pre_out, fc_out))

#Save the set differences
write.csv(Reduce(setdiff,list(pre_out,post_out,fc_out)),"unique.csv")
write.csv(Reduce(setdiff,list(post_out,pre_out,fc_out)),"unique.csv")
write.csv(Reduce(setdiff,list(fc_out,post_out,pre_out)), "unique.csv")
