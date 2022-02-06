#Data manipulation packs
library(readxl)
library(reshape2)
library(stringr)
library(boot)
library(openxlsx)

#Bioconductor packs
library(Biobase)
library(convert)

#Load main data
path <- "/Volumes/GoogleDrive-115111199924997198421/My Drive/Objectives/Code/"
setwd(path)
file <- 'Preprocessed_UntargetedMetabolomics-ReNormalizedData-DK.xlsx'
peak <- read_excel(file, sheet = 'Peak Area Data')
renorm <- read_excel(file, sheet = 'MBA-Renormalized')

#Divided by the median so that each metabolite has a median of one
#This also ensures that the metrics become unit-less, allowing an easier
#comparison across multiple studies
peak_med <- apply(as.matrix(apply(peak[,-1], 2, as.numeric)), 2, median, na.rm = T)
peak_norm <- sweep(as.matrix(apply(peak[,-1], 2, as.numeric)), 2, peak_med, '/')
rownames(peak_norm) <- peak$PARENT_SAMPLE_NAME

#Separate out the IDs
temp <- renorm[,1]

#Make everything else numeric
renorm <- as.matrix(apply(renorm[,-1],2,as.numeric))
rownames(renorm) <- temp$PARENT_SAMPLE_NAME

#Take out some erroneous datapoints based on previous observations
renorm <- renorm[-c(699,700),]

#Sanity check
exp(sum(log(colnames(peak_norm) == colnames(renorm))))
exp(sum(log(rownames(peak_norm) == rownames(renorm))))

#Reverse engineer the normalization factor based on data received from Metabolon
ratios <- peak_norm / renorm
med_ratios <- apply(ratios, 2, median, na.rm = T)

#Plot an instance of ratios
plot(ratios[,grep("482",colnames(ratios))], 
     ylab = "Peak Area/Renormalized Value",
     xlab = "Sample Index",
     main = "Lactate (HMDB0000190)
     Peak Area to Renormalized Value Ratio",
     cex.lab = 1.4,
     cex.main = 1.4)
abline(h=med_ratio_distilled[grep("482",colnames(ratios))], lty =2, col = "red")
legend("topright", lty = 2, col = "red", legend = "Median")

#Import athletes data that will be merged into the subjects dataset
athletes <- read_excel("overlapping_metabolites v0.2dk.xlsx", sheet = "athlete_191_chem_id")
temp <- athletes$Subject_ID
athletes <- as.matrix(apply(athletes[,-1],2,as.numeric))
rownames(athletes) <- temp

#Create the list of matching metabolite indices
idx <- match(colnames(athletes),names(med_ratios))
idx <- idx[!is.na(idx)]

#Again normalize the athletes data so that the median is just one for each feat
athletes_med <- apply(athletes, 2, median,na.rm=T)
athletes_norm <- sweep(athletes, 2, athletes_med, '/')

#Take into account the median normalizing factor from the Metabolon data
#And "batch-normalize" the athletes data accordingly
med_ratio_distilled <- med_ratios[idx]
athletes_batch_norm <- sweep(athletes_norm,2,med_ratio_distilled,'/')

#Then adjust by the normalization factor
ath_bn_med <- apply(athletes_batch_norm,2,median,na.rm=T)
athletes_imputed <- sweep(athletes_batch_norm,2, ath_bn_med,'/')

#Then impute missing values
min_ath_imputed <-  apply(athletes_imputed,2,min,na.rm=T)
min_ath_imputed <- replace(min_ath_imputed,min_ath_imputed==Inf,NA)
idx_na <- is.na(athletes_imputed)
athletes_imputed[idx_na] <- min_ath_imputed[col(athletes_imputed)][idx_na]

#Finally, take logs
log_athletes <- log(athletes_imputed)

#Save the data
workbook <- loadWorkbook("overlapping_metabolites v0.2dk.xlsx", xlsxFile = NULL, isUnzipped = FALSE)
addWorksheet(workbook, "athlete_191_log_transf_data")
writeData(workbook, sheet="athlete_191_log_transf_data", x = log_athletes)
saveWorkbook(workbook, file, overwrite = T)
write.csv(log_athletes,'log_norm_dat_ath.csv')
