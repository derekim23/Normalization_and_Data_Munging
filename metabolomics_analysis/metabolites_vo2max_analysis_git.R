#Data manipulation packs
library(readxl)
library(reshape2)
library(robustbase)
library(stringr)
library(boot)
library(splitstackshape)
library(grid)

#Bioconductor pack
library(Biobase)
library(convert)

#Additional packs
library(psych)
library(BSDA)
library(DescTools)

#Draw in VO2Max, hemoglobin, metabolites, etc. data
vo2path <- "/Volumes/GoogleDrive-115111199924997198421/My Drive/Objectives/Code/"
vo2_file <- "VO2max-InBody-Spring21-08142021.xlsx"
vo2 <- read_excel(paste(vo2path,vo2_file,sep=""),sheet="Basic InBody VO2 max")

hgb_path <- "/Users/Derek/Downloads/"
hgb_bd1 <- read_excel(paste(hgb_path, "iSTAT_Munged.xlsx", sep=""), sheet = 'S2BD1')
hgb_bd2 <- read_excel(paste(hgb_path, "iSTAT_Munged.xlsx", sep=""), sheet = "S2BD2")
hgb_bd1 <- hgb_bd1[,c('subject ID', 'Hgb (EC8+(4)) Result for Blood Draw 1')]
hgb_bd2 <- hgb_bd2[,c('subject ID', 'Hgb (EC8+(4)) Result for Blood Draw 2')]
colnames(hgb_bd2)[2] <- "hgb"
colnames(hgb_bd1)[2] <- "hgb"
cog <- read_excel(paste(hgb_path, "AMS-Spring-01052022-dk.xlsx", sep = ""), sheet = "Sheet2")

metab_file <- 'Preprocessed_UntargetedMetabolomics-ReNormalizedData-Spring21-09272021.xlsx'
metab <- read_excel(paste(vo2path, metab_file, sep = ""), 
                    sheet = "Log-Transformed and Merged Data")
short_list <- read_excel(paste(vo2path, metab_file, sep = ""), 
                         sheet = "short_list_metab")

#Focus just on Segment 2
metab_bd1 <- metab[metab$CUSTOM_ATTRIBUTE_1 == "Segment 2" & 
                     metab$CUSTOM_ATTRIBUTE_2 == "Blood Draw 1",]
metab_bd2 <- metab[metab$CUSTOM_ATTRIBUTE_1 == "Segment 2" & 
                     metab$CUSTOM_ATTRIBUTE_2 == "Blood Draw 2",]

#Dedupe by subject ID
metab_bd1 <- metab_bd1[!duplicated(metab_bd1$CLIENT_SAMPLE_ID),]
metab_bd2 <- metab_bd2[!duplicated(metab_bd2$CLIENT_SAMPLE_ID),]
metab_bd1 <- metab_bd1[order(metab_bd1$CLIENT_SAMPLE_ID),-1]
metab_bd2 <- metab_bd2[order(metab_bd2$CLIENT_SAMPLE_ID),-1]

#Calculate tyrosine fold change
df <- exp(metab_bd2[metab_bd2$CLIENT_SAMPLE_ID%in%metab_bd1$CLIENT_SAMPLE_ID,
                    c("815")] - 
  metab_bd1[metab_bd1$CLIENT_SAMPLE_ID%in%metab_bd2$CLIENT_SAMPLE_ID,
                                          c("815")])
df[,"subject ID"] <- metab_bd1$CLIENT_SAMPLE_ID[metab_bd1$CLIENT_SAMPLE_ID%in%
                                                    metab_bd2$CLIENT_SAMPLE_ID]
colnames(df)[1] <- "tyr_fc"

df <- merge(df, metab_bd1[,c("CLIENT_SAMPLE_ID","815")], 
                by.x = "subject ID", by.y = "CLIENT_SAMPLE_ID")
df <- merge(df, metab_bd2[,c("CLIENT_SAMPLE_ID","815")], 
                by.x = "subject ID", by.y = "CLIENT_SAMPLE_ID")
df <- merge(df, hgb_bd1, by = "subject ID")
df <- merge(df, hgb_bd2, by = "subject ID")

#Also merge pyruvate, lactate, tryptophan, creatine, n-acetyl glutamine, 
#phenylalanine, leucine, isoleucine, valine, picolinate, kynurenate, kynurenine,
#n-acetyl-aspartyl-glutamate (NAAG)
df <- merge(df, metab_bd1[,c('CLIENT_SAMPLE_ID', '823', '482', '565', 
                                     '1221', '100001253', '460', '397', '376',
                                     '566', '1022', '98', '100000265', '100001612',
                             as.character(short_list$CHEM_ID[1:264]))], 
                by.x = "subject ID",
                by.y = "CLIENT_SAMPLE_ID")
df <- merge(df, metab_bd2[,c('CLIENT_SAMPLE_ID', '823', '482', '565', 
                             '1221', '100001253', '460', '397', '376',
                             '566', '1022', '98', '100000265', '100001612', 
                             as.character(short_list$CHEM_ID[1:264]))], 
            by.x = "subject ID",
            by.y = "CLIENT_SAMPLE_ID")

#Calculate desired ratios
df[,c('trp/phen1', 'trp/leuc1', 'trp/isol1','trp/vali1')] <- 
  exp(df[,'565.x']-df[,c('460.x', '397.x', '376.x','566.x')]) 
df[,'trp/(t + p + l + i + v)1'] <- 
  exp(df[,'565.x'])/apply(exp(df[,c('460.x', '397.x', '376.x','566.x')]),1,sum) 
df[,'trp/tyr1'] <- exp(df[,'565.x'])/exp(df[,'815.x']) 

df[,c('trp/phen2', 'trp/leuc2', 'trp/isol2','trp/vali2')] <- 
  exp(df[,'565.y']-df[,c('460.y', '397.y', '376.y','566.y')]) 
df[,'trp/(t + p + l + i + v)2'] <- 
  exp(df[,'565.y'])/apply(exp(df[,c('460.y', '397.y', '376.y','566.y')]),1,sum)
df[,'trp/tyr2'] <- exp(df[,'565.y'])/exp(df[,'815.y']) 

df[,'glut/crea1'] <- exp(df[,'100001253.x'])/exp(df[,'1221.x']) 

df[,'glut/crea2'] <- exp(df[,'100001253.y'])/exp(df[,'1221.y'])

df[,'trp_fc'] <- exp(df[,'565.y'] - df[,'565.x'])

df[,'kyn/trp1'] <- exp(df[,'100000265.x'])/exp(df[,'565.x'])

df[,'kyn/trp2'] <- exp(df[,'100000265.y'])/exp(df[,'565.y'])

df[,'kyna/trp1'] <- exp(df[,'98.x'])/exp(df[,'565.x'])

df[,'kyna/trp2'] <- exp(df[,'98.y'])/exp(df[,'565.y'])

df[,'kyn/trp_fc'] <- df[,'kyn/trp2']/df[,'kyn/trp1']

#Take a short_list of just chem_ids
chem_list <- as.character(short_list$CHEM_ID)[1:264]

#Also, create fold change values for the short_list
df[,paste(chem_list,"_fc",sep="")] <- exp(df[,paste(chem_list,".y",sep="")]-df[,paste(chem_list,".x",sep="")])

#Merge VO2Max volume and power
df <- merge(df, vo2[,c("subject ID", "VO2 Max (mL/Kg/min)", "VO2 Rank (Bin)",
                       "Endurance Level on the Bike- VO2 max Test",
                       "Power @ VO2 max (From Pedals/Garmin)- IMPUTED", 
                       "Skeletal Muscle Mass (lbs)", "Dynamic")], 
                by = "subject ID")
df$`Power @ VO2 max (From Pedals/Garmin)- IMPUTED` <- 
  as.numeric(df$`Power @ VO2 max (From Pedals/Garmin)- IMPUTED`)

#Set infinite values to NA
df$`glut/crea1`[is.infinite(df$`glut/crea1`)] <- NA

#Merge cognition data into main dataframe
df <- merge(df, cog, by.x = "subject ID", by.y = "ID")
df$`subject Record Brief (Average msc)` <- as.numeric(df$`subject Record Brief (Average msc)`)
df$cps_Cum <- as.numeric(df$cps_Cum)
df$GPA_Cum <- as.numeric(df$GPA_Cum)
df$CPS <- as.numeric(df$CPS)
df$`Whole Candidate Score` <- as.numeric(df$`Whole Candidate Score`)
df$`strength score- 2021 (sc Cum)` <- as.numeric(df$`strength score- 2021 (sc Cum)`)
df$Dynamic <- as.numeric(df$Dynamic)

#Create two categories for each of power and VO2Max, which will serve as
#proxies for anaerobic/aerobic fitness
df$`VO2 Rank (Bin)`[df$`VO2 Rank (Bin)`=="Very Poor"]=1
df$`VO2 Rank (Bin)`[df$`VO2 Rank (Bin)`=="Poor"]=2
df$`VO2 Rank (Bin)`[df$`VO2 Rank (Bin)`=="Fair"]=3
df$`VO2 Rank (Bin)`[df$`VO2 Rank (Bin)`=="Good"]=4
df$`VO2 Rank (Bin)`[df$`VO2 Rank (Bin)`=="Excellent"]=5
df$`VO2 Rank (Bin)`[df$`VO2 Rank (Bin)`=="Superior"]=1
df$`VO2 Rank (Bin)` <- as.numeric(df$`VO2 Rank (Bin)`)

df[,'aerobic_cat'] <- df$`VO2 Rank (Bin)`>3
df[,'aerobic_cat'] <- df[,'aerobic_cat']*1

df[,'anaerobic_cat'] <- df$`Power @ VO2 max (From Pedals/Garmin)- IMPUTED`>
  median(df$`Power @ VO2 max (From Pedals/Garmin)- IMPUTED`, na.rm = T)
df[,'anaerobic_cat'] <- df[,'anaerobic_cat']*1

df[,'cps_cat'] <- df$CPS>median(df$CPS,na.rm=T)
df[,'cps_cat'] <- df[,'cps_cat']*1

#Also, bundle novice and intermediate endurance athletes as there are too few
#of novice athletes
df$`Endurance Level on the Bike- VO2 max Test`[
  df$`Endurance Level on the Bike- VO2 max Test` == "Intermediate"] <- 
  "Novice/Intermediate"
df$`Endurance Level on the Bike- VO2 max Test`[
  df$`Endurance Level on the Bike- VO2 max Test` == "Novice"] <- 
  "Novice/Intermediate"

#Calculate p_values for the low vs. high endurance pairs for the
#short_list metabolites
p_vals_x <- c()
cohen_x <- c()
errors_x <- c()
p_vals_y <- c()
cohen_y <- c()
errors_y <- c()
p_vals_fc <- c()
cohen_fc <- c()
errors_fc <- c()
for (i in 1:264){
  err_x <- tryCatch({
    p_vals_x <- c(p_vals_x, wilcox.test(df[df$Dynamic==1,paste(chem_list[i],".x",sep="")],
                                    df[df$Dynamic==2,paste(chem_list[i],".x",sep="")])$p.value)
    },
    error = function(w) w)
  if(inherits(err_x,"warning"))
     errors_x <- c(errors_x,i)
  err_y <- tryCatch({
    p_vals_y <- c(p_vals_y, wilcox.test(df[df$Dynamic==1,paste(chem_list[i],".y",sep="")],
                                      df[df$Dynamic==2,paste(chem_list[i],".y",sep="")])$p.value)
  },
  error = function(w) w)
  if(inherits(err_y,"warning"))
    errors_y <- c(errors_y,i)
  err_fc <- tryCatch({
    p_vals_fc <- c(p_vals_fc, wilcox.test(df[df$Dynamic==1,paste(chem_list[i],"_fc",sep="")],
                                       df[df$Dynamic==2,paste(chem_list[i],"_fc",sep="")])$p.value)
  },
  error = function(w) w)
  if(inherits(err_fc,"warning"))
    errors_fc <- c(errors_fc,i)
}

#Look at the min adjusted p-values:
min(p.adjust(p_vals_x, "fdr"), na.rm = T)
min(p.adjust(p_vals_y, "fdr"), na.rm = T)
min(p.adjust(p_vals_fc, "fdr"), na.rm = T)

#Look at which are significant pre-adjustment
raw_sig_x <- rm.na(chem_list[p_vals_x< 0.05])
raw_sig_y <- rm.na(chem_list[p_vals_y< 0.05])
raw_sig_fc <- rm.na(chem_list[p_vals_fc< 0.05])

#Create a matrix of hypothesis test results instead
stats <- matrix(NA, nrow = 264, ncol = 18)
colnames(stats) <- c("HMDB", "Chem_ID", "Name", "Pre_p_Val", "Pre_t_Val", 
                     "Pre_Cohen", "Pre_FWER", "Pre_FDR", "Post_p_Val", 
                     "Post_t_Val", "Post_Cohen", "Post_FWER", "Post_FDR",
                     "FC_p_Val", "FC_t_Val", "FC_Cohen", "FC_FWER", 
                     "FC_FDR")
                     
for (i in 1:264){
  stats[i,1] <- short_list$HMDB[grep(chem_list[i],short_list$CHEM_ID)[1]]
  stats[i,2] <- chem_list[i]
  stats[i,3] <- short_list$NAME[grep(chem_list[i],short_list$CHEM_ID)[1]]
  err_x <- tryCatch({
    stats[i,4] <- t.test(df[df$Dynamic==1,paste(chem_list[i],".x",sep="")],
                         df[df$Dynamic==2,paste(chem_list[i],".x",sep="")])$p.value
    stats[i,5] <- signif(t.test(df[df$Dynamic==1,paste(chem_list[i],".x",sep="")],
                         df[df$Dynamic==2,paste(chem_list[i],".x",sep="")])$statistic,2)
    stats[i,6] <- signif(cohen.d(df[,paste(chem_list[i],".x",sep="")], df$Dynamic)$cohen.d[2],2)
  },
  error = function(w) w)
  if(inherits(err_x,"warning"))
    errors_x <- c(errors_x,i)
  err_y <- tryCatch({
    stats[i,9] <- t.test(df[df$Dynamic==1,paste(chem_list[i],".y",sep="")],
                         df[df$Dynamic==2,paste(chem_list[i],".y",sep="")])$p.value
    stats[i,10] <- signif(t.test(df[df$Dynamic==1,paste(chem_list[i],".y",sep="")],
                         df[df$Dynamic==2,paste(chem_list[i],".y",sep="")])$statistic,2)
    stats[i,11] <- signif(cohen.d(df[,paste(chem_list[i],".y",sep="")], df$Dynamic)$cohen.d[2],2)
  },
  error = function(w) w)
  if(inherits(err_y,"warning"))
    errors_y <- c(errors_y,i)
  err_fc <- tryCatch({
    stats[i,14] <- t.test(df[df$Dynamic==1,paste(chem_list[i],"_fc",sep="")],
                         df[df$Dynamic==2,paste(chem_list[i],"_fc",sep="")])$p.value
    stats[i,15] <- signif(t.test(df[df$Dynamic==1,paste(chem_list[i],"_fc",sep="")],
                         df[df$Dynamic==2,paste(chem_list[i],"_fc",sep="")])$statistic,2)
    stats[i,16] <- signif(cohen.d(df[,paste(chem_list[i],"_fc",sep="")], df$Dynamic)$cohen.d[2],2)
  },
  error = function(w) w)
  if(inherits(err_fc,"warning"))
    errors_fc <- c(errors_fc,i)
}

#Find whether we have significance after adjustment
stats[,7] <- p.adjust(as.numeric(stats[,4]),"hochberg") < 0.05
stats[,8] <- p.adjust(as.numeric(stats[,4]),"fdr") < 0.05
stats[,12] <- p.adjust(as.numeric(stats[,9]),"hochberg") < 0.05
stats[,13] <- p.adjust(as.numeric(stats[,9]),"fdr") < 0.05
stats[,17] <- p.adjust(as.numeric(stats[,14]),"hochberg") < 0.05
stats[,18] <- p.adjust(as.numeric(stats[,14]),"fdr") < 0.05

#Now find significant digits of p-values
stats[,4] <- signif(as.numeric(stats[,4]),2)
stats[,9] <- signif(as.numeric(stats[,9]),2)
stats[,14] <- signif(as.numeric(stats[,14]),2)

write.csv(stats,"/Users/Derek/Downloads/p_vals.csv")

#Also, make the same comparison between pre- and post-stress 
#for Dynamic == 1 sub-cohort
stats3 <- matrix(NA, nrow = 264, ncol = 8)

for (i in 1:264){
  stats3[i,1] <- short_list$HMDB[grep(chem_list[i],short_list$CHEM_ID)[1]]
  stats3[i,2] <- chem_list[i]
  stats3[i,3] <- short_list$NAME[grep(chem_list[i],short_list$CHEM_ID)[1]]
  err_x <- tryCatch({
    stats3[i,4] <- t.test(df[df$Dynamic==1,paste(chem_list[i],".x",sep="")],
                         df[df$Dynamic==1,paste(chem_list[i],".y",sep="")], paired = T)$p.value
    stats3[i,5] <- signif(t.test(df[df$Dynamic==1,paste(chem_list[i],".x",sep="")],
                                df[df$Dynamic==1,paste(chem_list[i],".y",sep="")], paired = T)$statistic,2)
    stats3[i,6] <- signif(cohen.d(df[df$Dynamic==1,paste(chem_list[i],".y",sep="")], df$Dynamic)$cohen.d[2],2)
  },
  error = function(w) w)
  if(inherits(err_x,"warning"))
    errors_x <- c(errors_x,i)
}
colnames(stats3) <- c("HMDB", "Chem_ID", "Name", "p_Val", "t_Val", 
                     "Cohen", "FWER", "FDR")

stats3[,7] <- p.adjust(as.numeric(stats3[,4]),"hochberg") < 0.05
stats3[,8] <- p.adjust(as.numeric(stats3[,4]),"fdr") < 0.05
stats3[,4] <- signif(as.numeric(stats3[,4]),2)

write.csv(stats3,"/Users/Derek/Downloads/p_vals_low_endur.csv")

#Also, make the same comparison between pre- and post-stress 
#for Dynamic == 2 sub-cohort
stats4 <- matrix(NA, nrow = 264, ncol = 8)

for (i in 1:264){
  stats4[i,1] <- short_list$HMDB[grep(chem_list[i],short_list$CHEM_ID)[1]]
  stats4[i,2] <- chem_list[i]
  stats4[i,3] <- short_list$NAME[grep(chem_list[i],short_list$CHEM_ID)[1]]
  err_x <- tryCatch({
    stats4[i,4] <- t.test(df[df$Dynamic==2,paste(chem_list[i],".x",sep="")],
                          df[df$Dynamic==2,paste(chem_list[i],".y",sep="")], paired = T)$p.value
    stats4[i,5] <- signif(t.test(df[df$Dynamic==2,paste(chem_list[i],".x",sep="")],
                                 df[df$Dynamic==2,paste(chem_list[i],".y",sep="")], paired = T)$statistic,2)
    stats4[i,6] <- signif(cohen.d(df[df$Dynamic==2,paste(chem_list[i],".y",sep="")], df$Dynamic)$cohen.d[2],2)
  },
  error = function(w) w)
  if(inherits(err_x,"warning"))
    errors_x <- c(errors_x,i)
}
colnames(stats4) <- c("HMDB", "Chem_ID", "Name", "p_Val", "t_Val", 
                      "Cohen", "FWER", "FDR")

stats4[,7] <- p.adjust(as.numeric(stats4[,4]),"hochberg") < 0.05
stats4[,8] <- p.adjust(as.numeric(stats4[,4]),"fdr") < 0.05
stats4[,4] <- signif(as.numeric(stats4[,4]),2)

write.csv(stats4,"/Users/Derek/Downloads/p_vals_high_endur.csv")

#Create boxplots and save as jpegs
for (i in 1:length(raw_sig_x)){
  print(i)
  jpeg(file = paste('/Users/Derek/Downloads/pre_boxplots/',
                    paste(str_to_title(short_list$NAME[
                      grep(raw_sig_x[i],short_list$CHEM_ID)[i]]), " (", 
                          short_list$HMDB[grep(raw_sig_x[1],short_list$CHEM_ID)[1]], 
                      ")", sep = ""),
                          '.jpeg', sep = ''))
  boxplot(df[,paste(raw_sig_x[i],".x",sep="")]~df$Dynamic, 
          col = c('lightpink', 'lightblue'), horizontal = FALSE,
          names = c("Low", "High"),
          xlab = "Endurance",
          ylab = "Z-Scores (Pre-Stress)",
          main = paste(str_to_title(short_list$NAME[
            grep(raw_sig_x[i],short_list$CHEM_ID)[i]]), " (", 
            short_list$HMDB[grep(raw_sig_x[1],short_list$CHEM_ID)[1]], 
            ")", sep = ""),
          cex.lab =1.1)
  grid.text(paste("p-value = ", 
                  round(wilcox.test(df[df$Dynamic==1,paste(raw_sig_x[i],".x",sep="")],
                                    df[df$Dynamic==2,paste(raw_sig_x[i],".x",sep="")])$p.value, 4),
                  sep = ""), x=unit(2, "mm"), y=unit(1, "npc") - unit(2, "mm"),
            just=c("left", "top"), gp = gpar(fontsize = 14))
  dev.off()
}

#See if spermidine abundance is significantly different between two sub-cohorts w/ boxpots
boxplot(df[,paste(chem_list[1],".x",sep="")]~df$Dynamic, 
        col = c('lightpink', 'lightblue'), horizontal = FALSE,
        names = c("Low", "High"),
        xlab = "Endurance",
        ylab = "Z-Scores (Pre-Stress)",
        main = paste(str_to_title(short_list$NAME[1]), " (", 
                     short_list$HMDB[1], ")", sep = ""),
        cex.lab =1.1)

grid.text(paste("p-value = ", 
                round(wilcox.test(df[df$Dynamic==1,paste(chem_list[1],".x",sep="")],
                            df[df$Dynamic==2,paste(chem_list[1],".x",sep="")])$p.value, 4),
                sep = ""), x=unit(2, "mm"), y=unit(1, "npc") - unit(2, "mm"),
          just=c("left", "top"), gp = gpar(fontsize = 14))


#Try fitting linear models for molecules vs. glycolitic, anaerobic, and aerobic targets
lm_hgb_tyr_fc <- lm(hgb.x~`tyr_fc`,data = df)
slm_hgb_tyr_fc <- summary(lm_hgb_tyr)

lm_lac_tyr <- lm(`482.y`~`815.x`,data = df)
slm_lac_tyr <- summary(lm_lac_tyr)

lm_pyr_tyr <- lm(`823.y`~`815.y`,data = df)
slm_pyr_tyr <- summary(lm_pyr_tyr)

lm_pyr_trp_fc <- lm(`823.y`~trp_fc,data=df)
slm_pyr_trp_fc <- summary(lm_pyr_trp_fc)

lm_lac_trp_fc <- lm(`482.y`~trp_fc,data=df)
slm_lac_trp_fc <- summary(lm_lac_trp_fc)

lm_pyr_cre <- lm(`823.y`~`1221.y`,data=df)
slm_pyr_cre <- summary(lm_pyr_cre)

lm_lac_cre <- lm(`482.y`~`1221.y`,data=df)
slm_lac_cre <- summary(lm_lac_cre)

lm_vo2_cre <- lm(`VO2 Max (mL/Kg/min)`~`1221.y`,data=df)
slm_vo2_cre <- summary(lm_vo2_cre)

lm_pwr_cre <- lm(`Power @ VO2 max (From Pedals/Garmin)- IMPUTED`~`1221.y`,data=df)
slm_pwr_cre <- summary(lm_pwr_cre)

lm_vo2_glu_cre <-  lm(`VO2 Max (mL/Kg/min)`~`glut/crea1`,data=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,])
slm_vo2_glu_cre <- summary(lm_vo2_glu_cre)

lm_pyr_glu_cre <-  lm(`823.y`~`glut/crea1`,data=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,])
slm_pyr_glu_cre <- summary(lm_pyr_glu_cre)

lm_lac_glu_cre <-  lm(`482.y`~`glut/crea1`,data=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,])
slm_lac_glu_cre <- summary(lm_lac_glu_cre)

lm_pwr_glu_cre <-  lm(`Power @ VO2 max (From Pedals/Garmin)- IMPUTED`~`glut/crea1`,
                      data=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,])
slm_pwr_glu_cre <- summary(lm_pwr_glu_cre)

lm_mm_glu_cre <-  lm(`Skeletal Muscle Mass (lbs)`~`glut/crea1`,
                      data=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,])
slm_mm_glu_cre <- summary(lm_mm_glu_cre)

lm_mm_glu_cre <-  lm(`Skeletal Muscle Mass (lbs)`~`glut/crea1`,
                     data=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,])
slm_mm_glu_cre <- summary(lm_mm_glu_cre)

lm_hgb_glu_cre <-  lm(hgb.y~`glut/crea1`,
                     data=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,])
slm_hgb_glu_cre <- summary(lm_hgb_glu_cre)

lm_lac_trp_tyr <- lm(`482.y`~`trp/tyr1`,data=df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,])
slm_lac_trp_tyr <- summary(lm_lac_trp_tyr)

lm_pyr_trp_tyr <- lm(`823.y`~`trp/tyr1`,data=df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,])
slm_pyr_trp_tyr <- summary(lm_pyr_trp_tyr)

lm_hgb_trp_tyr <- lm(`hgb.y`~`trp/tyr1`,data=df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,])
slm_hgb_trp_tyr <- summary(lm_hgb_trp_tyr)

lm_hgb_trp_tpliv <- lm(`hgb.y`~`trp/(t + p + l + i + v)1`,data=df[df$`trp/(t + p + l + i + v)1`<25,])
slm_hgb_trp_tpliv <- summary(lm_hgb_trp_tpliv)

lm_pyr_trp_tpliv <- lm(`823.y`~`trp/(t + p + l + i + v)1`,data=df[df$`trp/(t + p + l + i + v)1`<25,])
slm_pyr_trp_tpliv <- summary(lm_pyr_trp_tpliv)

lm_lac_trp_tpliv <- lm(`482.y`~`trp/(t + p + l + i + v)1`,data=df[df$`trp/(t + p + l + i + v)1`<25,])
slm_lac_trp_tpliv <- summary(lm_lac_trp_tpliv)

lm_pic_trp_tpliv <- lm(`1022.y`~`trp/(t + p + l + i + v)1`,data=df[df$`trp/(t + p + l + i + v)1`<25,])
slm_pic_trp_tpliv <- summary(lm_pic_trp_tpliv)

lm_kyn_trp_tpliv <- lm(`98.y`~`trp/(t + p + l + i + v)1`,data=df[df$`trp/(t + p + l + i + v)1`<25,])
slm_kyn_trp_tpliv <- summary(lm_kyn_trp_tpliv)

lm_pic_trp_tyr <- lm(`1022.y`~`trp/tyr1`,data=df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,])
slm_pic_trp_tyr <- summary(lm_pic_trp_tyr)

lm_kyn_trp_tyr <- lm(`98.y`~`trp/tyr1`,data=df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,])
slm_kyn_trp_tyr <- summary(lm_kyn_trp_tyr)

lm_kyn_trp_tyr <- lm(`98.y`~`trp/tyr1`,data=df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,])
slm_kyn_trp_tyr <- summary(lm_kyn_trp_tyr)

lm_msc_kyn_trp <- lm(`subject Record Brief (Average msc)`~`kyn/trp1`,
                      data=df[df$`kyn/trp1`<60,])
slm_msc_kyn_trp <- summary(lm_msc_kyn_trp)

lm_cps_kyna_trp <- lm(CPS~`kyna/trp1`,
                      data=df[df$`kyna/trp1`<60,])
slm_cps_kyna_trp <- summary(lm_cps_kyna_trp)

lm_cps_kyn_trp <- lm(cps_Cum~`kyn/trp1`,
                      data=df[df$`kyn/trp1`<60,])
slm_cps_kyn_trp <- summary(lm_cps_kyn_trp)

lm_cps_kyn_trp <- lm(CPS~`kyn/trp1`,
                     data=df[df$`kyn/trp1`<60,])
slm_cps_kyn_trp <- summary(lm_cps_kyn_trp)

lm_cps_kyn_trp2 <- lm(CPS~`kyn/trp2`,
                     data=df[df$`kyn/trp2`<60,])
slm_cps_kyn_trp2 <- summary(lm_cps_kyn_trp2)

lm_cps_kyn_trp_fc <- lm(CPS~`kyn/trp_fc`,data = df)
slm_cps_kyn_trp_fc <- summary(lm_cps_kyn_trp_fc)

lm_pyr_kyn_trp1 <- lm(`823.y`~`kyn/trp1`,data = df)
slm_pyr_kyn_trp1 <- summary(lm_pyr_kyn_trp1)

lm_pyr_kyn_trp_fc <- lm(`823.y`~`kyn/trp_fc`,data = df)
slm_pyr_kyn_trp_fc <- summary(lm_pyr_kyn_trp_fc)

lm_lac_kyn_trp1 <- lm(`482.y`~`kyn/trp1`,data = df)
slm_lac_kyn_trp1 <- summary(lm_lac_kyn_trp1)

lm_lac_kyn_trp_fc <- lm(`482.y`~`kyn/trp_fc`,data = df)
slm_lac_kyn_trp_fc <- summary(lm_lac_kyn_trp_fc)

lm_cps_naag1 <- lm(CPS~`100001612.x`,data=df)
slm_cps_naag1 <- summary(lm_cps_naag1)

lm_cps_naag2 <- lm(CPS~`100001612.y`,data=df)
slm_cps_naag2 <- summary(lm_cps_naag2)

lm_wcs_naag1 <- lm(`Whole Candidate Score`~`100001612.x`,data=df)
slm_wcs_naag1 <- summary(lm_wcs_naag1)

#Plot results
par(mar = c(5,5,4,2) + 0.1)
par(mfrow=c(1,1))
boxplot(CPS~`Endurance Level on the Bike- VO2 max Test`, data=df, col = c('lightpink',
                                                                          'lightblue'))
boxplot(`Whole Candidate Score`~`Endurance Level on the Bike- VO2 max Test`, data=df, col = c('lightpink',
                                                                                              'lightblue'))

wilcox.test(df$`Whole Candidate Score`[df$`Endurance Level on the Bike- VO2 max Test`
                                       == "Novice/Intermediate"],
            df$`Whole Candidate Score`[df$`Endurance Level on the Bike- VO2 max Test`
                                       == "Advanced"], alternative = "greater")

plot(x=df[,]$`100001612.x`, 
     y=as.numeric(df[,]$`Whole Candidate Score`),
     xlab = "Pre-Stress N-Acetyl-Aspartyl-Glutamate",
     ylab = "Whole Candidate Score",
     main = "WCS vs. Segment 2 NAAG",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim =,
     col = as.factor(df$`Endurance Level on the Bike- VO2 max Test`)) #,
#col = df$cps_cat+1)
text(-1.6, 7400, 
     bquote("Correlation" ==
              .(round(cor(df[,]$`Whole Candidate Score`,
                          as.numeric(df[,]$`100001612.y`), 
                          use = 
                            "pairwise.complete.obs"),4))), 
     pos = 4)
text(-1.6, 7250, bquote(R^2 == .(round(slm_wcs_naag1$r.squared,5))), pos = 4)
text(-1.6, 7100, bquote("p-value" == .(0.50)), pos = 4)
abline(lm_wcs_naag1, lty=2, col = "red", lwd = 2)


plot(x=df[,]$`100001612.y`, 
     y=as.numeric(df[,]$CPS),
     xlab = "Post-Stress N-Acetyl-Aspartyl-Glutamate",
     ylab = "Cumulative Performance Score",
     main = "CPS vs. Segment 2 NAAG",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(2,4.5)) #,
     #col = df$cps_cat+1)
text(-1.55, 4.4, 
     bquote("Correlation" ==
              .(round(cor(df[,]$CPS,
                          as.numeric(df[,]$`100001612.y`), 
                          use = 
                            "pairwise.complete.obs"),4))), 
     pos = 4)
text(-1.55, 4.2, bquote(R^2 == .(round(slm_cps_naag2$r.squared,5))), pos = 4)
text(-1.55, 4, bquote("p-value" == .(0.96)), pos = 4)
abline(lm_cps_naag2, lty=2, col = "red", lwd = 2)


plot(x=df[,]$`kyn/trp_fc`, 
     y=as.numeric(df[,]$`482.y`),
     xlab = "Kynurenine/Tryptophan FC",
     ylab = "Post-Stress Lactate",
     main = "Segment 2 Lactate vs.
    Kynurenine/Tryptophan FC",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-0.5, 1.3),
     col = df$cps_cat+1)
text(0.82, 1.25, 
     bquote("Correlation" ==
              .(round(cor(df[,]$`kyn/trp_fc`,
                          as.numeric(df[,]$`482.y`), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.82, 1.125, bquote(R^2 == .(round(slm_lac_kyn_trp_fc$r.squared,2))), pos = 4)
text(0.82, 1, bquote("p-value" == .(0.00085)), pos = 4)
abline(lm_lac_kyn_trp_fc, lty=2, col = "red", lwd = 2)


plot(x=df[,]$`kyn/trp1`, 
     y=as.numeric(df[,]$`482.y`),
     xlab = "Pre-Stress Kynurenine/Tryptophan",
     ylab = "Post-Stress Lactate",
     main = "Segment 2 Lactate vs.
    Kynurenine/Tryptophan",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-0.5, 1.5),
     col = df$cps_cat+1)
text(0.4, 1.48, 
     bquote("Correlation" ==
              .(round(cor(df[,]$`kyn/trp1`,
                          as.numeric(df[,]$`823.y`), 
                          use = 
                            "pairwise.complete.obs"),3))), 
     pos = 4)
text(0.4, 1.32, bquote(R^2 == .(round(slm_lac_kyn_trp1$r.squared,5))), pos = 4)
text(0.4, 1.16, bquote("p-value" == .(0.97)), pos = 4)
abline(lm_lac_kyn_trp1, lty=2, col = "red", lwd = 2)


plot(x=df[,]$`kyn/trp_fc`, 
     y=as.numeric(df[,]$`823.y`),
     xlab = "Kynurenine/Tryptophan FC",
     ylab = "Post-Stress Pyruvate",
     main = "Segment 2 Pyruvate vs.
    Kynurenine/Tryptophan FC",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-0.5, 1.2),
     col = df$cps_cat+1)
text(0.82, 1.15, 
     bquote("Correlation" ==
              .(round(cor(df[,]$`kyn/trp_fc`,
                          as.numeric(df[,]$`823.y`), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.82, 1.025, bquote(R^2 == .(round(slm_pyr_kyn_trp_fc$r.squared,2))), pos = 4)
text(0.82, 0.9, bquote("p-value" == .(0.0013)), pos = 4)
abline(lm_pyr_kyn_trp_fc, lty=2, col = "red", lwd = 2)


plot(x=df[,]$`kyn/trp1`, 
     y=as.numeric(df[,]$`823.y`),
     xlab = "Pre-Stress Kynurenine/Tryptophan",
     ylab = "Post-Stress Pyruvate",
     main = "Segment 2 Pyruvate vs.
    Kynurenine/Tryptophan",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-0.5, 1.2),
     col = df$cps_cat+1)
text(0.4, 1.15, 
     bquote("Correlation" ==
              .(round(cor(df[,]$`kyn/trp1`,
                          as.numeric(df[,]$`823.y`), 
                          use = 
                            "pairwise.complete.obs"),3))), 
     pos = 4)
text(0.4, 1.025, bquote(R^2 == .(round(slm_pyr_kyn_trp1$r.squared,4))), pos = 4)
text(0.4, 0.9, bquote("p-value" == .(0.81)), pos = 4)
abline(lm_pyr_kyn_trp1, lty=2, col = "red", lwd = 2)


plot(x=df[df$`kyn/trp2`<60,]$`kyn/trp_fc`, 
     y=as.numeric(df[df$`kyn/trp2`<60,]$CPS),
     xlab = "Kynurenine/Tryptophan FC",
     ylab = "Cumulative Performance Score",
     main = "Cumulative Performance Score vs. 
     Kynurenine/Tryptophan FC",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim =c(1.5, 5.2))
text(0.82, 5.1, 
     bquote("Correlation" ==
              .(round(cor(df[df$`kyn/trp2`<60,]$`kyn/trp_fc`,
                          as.numeric(df[df$`kyn/trp1`<60,]$CPS), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.82, 4.84, bquote(R^2 == .(round(slm_cps_kyn_trp_fc$r.squared,3))), pos = 4)
text(0.82, 4.58, bquote("p-value" == .(0.52)), pos = 4)
abline(lm_cps_kyn_trp_fc, lty=2, col = "red", lwd = 2)


plot(x=df[df$`kyn/trp2`<60,]$`kyn/trp2`, 
     y=as.numeric(df[df$`kyn/trp1`<60,]$CPS),
     xlab = "Post-Stress Kynurenine/Tryptophan",
     ylab = "Cumulative Performance Score",
     main = "Cumulative Performance Score vs. 
     Segment 2 Kynurenine/Tryptophan",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim =c(1.5, 5.2))
text(0.57, 5.1, 
     bquote("Correlation" ==
              .(round(cor(df[df$`kyn/trp2`<60,]$`kyn/trp2`,
                          as.numeric(df[df$`kyn/trp1`<60,]$CPS), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.57, 4.84, bquote(R^2 == .(round(slm_cps_kyn_trp2$r.squared,3))), pos = 4)
text(0.57, 4.58, bquote("p-value" == .(0.65)), pos = 4)
abline(lm_cps_kyn_trp2, lty=2, col = "red", lwd = 2)


plot(x=df[df$`kyn/trp1`<60,]$`kyn/trp1`, 
     y=as.numeric(df[df$`kyn/trp1`<60,]$cps_Cum),
     xlab = "Pre-Stress Kynurenine/Tryptophan",
     ylab = "Cognitive performance Score",
     main = "Cognitive performance Score vs. 
     Segment 2 Kynurenine/Tryptophan",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim =c(1.5, 5.2))
text(0.39, 5.1, 
     bquote("Correlation" ==
              .(round(cor(df[df$`kyn/trp1`<60,]$`kyn/trp1`,
                          as.numeric(df[df$`kyn/trp1`<60,]$cps_Cum), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.39, 4.84, bquote(R^2 == .(round(slm_cps_kyn_trp$r.squared,3))), pos = 4)
text(0.39, 4.58, bquote("p-value" == .(0.74)), pos = 4)
abline(lm_cps_kyn_trp, lty=2, col = "red", lwd = 2)


plot(x=df[df$`kyn/trp1`<60,]$`kyn/trp1`, 
     y=as.numeric(df[df$`kyn/trp1`<60,]$CPS),
     xlab = "Pre-Stress Kynurenine/Tryptophan",
     ylab = "Cumulative Performance Score",
     main = "Cumulative Performance Score vs. 
     Segment 2 Kynurenine/Tryptophan",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim =c(1.5, 5.2))
text(0.39, 5.1, 
     bquote("Correlation" ==
              .(round(cor(df[df$`kyn/trp1`<60,]$`kyn/trp1`,
                          as.numeric(df[df$`kyn/trp1`<60,]$CPS), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.39, 4.84, bquote(R^2 == .(round(slm_cps_kyn_trp$r.squared,3))), pos = 4)
text(0.39, 4.58, bquote("p-value" == .(0.48)), pos = 4)
abline(lm_cps_kyn_trp, lty=2, col = "red", lwd = 2)


plot(x=df[df$`kyn/trp1`<60,]$`kyn/trp1`, 
     y=as.numeric(df[df$`kyn/trp1`<60,]$cps_Cum),
     xlab = "Pre-Stress Kynurenate/Tryptophan",
     ylab = "Cognitive performance Score",
     main = "Cognitive performance Score vs. 
     Segment 2 Kynurenate/Tryptophan",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim =c(1.5, 5.2))
text(0.25, 5.1, 
     bquote("Correlation" ==
              .(round(cor(df[df$`kyn/trp1`<60,]$`kyn/trp1`,
                          as.numeric(df[df$`kyn/trp1`<60,]$cps_Cum), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.25, 4.84, bquote(R^2 == .(round(slm_cps_kyn_trp$r.squared,3))), pos = 4)
text(0.25, 4.58, bquote("p-value" == .(0.32)), pos = 4)
abline(lm_cps_kyn_trp, lty=2, col = "red", lwd = 2)


plot(x=df[df$`kyna/trp1`<60,]$`kyna/trp1`, 
     y=as.numeric(df[df$`kyna/trp1`<60,]$CPS),
     xlab = "Pre-Stress Kynurenate/Tryptophan",
     ylab = "Cumulative Performance Score",
     main = "Cumulative Performance Score vs. 
     Segment 2 Kynurenate/Tryptophan",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim =c(2, 5.2))
text(0.25, 5.1, 
     bquote("Correlation" ==
              .(round(cor(df[df$`kyn/trp1`<60,]$`kyna/trp1`,
                          as.numeric(df[df$`kyna/trp1`<60,]$CPS), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.25, 4.85, bquote(R^2 == .(round(slm_cps_kyna_trp$r.squared,3))), pos = 4)
text(0.25, 4.6, bquote("p-value" == .(0.40)), pos = 4)
abline(lm_cps_kyna_trp, lty=2, col = "red", lwd = 2)


plot(x=df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`trp/tyr1`, 
     y=as.numeric(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`98.y`),
     xlab = "Pre-Stress Tryptophan/Tyrosine",
     ylab = "Post-Stress Kynurenate",
     main = "Segment 2 Kynurenate vs.
     Tryptophan/Tyrosine",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim =c(-1.2,2.1))
text(0.65, 2, 
     bquote("Correlation" ==
              .(round(cor(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`trp/tyr1`,
                          as.numeric(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`98.y`), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.65, 1.74, bquote(R^2 == .(round(slm_kyn_trp_tyr$r.squared,3))), pos = 4)
text(0.65, 1.48, bquote("p-value" == .(0.0022)), pos = 4)
abline(lm_kyn_trp_tyr, lty=2, col = "red", lwd = 2)


plot(x=df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`trp/tyr1`, 
     y=as.numeric(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`1022.y`),
     xlab = "Pre-Stress Tryptophan/Tyrosine",
     ylab = "Post-Stress Picolinate",
     main = "Segment 2 Picolinate vs.
     Tryptophan/Tyrosine",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim =c(-1.2,3))
text(0.6, 2.9, 
     bquote("Correlation" ==
              .(round(cor(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`trp/tyr1`,
                          as.numeric(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`1022.y`), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.6, 2.6, bquote(R^2 == .(round(slm_pic_trp_tyr$r.squared,3))), pos = 4)
text(0.6, 2.3, bquote("p-value" == .(0.18)), pos = 4)
abline(lm_pic_trp_tyr, lty=2, col = "red", lwd = 2)


plot(x=df[df$`trp/(t + p + l + i + v)1`<25,]$`trp/(t + p + l + i + v)1`, 
     y=as.numeric(df[df$`trp/(t + p + l + i + v)1`<25,]$`98.y`),
     xlab = "Pre-Stress Tryptophan/(T+P+L+I+V)",
     ylab = "Post-Stress Kynurenate",
     main = "Segment 2 Kynurenate vs.
    Tryptophan/(T+P+L+I+V)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-1.25, 1.77))
text(0.175, 1.75, 
     bquote("Correlation" ==
              .(round(cor(df[df$`trp/(t + p + l + i + v)1`<25,]$`trp/(t + p + l + i + v)1`,
                          as.numeric(df[df$`trp/(t + p + l + i + v)1`<25,]$`98.y`), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.175, 1.55, bquote(R^2 == .(round(slm_kyn_trp_tpliv$r.squared,3))), pos = 4)
text(0.175, 1.33, bquote("p-value" == .(0.077)), pos = 4)
abline(lm_kyn_trp_tpliv, lty=2, col = "red", lwd = 2)


plot(x=df[df$`trp/(t + p + l + i + v)1`<25,]$`trp/(t + p + l + i + v)1`, 
     y=as.numeric(df[df$`trp/(t + p + l + i + v)1`<25,]$`1022.y`),
     xlab = "Pre-Stress Tryptophan/(T+P+L+I+V)",
     ylab = "Post-Stress Picolinate",
     main = "Segment 2 Picolinate vs.
    Tryptophan/(T+P+L+I+V)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-1.25, 2))
text(0.175, 1.9, 
     bquote("Correlation" ==
              .(round(cor(df[df$`trp/(t + p + l + i + v)1`<25,]$`trp/(t + p + l + i + v)1`,
                          as.numeric(df[df$`trp/(t + p + l + i + v)1`<25,]$`1022.y`), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.175, 1.65, bquote(R^2 == .(round(slm_pic_trp_tpliv$r.squared,3))), pos = 4)
text(0.175, 1.4, bquote("p-value" == .(0.070)), pos = 4)
abline(lm_pic_trp_tpliv, lty=2, col = "red", lwd = 2)

write.csv(df,paste(hgb_path,'df.csv',sep=""))

plot(x=df[df$`trp/(t + p + l + i + v)1`<25,]$`trp/(t + p + l + i + v)1`, 
     y=as.numeric(df[df$`trp/(t + p + l + i + v)1`<25,]$`482.y`),
     xlab = "Pre-Stress Tryptophan/(T+P+L+I+V)",
     ylab = "Post-Stress Lactate",
     main = "Segment 2 Lactate vs.
    Tryptophan/(T+P+L+I+V)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-0.5, 1.53))
text(0.165, 1.5, 
     bquote("Correlation" ==
              .(round(cor(df[df$`trp/(t + p + l + i + v)1`<25,]$`trp/(t + p + l + i + v)1`,
                          as.numeric(df[df$`trp/(t + p + l + i + v)1`<25,]$`482.y`), 
                          use = 
                            "pairwise.complete.obs"),3))), 
     pos = 4)
text(0.165, 1.35, bquote(R^2 == .(round(slm_lac_trp_tpliv$r.squared,4))), pos = 4)
text(0.165, 1.2, bquote("p-value" == .(0.46)), pos = 4)
abline(lm_lac_trp_tpliv, lty=2, col = "red", lwd = 2)


plot(x=df[df$`trp/(t + p + l + i + v)1`<25,]$`trp/(t + p + l + i + v)1`, 
     y=as.numeric(df[df$`trp/(t + p + l + i + v)1`<25,]$`823.y`),
     xlab = "Pre-Stress Tryptophan/(T+P+L+I+V)",
     ylab = "Post-Stress Pyruvate",
     main = "Segment 2 Pyruvate vs.
    Tryptophan/(T+P+L+I+V)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-0.5, 1.2))
text(0.165, 1.15, 
     bquote("Correlation" ==
              .(round(cor(df[df$`trp/(t + p + l + i + v)1`<25,]$`trp/(t + p + l + i + v)1`,
                          as.numeric(df[df$`trp/(t + p + l + i + v)1`<25,]$`823.y`), 
                          use = 
                            "pairwise.complete.obs"),3))), 
     pos = 4)
text(0.165, 1.025, bquote(R^2 == .(round(slm_pyr_trp_tpliv$r.squared,4))), pos = 4)
text(0.165, 0.9, bquote("p-value" == .(0.66)), pos = 4)
abline(lm_pyr_trp_tpliv, lty=2, col = "red", lwd = 2)


plot(x=df[df$`trp/(t + p + l + i + v)1`<25,]$`trp/(t + p + l + i + v)1`, 
     y=as.numeric(df[df$`trp/(t + p + l + i + v)1`<25,]$`hgb.y`),
     xlab = "Pre-Stress Tryptophan/(T+P+L+I+V)",
     ylab = "Post-Stress Hemoglobin (g/dL)",
     main = "Segment 2 Hemoglobin vs.
    Tryptophan/(T+P+L+I+V)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(12,20.1))
text(0.165, 20, 
     bquote("Correlation" ==
              .(round(cor(df[df$`trp/(t + p + l + i + v)1`<25,]$`trp/(t + p + l + i + v)1`,
                          as.numeric(df[df$`trp/(t + p + l + i + v)1`<25,]$`hgb.y`), 
                          use = 
                            "pairwise.complete.obs"),4))), 
     pos = 4)
text(0.165, 19.4, bquote(R^2 == .(round(slm_hgb_trp_tpliv$r.squared,4))), pos = 4)
text(0.165, 18.8, bquote("p-value" == .(0.94)), pos = 4)
abline(lm_hgb_trp_tpliv, lty=2, col = "red", lwd = 2)


plot(x=df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`trp/tyr1`, 
     y=as.numeric(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`hgb.y`),
     xlab = "Pre-Stress Tryptophan/Tyrosine",
     ylab = "Post-Stress Hemoglobin (g/dL)",
     main = "Segment 2 Hemoglobin vs.
     Tryptophan/Tyrosine",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim =c(13, 21))
text(0.6, 20.5, 
     bquote("Correlation" ==
              .(round(cor(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`trp/tyr1`,
                          as.numeric(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`hgb.y`), 
                          use = 
                            "pairwise.complete.obs"),3))), 
     pos = 4)
text(0.6, 19.9, bquote(R^2 == .(round(slm_hgb_trp_tyr$r.squared,4))), pos = 4)
text(0.6, 19.3, bquote("p-value" == .(0.69)), pos = 4)
abline(lm_hgb_trp_tyr, lty=2, col = "red", lwd = 2)


plot(x=df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`trp/tyr1`, 
     y=as.numeric(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`823.y`),
     xlab = "Pre-Stress Tryptophan/Tyrosine",
     ylab = "Post-Stress Pyruvate",
     main = "Segment 2 Pyruvate vs.
     Tryptophan/Tyrosine",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-0.5, 1.5))
text(0.6, 1.48, 
     bquote("Correlation" ==
              .(round(cor(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`trp/tyr1`,
                          as.numeric(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`823.y`), 
                          use = 
                            "pairwise.complete.obs"),3))), 
     pos = 4)
text(0.6, 1.33, bquote(R^2 == .(round(slm_pyr_trp_tyr$r.squared,4))), pos = 4)
text(0.6, 1.18, bquote("p-value" == .(0.79)), pos = 4)
abline(lm_pyr_trp_tyr, lty=2, col = "red", lwd = 2)


plot(x=df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`trp/tyr1`, 
     y=as.numeric(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`482.y`),
     xlab = "Pre-Stress Tryptophan/Tyrosine",
     ylab = "Post-Stress Lactate",
     main = "Segment 2 Lactate vs.
     Tryptophan/Tyrosine",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-1, 1.8))
text(0.6, 1.75, 
     bquote("Correlation" ==
              .(round(cor(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`trp/tyr1`,
                          as.numeric(df[df$`trp/tyr1`>-5 & df$`trp/tyr1`<10,]$`482.y`), 
                          use = 
                            "pairwise.complete.obs"),4))), 
     pos = 4)
text(0.6, 1.55, bquote(R^2 == .(round(slm_lac_trp_tyr$r.squared,5))), pos = 4)
text(-3.5, 1.35, bquote("p-value" == .(0.96)), pos = 4)
abline(lm_lac_trp_tyr, lty=2, col = "red", lwd = 2)


plot(x=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`glut/crea1`, 
     y=as.numeric(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$hgb.y),
     xlab = "Pre-Stress N-Acetyl Glutamine/Creatine",
     ylab = "Post-Stress Hemoglobin (g/dL)",
     main = "Segment 2 Hemoglobin vs.
     N-Acetyl Glutamine/Creatine",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(11,20))
text(0.05, 19.8, 
     bquote("Correlation" ==
              .(round(cor(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`glut/crea1`,
                          as.numeric(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$hgb.y), 
                          use = 
                            "pairwise.complete.obs"),3))), 
     pos = 4)
text(0.05, 19.1, bquote(R^2 == .(round(slm_hgb_glu_cre$r.squared,4))), pos = 4)
text(0.05, 18.3, bquote("p-value" == .(0.59)), pos = 4)
abline(lm_hgb_glu_cre, lty=2, col = "red", lwd = 2)


plot(x=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`glut/crea1`, 
     y=as.numeric(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`Skeletal Muscle Mass (lbs)`),
     xlab = "Pre-Stress N-Acetyl Glutamine/Creatine",
     ylab = "Skeletal Muscle Mass (lbs)",
     main = "Skeletal Muscle Mass vs. Segment 2 
     N-Acetyl Glutamine/Creatine",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(40,155))
text(0.1, 152, 
     bquote("Correlation" ==
              .(round(cor(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`glut/crea1`,
                          as.numeric(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`Skeletal Muscle Mass (lbs)`), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.1, 143, bquote(R^2 == .(round(slm_mm_glu_cre$r.squared,3))), pos = 4)
text(0.1, 134, bquote("p-value" == .(0.25)), pos = 4)
abline(lm_mm_glu_cre, lty=2, col = "red", lwd = 2)


plot(x=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`glut/crea1`, 
     y=as.numeric(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`Power @ VO2 max (From Pedals/Garmin)- IMPUTED`),
     xlab = "Pre-Stress N-Acetyl Glutamine/Creatine",
     ylab = "Power @ VO2Max (W)",
     main = "Power @ VO2Max vs. Segment 2 
     N-Acetyl Glutamine/Creatine",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(190, 450))
text(0.1, 445, 
     bquote("Correlation" ==
              .(round(cor(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`glut/crea1`,
                          as.numeric(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`Power @ VO2 max (From Pedals/Garmin)- IMPUTED`), 
                          use = 
                            "pairwise.complete.obs"),3))), 
     pos = 4)
text(0.1, 425, bquote(R^2 == .(round(slm_pwr_glu_cre$r.squared,4))), pos = 4)
text(0.1, 405, bquote("p-value" == .(0.49)), pos = 4)
abline(lm_pwr_glu_cre, lty=2, col = "red", lwd = 2)


plot(x=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`glut/crea1`, 
     y=as.numeric(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`482.y`),
     xlab = "Pre-Stress N-Acetyl Glutamine/Creatine",
     ylab = "Post-Stress Lactate",
     main = "Segment 2 Lactate vs. 
     N-Acetyl Glutamine/Creatine
     (log Normalized Values)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-1,2))
text(0.05, 1.9, 
     bquote("Correlation" ==
              .(round(cor(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`glut/crea1`,
                          as.numeric(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`482.y`), 
                          use = 
                            "pairwise.complete.obs"),3))), 
     pos = 4)
text(0.05, 1.65, bquote(R^2 == .(round(slm_lac_glu_cre$r.squared,4))), pos = 4)
text(0.05, 1.4, bquote("p-value" == .(0.37)), pos = 4)
abline(lm_lac_glu_cre, lty=2, col = "red", lwd = 2)


plot(x=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`glut/crea1`, 
     y=as.numeric(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`823.y`),
     xlab = "Pre-Stress N-Acetyl Glutamine/Creatine",
     ylab = "Post-Stress Pyruvate",
     main = "Segment 2 Pyruvate vs. 
     N-Acetyl Glutamine/Creatine
     (log Normalized Values)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-0.5,1.2))
text(0.1, 1.1, 
     bquote("Correlation" ==
              .(round(cor(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`glut/crea1`,
                          as.numeric(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`823.y`), 
                          use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.1, 0.95, bquote(R^2 == .(round(slm_pyr_glu_cre$r.squared,3))), pos = 4)
text(0.1, 0.8, bquote("p-value" == .(0.34)), pos = 4)
abline(lm_pyr_glu_cre, lty=2, col = "red", lwd = 2)


plot(x=df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`glut/crea1`, 
     y=as.numeric(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`VO2 Max (mL/Kg/min)`),
     xlab = "Pre-Stress N-Acetyl Glutamine/Creatine",
     ylab = "VO2Max (mL/Kg/min)",
     main = "VO2Max vs. Segment 2 
     N-Acetyl Glutamine/Creatine",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(29, 70))
text(0.05, 69, 
     bquote("Correlation" ==
              .(round(cor(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`glut/crea1`,
                          as.numeric(df[df$`glut/crea1`>-20 & df$`glut/crea1`<15,]$`VO2 Max (mL/Kg/min)`), 
                          use = 
                            "pairwise.complete.obs"),3))), 
     pos = 4)
text(0.05, 65.5, bquote(R^2 == .(round(slm_vo2_glu_cre$r.squared,4))), pos = 4)
text(0.05, 62, bquote("p-value" == .(0.33)), pos = 4)
abline(lm_vo2_glu_cre, lty=2, col = "red", lwd = 2)


plot(x=df$`1221.y`, y=df$`Power @ VO2 max (From Pedals/Garmin)- IMPUTED`,
     xlab = "Post-Stress Creatine",
     ylab = "Power @ VO2Max (W)",
     main = "Power @ VO2Max vs. Segment 2 
     Creatine (log Normalized Vals)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(190,425))
text(-0.7, 423, 
     bquote("Correlation" ==
              .(round(cor(df$`1221.y`,as.numeric(df$`Power @ VO2 max (From Pedals/Garmin)- IMPUTED`), use = 
                            "pairwise.complete.obs"),4))), 
     pos = 4)
text(-0.7, 407, bquote(R^2 == .(round(slm_pwr_cre$r.squared,6))), pos = 4)
text(-0.7, 389, bquote("p-value" == .(0.98)), pos = 4)
abline(lm_pwr_cre, lty=2, col = "red", lwd = 2)


plot(x=df$`1221.y`, y=df$`VO2 Max (mL/Kg/min)`,
     xlab = "Post-Stress Creatine",
     ylab = "VO2Max (mL/Kg/min)",
     main = "VO2Max vs. Segment 2 Creatine 
     (log Normalized Vals)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(30, 65))
text(-0.7, 64.5, 
     bquote("Correlation" ==
              .(round(cor(df$`1221.y`,as.numeric(df$`VO2 Max (mL/Kg/min)`), use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(-0.7, 61.75, bquote(R^2 == .(round(slm_vo2_cre$r.squared,4))), pos = 4)
text(-0.7, 59, bquote("p-value" == .(0.80)), pos = 4)
abline(lm_vo2_cre, lty=2, col = "red", lwd = 2)


plot(x=df$`1221.y`, y=df$`482.y`,
     xlab = "Post-Stress Creatine",
     ylab = "Post-Stress Lactate",
     main = "Segment 2 Lactate vs. Creatine 
     (log Normalized Vals)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-1, 1.75))
text(-0.7, 1.65, 
     bquote("Correlation" ==
              .(round(cor(df$`1221.y`,df$`482.y`, use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(-0.7, 1.44, bquote(R^2 == .(round(slm_lac_cre$r.squared,4))), pos = 4)
text(-0.7, 1.23, bquote("p-value" == .(0.58)), pos = 4)
abline(lm_lac_cre, lty=2, col = "red", lwd = 2)


plot(x=df$`1221.y`, y=df$`823.y`,
     xlab = "Post-Stress Creatine",
     ylab = "Post-Stress Pyruvate",
     main = "Segment 2 Pyruvate vs. Creatine 
     (log Normalized Vals)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-1, 1.75))
text(-0.7, 1.65, 
     bquote("Correlation" ==
              .(round(cor(df$`1221.y`,df$`823.y`, use = 
                            "pairwise.complete.obs"),3))), 
     pos = 4)
text(-0.7, 1.44, bquote(R^2 == .(round(slm_pyr_cre$r.squared,4))), pos = 4)
text(-0.7, 1.23, bquote("p-value" == .(0.91)), pos = 4)
abline(lm_pyr_cre, lty=2, col = "red", lwd = 2)


plot(x=df$trp_fc, y=df$`482.y`,
     xlab = "Tryptophan BD2/BD1 Fold Change",
     ylab = "Post-Stress Lactate",
     main = "Segment 2 Lactate vs. Tryptophan 
     Fold Change (log Normalized Vals)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-1, 1.75),
     col = df$cps_cat+1)
text(0.5, 1.65, 
     bquote("Correlation" ==
              .(round(cor(df$trp_fc,df$`482.y`, use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.5, 1.44, bquote(R^2 == .(round(slm_lac_trp_fc$r.squared,3))), pos = 4)
text(0.5, 1.23, bquote("p-value" == .(0.026)), pos = 4)
abline(lm_lac_trp_fc, lty=2, col = "red", lwd = 2)

h_lac_trp1 <- hist(df$trp_fc[df$anaerobic_cat==0])
h_lac_trp2 <- hist(df$trp_fc[df$anaerobic_cat==1])
plot(h_lac_trp1, col = rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue"), 
     xlim = c(0.4, 1.4), ylim = c(0, 15))
plot(h_lac_trp2, col = rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink"), add = T)

h_lac_trp1 <- hist(df$`482.y`[df$anaerobic_cat==0])
h_lac_trp2 <- hist(df$`482.y`[df$anaerobic_cat==1])
plot(h_lac_trp1, col = rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue"), 
     xlim = c(-2, 2), ylim = c(0, 20))
plot(h_lac_trp2, col = rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink"), add = T)


wilcox.test(df$trp_fc[df$anaerobic_cat==0], df$trp_fc[df$anaerobic_cat==1], paired = F)

plot(x=df$trp_fc, y=df$`823.y`,
     xlab = "Tryptophan BD2/BD1 Fold Change",
     ylab = "Post-Stress Pyruvate",
     main = "Segment 2 Pyruvate vs. Tryptophan 
     Fold Change (log Normalized Vals)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-1, 1.75))
text(0.5, 1.61, 
     bquote("Correlation" ==
              .(round(cor(df$trp_fc,df$`823.y`, use = 
                            "pairwise.complete.obs"),3))), 
     pos = 4)
text(0.5, 1.4, bquote(R^2 == .(round(slm_pyr_trp_fc$r.squared,3))), pos = 4)
text(0.5, 1.19, bquote("p-value" == .(0.53)), pos = 4)
abline(lm_pyr_trp_fc, lty=2, col = "red", lwd = 2)


plot(x=df$`815.y`, y=df$`823.y`,
     xlab = "Post-Stress Tyrosine",
     ylab = "Post-Stress Pyruvate",
     main = "Segment 2 Pyruvate vs. Tyrosine 
     (log of Normalized Values)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-1, 1.75),
     col = df$cps_cat + 1)
text(-0.95, 1.61, 
     bquote("Correlation" ==
              .(round(cor(df$`815.y`,df$`823.y`, use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(-0.95, 1.4, bquote(R^2 == .(round(slm_pyr_tyr$r.squared,2))), pos = 4)
text(-0.95, 1.19, bquote("p-value" == .(0.003)), pos = 4)
abline(lm_pyr_tyr, lty=2, col = "red", lwd = 2)


plot(x=df$`815.x`, y=df$`482.y`,
     xlab = "Pre-Stress Tyrosine",
     ylab = "Post-Stress Lactate",
     main = "Segment 2 Lactate vs. Tyrosine 
     (log of Normalized Values)",
     cex.lab = 1.4,
     cex.main = 1.4,
     ylim = c(-1, 1.75))
text(-0.6, 1.61, 
     bquote("Correlation" ==
              .(round(cor(df$`815.x`,df$`482.y`, use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(-0.6, 1.4, bquote(R^2 == .(round(slm_lac_tyr$r.squared,3))), pos = 4)
text(-0.6, 1.19, bquote("p-value" == .(0.38)), pos = 4)
abline(lm_lac_tyr, lty=2, col = "red", lwd = 2)


plot(x=df$`tyr_fc`, y=df$hgb.x,
     xlab = "Tyrosine BD2/BD1 Fold Change",
     ylab = "BD1 Hemoglobin (g/dL)",
     main = "Segment 2 BD1 Hemoglobin vs.
     Tyrosine BD2/BD1 Fold Change",
     cex.lab = 1.4,
     cex.main = 1.4)
text(0.41, 16.9, 
     bquote("Correlation" ==
              .(round(cor(df$`tyr_fc`,df$hgb.y, use = 
                            "pairwise.complete.obs"),2))), 
     pos = 4)
text(0.41, 16.4, bquote(R^2 == .(round(slm_hgb_tyr$r.squared,3))), pos = 4)
text(0.41, 15.9, bquote("p-value" == .(0.39)), pos = 4)
abline(lm_hgb_tyr_fc, lty=2, col = "red", lwd = 2)