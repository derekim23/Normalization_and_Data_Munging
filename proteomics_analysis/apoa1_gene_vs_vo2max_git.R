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

#Draw in VO2Max and APOA1 data
vo2path <- "/Volumes/GoogleDrive-115111199924997198421/My Drive/Objectives/Code/"
vo2 <- read_excel(paste(vo2path,"DARPA-MBA-VO2max-InBody-Spring21.xlsx",sep=""),sheet="Basic InBody VO2 max")
apoa1path <- "/Users/Derek/omics-integrator-master/notebooks/sample_data/proteomics/"

#There are multiple blood draws for the APOA1 data
bd1apoa1 <- read.csv(paste(apoa1path, "s2bd1.csv", sep = ""))[,c("subject_id","APOA1")]
bd2apoa1 <- read.csv(paste(apoa1path, "s2bd2.csv", sep = ""))[,c("subject_id","APOA1")]
bd3apoa1 <- read.csv(paste(apoa1path, "s2bd3.csv", sep = ""))[,c("subject_id","APOA1")]

#Merge the datasets
vo2bd1apoa1 <- merge(bd1apoa1, vo2[,c("subject ID","Gender","VO2 Max (mL/Kg/min)",
                                      "VO2 Rank (Bin)")],
                  by.x = "subject_id", by.y = "subject ID")
vo2bd2apoa1 <- merge(bd2apoa1, vo2[,c("subject ID","Gender","VO2 Max (mL/Kg/min)",
                                      "VO2 Rank (Bin)")],
                     by.x = "subject_id", by.y = "subject ID")
vo2bd3apoa1 <- merge(bd3apoa1, vo2[,c("subject ID","Gender","VO2 Max (mL/Kg/min)",
                                      "VO2 Rank (Bin)")],
                     by.x = "subject_id", by.y = "subject ID")

#Order data by subjects
vo2bd1apoa1 <- vo2bd1apoa1[order(vo2bd1apoa1$subject_id),]
vo2bd2apoa1 <- vo2bd2apoa1[order(vo2bd2apoa1$subject_id),]
vo2bd3apoa1 <- vo2bd3apoa1[order(vo2bd3apoa1$subject_id),]

#Take z-scores of APOA1 abundance
vo2bd1apoa1$APOA1 <- scale(as.numeric(vo2bd1apoa1$APOA1))
vo2bd2apoa1$APOA1 <- scale(as.numeric(vo2bd2apoa1$APOA1))
vo2bd3apoa1$APOA1 <- scale(as.numeric(vo2bd3apoa1$APOA1))

#Change variable type just in case
vo2bd1apoa1$`VO2 Max (mL/Kg/min)` <- as.numeric(vo2bd1apoa1$`VO2 Max (mL/Kg/min)`)
vo2bd2apoa1$`VO2 Max (mL/Kg/min)` <- as.numeric(vo2bd2apoa1$`VO2 Max (mL/Kg/min)`)
vo2bd3apoa1$`VO2 Max (mL/Kg/min)` <- as.numeric(vo2bd3apoa1$`VO2 Max (mL/Kg/min)`)

#Try ratio of APOA1 between BD1 and BD2
vo2bd2v1apoa1 <- vo2bd2apoa1[vo2bd2apoa1$subject_id%in%
                               vo2bd1apoa1$subject_id,]
vo2bd2v1apoa1$APOA1 <- vo2bd2v1apoa1$APOA1/
  vo2bd1apoa1$APOA1[vo2bd1apoa1$subject_id%in%
                      vo2bd2apoa1$subject_id]

#Scale the new ratio data
vo2bd2v1apoa1$APOA1 <- scale(vo2bd2v1apoa1$APOA1)

#Plots for Segment 2 Blood Draw 1
par(mar = c(5,5,4,2) + 0.1)
plot(vo2bd2v1apoa1[vo2bd2v1apoa1$`VO2 Rank (Bin)`!="a"
                   & vo2bd2v1apoa1$APOA1 < 1,]$`VO2 Max (mL/Kg/min)`, 
     vo2bd2v1apoa1[vo2bd2v1apoa1$`VO2 Rank (Bin)`!="a"
                   & vo2bd2v1apoa1$APOA1 < 1,]$APOA1,
     main = "Segment 2 BD2/BD1
     APOA1 Abundance Ratio Z-Score vs. VO2Max (mL/Kg/min)",
     ylab = "APOA1 Abundance Ratio Z-Score",
     xlab = "VO2Max (mL/Kg/min)")
lm2v1 <- lm(APOA1~`VO2 Max (mL/Kg/min)`,data=
              vo2bd2v1apoa1[vo2bd2v1apoa1$`VO2 Rank (Bin)`!="a"
                            & vo2bd2v1apoa1$APOA1 < 1,])
abline(lm2v1)
text(32, 0.19, bquote(R^2 == .(round(summary(lm2v1)$r.squared,3))))

plot(vo2bd1apoa1[vo2bd1apoa1$`VO2 Rank (Bin)`!="a",]$`VO2 Max (mL/Kg/min)`, 
     vo2bd1apoa1[vo2bd1apoa1$`VO2 Rank (Bin)`!="a",]$APOA1,
     main = "Segment 2 Blood Draw 1
     APOA1 Abundance Z-Score vs. VO2Max (mL/Kg/min)",
     ylab = "APOA1 Abundance Z-Score",
     xlab = "VO2Max (mL/Kg/min)")
abline(lm1)
summary(lm1)
slm1 <- summary(lm1)
text(33,1.5, bquote( R^2 == .(round(slm1$r.squared,3))))
text(32,1.5, bquote(R^2 == .(round(slm1$adj.r.squared,3))))

plot(vo2bd1apoa1[vo2bd1apoa1$`VO2 Rank (Bin)`!="a",]$`VO2 Max (mL/Kg/min)`, 
     vo2bd1apoa1[vo2bd1apoa1$`VO2 Rank (Bin)`!="a",]$APOA1,
     main = "Segment 2 Blood Draw 1
     APOA1 Abundance Z-Score vs. VO2Max (mL/Kg/min)",
     ylab = "APOA1 Abundance Z-Score",
     xlab = "VO2Max (mL/Kg/min)")

plot(vo2bd1apoa1[vo2bd1apoa1$Gender=="M" & vo2bd1apoa1$`VO2 Rank (Bin)`=="Excellent",]$`VO2 Max (mL/Kg/min)`, 
     vo2bd1apoa1[vo2bd1apoa1$Gender=="M" & vo2bd1apoa1$`VO2 Rank (Bin)`=="Excellent",]$APOA1,
     main = "Segment 2 Blood Draw 1
     APOA1 Abundance Z-Score vs. VO2Max (mL/Kg/min)
     Males Only",
     ylab = "APOA1 Abundance Z-Score",
     xlab = "VO2Max (mL/Kg/min)")

plot(vo2bd1apoa1[vo2bd1apoa1$Gender=="F" & vo2bd1apoa1$`VO2 Rank (Bin)`=="Excellent",]$`VO2 Max (mL/Kg/min)`, 
     vo2bd1apoa1[vo2bd1apoa1$Gender=="F" & vo2bd1apoa1$`VO2 Rank (Bin)`=="Excellent",]$APOA1,
     main = "Segment 2 Blood Draw 1
     APOA1 Abundance Z-Score vs. VO2Max (mL/Kg/min)
     Females Only",
     ylab = "APOA1 Abundance Z-Score",
     xlab = "VO2Max (mL/Kg/min)")

#Plots for Segment 2 Blood Draw 2
plot(vo2bd2apoa1[vo2bd2apoa1$`VO2 Rank (Bin)`!="Excellent",]$`VO2 Max (mL/Kg/min)`, 
     vo2bd2apoa1[vo2bd2apoa1$`VO2 Rank (Bin)`!="Excellent",]$APOA1,
     main = "Segment 2 Blood Draw 2
     APOA1 Abundance Z-Score vs. VO2Max (mL/Kg/min)",
     ylab = "APOA1 Abundance Z-Score",
     xlab = "VO2Max (mL/Kg/min)")
abline(lm2)
summary(lm2)
slm2 <- summary(lm2)
text(31,1.75, bquote(R^2 == .(round(slm2$r.squared,3))))
text(31,1.25, bquote(F_stat == .(round(slm2$fstatistic,3))))

plot(vo2bd2apoa1[vo2bd2apoa1$Gender=="M" & vo2bd2apoa1$`VO2 Rank (Bin)`=="Excellent",]$`VO2 Max (mL/Kg/min)`, 
     vo2bd2apoa1[vo2bd2apoa1$Gender=="M" & vo2bd2apoa1$`VO2 Rank (Bin)`=="Excellent",]$APOA1,
     main = "Segment 2 Blood Draw 2
     APOA1 Abundance Z-Score vs. VO2Max (mL/Kg/min)
     Males Only",
     ylab = "APOA1 Abundance Z-Score",
     xlab = "VO2Max (mL/Kg/min)")

plot(vo2bd2apoa1[vo2bd2apoa1$Gender=="F" & vo2bd2apoa1$`VO2 Rank (Bin)`=="Excellent",]$`VO2 Max (mL/Kg/min)`, 
     vo2bd2apoa1[vo2bd2apoa1$Gender=="F" & vo2bd2apoa1$`VO2 Rank (Bin)`=="Excellent",]$APOA1,
     main = "Segment 2 Blood Draw 2
     APOA1 Abundance Z-Score vs. VO2Max (mL/Kg/min)
     Females Only",
     ylab = "APOA1 Abundance Z-Score",
     xlab = "VO2Max (mL/Kg/min)")

#Plots for Segment 2 Blood Draw 3
plot(vo2bd3apoa1[vo2bd3apoa1$`VO2 Rank (Bin)`!="Superior",]$`VO2 Max (mL/Kg/min)`, 
     vo2bd3apoa1[vo2bd3apoa1$`VO2 Rank (Bin)`!="Superior",]$APOA1,
     main = "Segment 2 Blood Draw 3
     APOA1 Abundance Z-Score vs. VO2Max (mL/Kg/min)",
     ylab = "APOA1 Abundance Z-Score",
     xlab = "VO2Max (mL/Kg/min)")
abline(lm3)
summary(lm3)
slm3 <- summary(lm3)
text(33,1.75, bquote(R^2 == .(round(slm3$r.squared,3))))

plot(vo2bd3apoa1[vo2bd3apoa1$Gender=="M" & vo2bd3apoa1$`VO2 Rank (Bin)`=="Excellent",]$`VO2 Max (mL/Kg/min)`, 
     vo2bd3apoa1[vo2bd3apoa1$Gender=="M" & vo2bd3apoa1$`VO2 Rank (Bin)`=="Excellent",]$APOA1,
     main = "Segment 2 Blood Draw 3
     APOA1 Abundance Z-Score vs. VO2Max (mL/Kg/min)
     Males Only",
     ylab = "APOA1 Abundance Z-Score",
     xlab = "VO2Max (mL/Kg/min)")

plot(vo2bd3apoa1[vo2bd3apoa1$Gender=="F" & vo2bd3apoa1$`VO2 Rank (Bin)`=="Excellent",]$`VO2 Max (mL/Kg/min)`, 
     vo2bd3apoa1[vo2bd3apoa1$Gender=="F" & vo2bd3apoa1$`VO2 Rank (Bin)`=="Excellent",]$APOA1,
     main = "Segment 2 Blood Draw 3
     APOA1 Abundance Z-Score vs. VO2Max (mL/Kg/min)
     Females Only",
     ylab = "APOA1 Abundance Z-Score",
     xlab = "VO2Max (mL/Kg/min)")

#Fit linear models by BDs
#Try also controlling for biological gender
lm1 <- lm(APOA1~`VO2 Max (mL/Kg/min)`,data = vo2bd1apoa1[,])
lm1m <- lm(APOA1~`VO2 Max (mL/Kg/min)`,data = vo2bd1apoa1[vo2bd1apoa1$Gender=="M",])
lm1f <- lm(APOA1~`VO2 Max (mL/Kg/min)`,data = vo2bd1apoa1[vo2bd1apoa1$Gender=="F",])

lm2 <- lm(APOA1~`VO2 Max (mL/Kg/min)`,data = vo2bd2apoa1[,])
lm2m <- lm(APOA1~`VO2 Max (mL/Kg/min)`,data = vo2bd2apoa1[vo2bd2apoa1$Gender=="M",])
lm2f <- lm(APOA1~`VO2 Max (mL/Kg/min)`,data = vo2bd2apoa1[vo2bd2apoa1$Gender=="F",])

lm3 <- lm(APOA1~`VO2 Max (mL/Kg/min)`,data = vo2bd3apoa1[,])
lm3m <- lm(APOA1~`VO2 Max (mL/Kg/min)`,data = vo2bd3apoa1[vo2bd3apoa1$Gender=="M",])
lm3f <- lm(APOA1~`VO2 Max (mL/Kg/min)`,data = vo2bd3apoa1[vo2bd3apoa1$Gender=="F",])