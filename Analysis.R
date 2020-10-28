#import data
genedata_orig = read.csv(file.choose())
clinicaldata_orig = read.csv(file.choose())

#cleaning packages
install.packages("tidyverse", dep = T)
library(tidyverse)

#gene data cleaning
genedata_clean = as.data.frame(t(genedata_orig), header = T)
names(genedata_clean) = genedata_clean %>% slice(1) %>% unlist()
genedata_clean = genedata_clean %>% slice(-1)
genedata_clean = mutate_all(genedata_clean, function(x) as.numeric(as.character(x)))

#scaling the gene data
genedata_scaled = as.data.frame(scale(genedata_clean, scale = T, center =  T))

#clinical data cleaning
clinicaldata_clean = clinicaldata_orig
clinicaldata_clean = select(clinicaldata_clean,-c(ID))
clinicaldata_clean = select(clinicaldata_clean,-c(X.Sample_characteristics_ch1 ,X.Sample_characteristics_ch1.14))
clinicaldata_clean = as_tibble(clinicaldata_clean) %>% rename(Gender = X.Sample_characteristics_ch1.1, 
                                                              Age = X.Sample_characteristics_ch1.2,
                                                              Ethnicity = X.Sample_characteristics_ch1.3, 
                                                              Vital_Status = X.Sample_characteristics_ch1.4,
                                                              Adjuvant_Chemo = X.Sample_characteristics_ch1.5, 
                                                              Adjuvant_RT = X.Sample_characteristics_ch1.6,
                                                              Disease_stage = X.Sample_characteristics_ch1.7, 
                                                              Relapse = X.Sample_characteristics_ch1.8,
                                                              PFS = X.Sample_characteristics_ch1.9, 
                                                              Months_to_last_cli_assess = X.Sample_characteristics_ch1.10,
                                                              OS =  X.Sample_characteristics_ch1.11, 
                                                              Smoking_history = X.Sample_characteristics_ch1.12,
                                                              Surgical_Margin = X.Sample_characteristics_ch1.13, 
                                                              Histological_grade =X.Sample_characteristics_ch1.15)

library(tidyr)
clinicaldata_clean = clinicaldata_clean %>% separate(Gender,c("Old1","Gender"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(Age,c("Old2","Age"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(Ethnicity,c("Old3","Ethnicity"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(Vital_Status,c("Old4","Vital_Status"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(Adjuvant_Chemo,c("Old5","Adjuvant_Chemo"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(Adjuvant_RT,c("Old6","Adjuvant_RT"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(Disease_stage,c("Old7","Disease_stage"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(Relapse,c("Old8","Relapse"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(PFS,c("Old9","PFS"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(Months_to_last_cli_assess,c("Old10","Months_to_last_cli_assess"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(OS,c("Old11","OS"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(Smoking_history,c("Old12","Smoking_history"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(Surgical_Margin,c("Old13","Surgical_Margin"),": ")
clinicaldata_clean = clinicaldata_clean %>% separate(Histological_grade,c("Old14","Histological_grade"),": ")
clinicaldata_clean = select(clinicaldata_clean,-c(Old1,Old2,Old3,Old4,Old5,Old6,Old7,Old8,Old9,Old10,Old11,Old12,Old13,Old14,Months_to_last_cli_assess))

clinicaldata_clean = clinicaldata_clean %>% mutate(Relapse, Relapse = ifelse(Relapse != "Yes","No","Yes"))
clinicaldata_clean = clinicaldata_clean %>% mutate(Smoking_history, Smoking_history = ifelse(Smoking_history == "--","Never smoked",ifelse(Smoking_history == "Unknown","Never smoked",Smoking_history)))
clinicaldata_clean = clinicaldata_clean %>% mutate(Surgical_Margin, Surgical_Margin = ifelse(Surgical_Margin == "--","Unknown",Surgical_Margin))
clinicaldata_clean = clinicaldata_clean %>% mutate(Histological_grade, Histological_grade = ifelse(Histological_grade == "--","Unknown",Histological_grade))

#creating a competing event variable
clinicaldata_clean$Event = NA
clinicaldata_clean$Event[clinicaldata_clean$Relapse == "No" & clinicaldata_clean$Vital_Status == "Alive" ] ="Alive"
clinicaldata_clean$Event[clinicaldata_clean$Relapse == "No" & clinicaldata_clean$Vital_Status == "Dead" ] ="Competing Risk"
clinicaldata_clean$Event[clinicaldata_clean$Relapse == "Yes" & clinicaldata_clean$Vital_Status == "Alive" ] ="Alive"
clinicaldata_clean$Event[clinicaldata_clean$Relapse == "Yes" & clinicaldata_clean$Vital_Status == "Dead" ] ="Cancer"

#final touch up cleanings of the clinical data
clinicaldata_clean$Age = as.numeric(clinicaldata_clean$Age)
clinicaldata_clean$PFS = as.numeric(clinicaldata_clean$PFS) 
clinicaldata_clean$OS = as.numeric(clinicaldata_clean$OS) 
clinicaldata_clean = clinicaldata_clean %>% mutate(OS, OS = ifelse(is.na(OS) == T,204.00,OS))
clinicaldata_clean = clinicaldata_clean %>% mutate(PFS, PFS = ifelse(is.na(PFS)==T,OS,PFS))
clinicaldata_clean$PFS = as.numeric(clinicaldata_clean$PFS) 
clinicaldata_clean$OS = as.numeric(clinicaldata_clean$OS) 

clinicaldata_clean$Gender = as.factor(clinicaldata_clean$Gender)
clinicaldata_clean$Ethnicity = as.factor(clinicaldata_clean$Ethnicity)
clinicaldata_clean$Vital_Status = as.factor(clinicaldata_clean$Vital_Status)
clinicaldata_clean$Adjuvant_Chemo = as.factor(clinicaldata_clean$Adjuvant_Chemo)
clinicaldata_clean$Adjuvant_RT = as.factor(clinicaldata_clean$Adjuvant_RT)
clinicaldata_clean$Disease_stage = as.factor(clinicaldata_clean$Disease_stage)
clinicaldata_clean$Relapse = as.factor(clinicaldata_clean$Relapse)
clinicaldata_clean$Smoking_history = as.factor(clinicaldata_clean$Smoking_history)
clinicaldata_clean$Surgical_Margin = as.factor(clinicaldata_clean$Surgical_Margin)
clinicaldata_clean$Histological_grade = as.factor(clinicaldata_clean$Histological_grade)
clinicaldata_clean$Event = as.factor(clinicaldata_clean$Event)
clinicaldata_clean$Relapse = as.numeric(c("No" = "0", "Yes" = "1")[clinicaldata_clean$Relapse])
clinicaldata_clean$Vital_Status = as.numeric(c("Alive" = "0", "Dead" = "1")[clinicaldata_clean$Vital_Status])

#deleting one patient to make the number of patients non-prime
clinicaldata_clean = clinicaldata_clean[-53,]
genedata_clean = genedata_clean[-53,]
genedata_scaled = genedata_scaled[-53,]

#merging both data sets
clinicaldata_clean$ID = seq(1:442)
genedata_clean$ID = seq(1:442)
genedata_scaled$ID = seq(1:442)
mergeddata = merge(clinicaldata_clean,genedata_clean, by = "ID")
mergedata_scaled = merge(clinicaldata_clean,genedata_scaled, by = "ID")

mergeddata = mergeddata[,-1]
mergedata_scaled = mergedata_scaled[,-1]
clinicaldata_clean = clinicaldata_clean[,-1]
genedata_clean = genedata_clean[,-1]

#FDR
install.packages("fuzzySim", dependencies = T)
library(fuzzySim)
fdr.fit = FDR(data = mergedata_scaled, sp.cols = 4, var.cols = 15:22229,family = "auto", correction = "fdr", q = 0.05, verbose = TRUE, simplif = FALSE)

hist(fdr.fit$exclude$p.value)
hist(fdr.fit$select$p.value)

#make subsets of data as per the event
cancerdata = subset(mergeddata,mergeddata$Event != "Competing Risk")
cancerdata$Event = as.numeric(c("Alive" = "0", "Cancer" = "1")[cancerdata$Event])
cmprskdata = subset(mergeddata,mergeddata$Event != "Cancer")
cmprskdata$Event = as.numeric(c("Alive" = "0", "Competing Risk" = "1")[cmprskdata$Event])

cancerdata_scaled = subset(mergedata_scaled,mergedata_scaled$Event != "Competing Risk")
cancerdata_scaled$Event = as.numeric(c("Alive" = "0", "Cancer" = "1")[cancerdata_scaled$Event])
cmprskdata_scaled = subset(mergedata_scaled,mergedata_scaled$Event != "Cancer")
cmprskdata_scaled$Event = as.numeric(c("Alive" = "0", "Competing Risk" = "1")[cmprskdata_scaled$Event])

#define survival objects
library(survival)
so1 = Surv(mergeddata$PFS,mergeddata$Relapse)
so2 = Surv(mergeddata$OS,mergeddata$Vital_Status)
so3 = Surv(cancerdata$OS,cancerdata$Event)
so4 = Surv(cmprskdata$OS,cmprskdata$Event)

#genematrix
gmatrix = data.matrix(mergeddata[,15:22229])
gmatrixe1 = data.matrix(cancerdata[,15:22229])
gmatrixe2 = data.matrix(cmprskdata[,15:22229])

s_gmatrix = data.matrix(mergedata_scaled[,15:22229])
s_gmatrixe1 = data.matrix(cancerdata_scaled[,15:22229])
s_gmatrixe2 = data.matrix(cmprskdata_scaled[,15:22229])

#Variable selections under varying survival scenarios
install.packages("glmnet", dep = T)
library(glmnet)

#setseed function for exact replication later
#set.seed(999)
#creating foldid for replication feasibility
#foldid = sample(1:15,size = length(so1),replace = T)

#lasso ss1
lr1 = cv.glmnet(gmatrix, so1, alpha=1, family = "cox", nfolds = 10)

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/lr1.1.pdf", width = 15, height = 10) 
plot(lr1)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/lr1.2.pdf", width = 15, height = 10)
plot(lr1$glmnet.fit, "lambda", label = T)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/lr1.3.pdf", width = 15, height = 10) 
plot(lr1$glmnet.fit, "norm", label = T)
dev.off()

lr1$lambda.min
log(lr1$lambda.min)
lr1c = coef(lr1, s = lr1$lambda.min)
lr1lam = which(lr1c!=0)
lr1lam
lr1lamval = lr1c[lr1lam]
lr1lamval

#lasso ss2
lr2 = cv.glmnet(gmatrix, so2, alpha=1, family = "cox", nfolds = 10)

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/lr2.1.pdf", width = 15, height = 10) 
plot(lr2)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/lr2.2.pdf", width = 15, height = 10)
plot(lr2$glmnet.fit, "lambda", label = T)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/lr2.3.pdf", width = 15, height = 10) 
plot(lr2$glmnet.fit, "norm", label = T)
dev.off()

lr2$lambda.min
log(lr2$lambda.min)
lr2c = coef(lr2, s = lr2$lambda.min)
lr2lam = which(lr2c!=0)
lr2lam
lr2lamval = lr2c[lr2lam]
lr2lamval

#lasso ss3
lr3 = cv.glmnet(gmatrixe1, so3, alpha=1, family = "cox", nfolds = 10)

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/lr3.1.pdf", width = 15, height = 10) 
plot(lr3)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/lr3.2.pdf", width = 15, height = 10)
plot(lr3$glmnet.fit, "lambda", label = T)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/lr3.3.pdf", width = 15, height = 10) 
plot(lr3$glmnet.fit, "norm", label = T)
dev.off()

lr3$lambda.min
log(lr3$lambda.min)
lr3c = coef(lr3, s = lr3$lambda.min)
lr3lam = which(lr3c!=0)
lr3lam
lr3lamval = lr3c[lr3lam]
lr3lamval

#lasso ss4
lr4 = cv.glmnet(gmatrixe2, so4, alpha=1, family = "cox", nfolds = 10)

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/lr4.1.pdf", width = 15, height = 10) 
plot(lr4)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/lr4.2.pdf", width = 15, height = 10)
plot(lr4$glmnet.fit, "lambda", label = T)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/lr4.3.pdf", width = 15, height = 10) 
plot(lr4$glmnet.fit, "norm", label = T)
dev.off()

lr4$lambda.min
log(lr4$lambda.min)
lr4c = coef(lr4, s = lr4$lambda.min)
lr4lam = which(lr4c!=0)
lr4lam
lr4lamval = lr4c[lr4lam]
lr4lamval

#enet ss1
er1 = cv.glmnet(gmatrix, so1, alpha=0.5, family = "cox", nfolds = 10)

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/er1.1.pdf", width = 15, height = 10) 
plot(er1)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/er1.2.pdf", width = 15, height = 10)
plot(er1$glmnet.fit, "lambda", label = T)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/er1.3.pdf", width = 15, height = 10) 
plot(er1$glmnet.fit, "norm", label = T)
dev.off()

er1$lambda.min
log(er1$lambda.min)
er1c = coef(er1, s = er1$lambda.min)
er1lam = which(er1c!=0)
er1lam
er1lamval = er1c[er1lam]
er1lamval

#enet ss2
er2 = cv.glmnet(gmatrix, so2,  alpha=0.5, family = "cox", nfolds = 10)

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/er2.1.pdf", width = 15, height = 10) 
plot(er2)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/er2.2.pdf", width = 15, height = 10)
plot(er2$glmnet.fit, "lambda", label = T)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/er2.3.pdf", width = 15, height = 10) 
plot(er2$glmnet.fit, "norm", label = T)
dev.off()

er2$lambda.min
log(er2$lambda.min)
er2c = coef(er2, s = er2$lambda.min)
er2lam = which(er2c!=0)
er2lam
er2lamval = er2c[er2lam]
er2lamval

#enet ss3
er3 = cv.glmnet(gmatrixe1, so3,  alpha=0.5, family = "cox", nfolds = 10)

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/er3.1.pdf", width = 15, height = 10) 
plot(er3)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/er3.2.pdf", width = 15, height = 10)
plot(er3$glmnet.fit, "lambda", label = T)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/er3.3.pdf", width = 15, height = 10) 
plot(er3$glmnet.fit, "norm", label = T)
dev.off()

er3$lambda.min
log(er3$lambda.min)
er3c = coef(er3, s = er3$lambda.min)
er3lam = which(er3c!=0)
er3lam
er3lamval = er3c[er3lam]
er3lamval

#enet ss4
er4 = cv.glmnet(gmatrixe2, so4, alpha=0.5, family = "cox", nfolds = 10)

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/er4.1.pdf", width = 15, height = 10) 
plot(er4)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/er4.2.pdf", width = 15, height = 10)
plot(er4$glmnet.fit, "lambda", label = T)
dev.off()

pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/er4.3.pdf", width = 15, height = 10) 
plot(er4$glmnet.fit, "norm", label = T)
dev.off()

er4$lambda.min
log(er4$lambda.min)
er4c = coef(er4, s = er4$lambda.min)
er4lam = which(er4c!=0)
er4lam
er4lamval = er4c[er4lam]
er4lamval

#qpcR package helps cbind unequal length vectors
install.packages("qpcR", dep = T)
library(qpcR)

#genes selected for survival scenario 1
tempdf = as.data.frame(qpcR:::cbind.na(lr1lam,er1lam))
allgene = stack(tempdf)
allgene = allgene[,-2]
ss1gene_index = unique(allgene)
ss1gene_index = na.omit(ss1gene_index)
rm(tempdf,allgene)

#genes selected for survival scenario 2
tempdf1 = as.data.frame(qpcR:::cbind.na(lr2lam,er2lam))
allgene1 = stack(tempdf1)
allgene1 = allgene1[,-2]
ss2gene_index = unique(allgene1)
ss2gene_index = na.omit(ss2gene_index)
rm(tempdf1,allgene1)

#genes selected for survival scenario 3
tempdf2 = as.data.frame(qpcR:::cbind.na(lr3lam,er3lam))
allgene2 = stack(tempdf2)
allgene2 = allgene2[,-2]
ss3gene_index = unique(allgene2)
ss3gene_index = na.omit(ss3gene_index)
rm(tempdf2,allgene2)

#genes selected for survival scenario 4
tempdf3 = as.data.frame(qpcR:::cbind.na(lr4lam,er4lam))
allgene3 = stack(tempdf3)
allgene3 = allgene3[,-2]
ss4gene_index = unique(allgene3)
ss4gene_index = na.omit(ss4gene_index)
rm(tempdf3,allgene3)

#s1 df, s2 df, s3 df and s4 df
s1df = as.data.frame(gmatrix[,ss1gene_index])
s2df = as.data.frame(gmatrix[,ss2gene_index])
s3df = as.data.frame(gmatrixe1[,ss3gene_index])
s4df = as.data.frame(gmatrixe2[,ss4gene_index])

s1df = cbind(clinicaldata_clean$Relapse,clinicaldata_clean$PFS,s1df)
colnames(s1df) = c("Relapse","PFS","GAPDHS","PSIP1","GM2A","DAB1","LRFN4","TRIM45","MLNR")

s2df = cbind(clinicaldata_clean$Vital_Status,clinicaldata_clean$OS,s2df)
colnames(s2df) = c("Vital_status","OS","NDST1","GALNT3","PDPK1","WDHD1","ZC2HC1A"
                   ,"XPNPEP1","SWAP70","MFHAS1","HILPDA","FAM117A")

s3df = cbind(cancerdata$Event,cancerdata$OS,s3df)
colnames(s3df) = c("Event","OS","GALNT3","XPNPEP1","SWAP70","PSIP1","MRPL9","NUS1P3",
                   "RCAN1","ITGB1","IL23A","PIAS1","HILPDA","AFF4","LRFN4","TESMIN","DAB1",
                   "FAM117A","H2AW","STC2")

s4df = cbind(cmprskdata$Event,cmprskdata$OS,s4df)
colnames(s4df) = c("Event","OS","TRIP12","BRD4","CHMP2B","NEBL","STEAP1","SLK","SOD2","VPS51",
                   "ERO1A","PTER","IFIH1","C1GALT1C1","STYK1","PPARD","GALNT3","UBA3","CCDC93","MIIP")

library(survival)
install.packages("survminer", dep = T)
library(survminer)
library(dyplr)
#coxregression on ss1
temp = s1df
temp$Relapse = as.numeric(temp$Relapse)
cr1 = coxph(Surv(PFS, Relapse) ~ GAPDHS+PSIP1+GM2A+DAB1+LRFN4+TRIM45+MLNR, data = temp)
summary(cr1)
pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/forrest1.pdf", width = 12, height = 12) 
ggforest(cr1)
dev.off()
rm(temp)

#coxregression on ss2
temp = s2df
temp$Vital_status = as.numeric(temp$Vital_status)
cr2 = coxph(Surv(OS, Vital_status) ~ NDST1+GALNT3+PDPK1+WDHD1+ZC2HC1A+XPNPEP1+SWAP70+MFHAS1+HILPDA+FAM117A, data = temp)
summary(cr2)
pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/forrest2.pdf", width = 12, height = 12) 
ggforest(cr2)
dev.off()
rm(temp)

#coxregression on ss3
temp = s3df
temp$Event = as.numeric(temp$Event)
cr3 = coxph(Surv(OS, Event) ~ GALNT3+XPNPEP1+SWAP70+PSIP1+MRPL9+NUS1P3+RCAN1+ITGB1+IL23A+PIAS1+HILPDA+AFF4+
                              LRFN4+TESMIN+DAB1+FAM117A+H2AW+STC2, data = temp)
summary(cr3)
pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/forrest3.pdf", width = 12, height = 12) 
ggforest(cr3)
dev.off()
rm(temp)

#coxregression on ss4
temp = s4df
temp$Event = as.numeric(temp$Event)
temp$Event[temp$Event == 3] = 2
cr4 = coxph(Surv(OS, Event) ~ TRIP12+BRD4+CHMP2B+NEBL+STEAP1+SLK+SOD2+VPS51+
            ERO1A+PTER+IFIH1+C1GALT1C1+STYK1+PPARD+GALNT3+UBA3+CCDC93+MIIP, data = temp)
summary(cr4)
pdf(file = "/Users/mdsayeefalam/Desktop/finalgraphs/forrest4.pdf", width = 12, height = 12) 
ggforest(cr4)
dev.off()
rm(temp)

t1 =as.data.frame(tidy(cr1))
t2 =as.data.frame(tidy(cr2))
t3 =as.data.frame(tidy(cr3))
t4 =as.data.frame(tidy(cr4))

t1 = subset(t1,t1$p.value<= 0.05)
t2 = subset(t2,t2$p.value<= 0.05)
t3 = subset(t3,t3$p.value<= 0.05)
t4 = subset(t4,t4$p.value<= 0.05)

t1 = t1[,-c(2:5)]
t2 = t2[,-c(2:5)]
t3 = t3[,-c(2:5)]

allg = as.data.frame(qpcR:::cbind.na(t1,t2,t3))
uni = stack(allg)
uni = uni[,-2]
uni = unique(uni)
uni = na.omit(uni)

genecsv = as.data.frame(qpcR:::cbind.na(uni,ss1gene_index,t1,ss2gene_index,t2,ss3gene_index,t3))
write.csv(genecsv,"/Users/mdsayeefalam/Desktop/finalgraphs/genecsv.csv")
geneindex = c(6642,10219,15263,19975,18855,19287,20728,2135,2924,4051,4255,4835,7951,9847,10219,10999,14581,14628,17229)

parsigenedf = as.data.frame(gmatrix[,geneindex])
colnames(parsigenedf) = c("GAPDHS","PSIP1","GM2A","DAB1","LRFN4","TRIM45","MLNR","NDST1","GALNT3","WDHD1","ZC2HC1A","XPNPEP1","HILPDA","SWAP70","RCAN1","ITGB1","IL23A","PIAS1","H2AW")
parsi_scaled = scale(parsigenedf,center = T,scale = T)

#heatmap
heatmap(as.matrix(parsi_scaled))

#correlation network map
install.packages("qgraph",dependencies = T)
install.packages("psych",dependencies = T)
library(qgraph)
library(psych)
CorMat1 = cor_auto(parsi_scaled)
Graph_pcor = qgraph(CorMat1,graph = "pcor",layout = "spring")

#brierplot
library(riskRegression)
brierdata = cbind(clinicaldata_clean,parsigenedf)
brierdata$Vital_Status = as.numeric(brierdata$Vital_Status)
brierdata$Relapse = as.numeric(brierdata$Relapse)
brierdata = cbind(brierdata,gs1 = rowSums(parsigenedf[1:19]))
f1 <- coxph(Surv(OS,Vital_Status)~gs1,data = brierdata,x=TRUE)
f2 <- coxph(Surv(PFS,Relapse)~gs1,data = brierdata,x=TRUE)
xscore <- Score(list("Gene Signature 1"=f1),formula=Surv(OS,Vital_Status)~1,data=brierdata,times=0:110,metrics="brier",null.model = F)
plotBrier(xscore)
mean(xscore$Brier$score$Brier)
mean(xscore$Brier$score$lower)
mean(xscore$Brier$score$upper)

