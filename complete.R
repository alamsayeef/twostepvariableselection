#import data
genedata_orig = read.csv(file.choose())
clinicaldata_orig = read.csv(file.choose())

#display data
head(genedata_orig)
head(clinicaldata_orig)

#cleaning packages
install.packages("tidyverse", repos = "http://cran.us.r-project.org")

library(tidyverse)

#gene data cleaning
genedata_clean = as.data.frame(t(genedata_orig), header = T)
names(genedata_clean) = genedata_clean %>% slice(1) %>% unlist()
genedata_clean = genedata_clean %>% slice(-1)
genedata_clean = mutate_all(genedata_clean, function(x) as.numeric(as.character(x)))

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

#after cleaning clinical data
summary(clinicaldata_clean)

#merging both data sets
clinicaldata_clean$ID = seq(1:443)
genedata_clean$ID = seq(1:443)
mergeddata = merge(clinicaldata_clean,genedata_clean, by = "ID")

#make subsets of data as per the event
cancerdata = subset(mergeddata,mergeddata$Event != "Competing Risk")
cancerdata$Event = as.numeric(c("Alive" = "0", "Cancer" = "1")[cancerdata$Event])

cmprskdata = subset(mergeddata,mergeddata$Event != "Cancer")
cmprskdata$Event = as.numeric(c("Alive" = "0", "Competing Risk" = "1")[cmprskdata$Event])

#define survival objects
library(survival)
so1 = Surv(mergeddata$PFS,mergeddata$Relapse)
so2 = Surv(mergeddata$OS,mergeddata$Vital_Status)
so3 = Surv(cancerdata$OS,cancerdata$Event)
so4 = Surv(cmprskdata$OS,cmprskdata$Event)

#genematrix
gmatrix = data.matrix(mergeddata[,16:22230])
gmatrixe1 = data.matrix(cancerdata[,16:22230])
gmatrixe2 = data.matrix(cmprskdata[,16:22230])

#Variable selections under varying survival scenarios
install.packages("glmnet", repos = "http://cran.us.r-project.org")
library(glmnet)

#setseed function for exact replication later
set.seed(98765)

#lasso ss1
lr1 = cv.glmnet(gmatrix,so1,alpha=1,family = "cox",nfolds = 10)
plot(lr1)
lr1$lambda.min
log(lr1$lambda.min)
lr1c = coef(lr1, s = lr1$lambda.min)
lr1lam = which(lr1c!=0)
lr1lam
lr1lamval = lr1c[lr1lam]
lr1lamval

#lasso ss2
lr2 = cv.glmnet(gmatrix,so2,alpha=1,family = "cox",nfolds = 10)
plot(lr2)
lr2$lambda.min
log(lr2$lambda.min)
lr2c = coef(lr2, s = lr2$lambda.min)
lr2lam = which(lr2c!=0)
lr2lam
lr2lamval = lr2c[lr2lam]
lr2lamval

#lasso ss3
lr3 = cv.glmnet(gmatrixe1,so3,alpha=1,family = "cox",nfolds = 10)
plot(lr3)
lr3$lambda.min
log(lr3$lambda.min)
lr3c = coef(lr3, s = lr3$lambda.min)
lr3lam = which(lr3c!=0)
lr3lam
lr3lamval = lr3c[lr3lam]
lr3lamval

#lasso ss4
lr4 = cv.glmnet(gmatrixe2,so4,alpha=1,family = "cox",nfolds = 10)
plot(lr4)
lr4$lambda.min
log(lr4$lambda.min)
lr4c = coef(lr4, s = lr4$lambda.min)
lr4lam = which(lr4c!=0)
lr4lam
lr4lamval = lr4c[lr4lam]
lr4lamval

#enet ss1
er1 = cv.glmnet(gmatrix,so1,alpha=0.5,family = "cox", nfolds = 10)
plot(er1)
er1$lambda.min
log(er1$lambda.min)
er1c = coef(er1, s = er1$lambda.min)
er1lam = which(er1c!=0)
er1lam
er1lamval = er1c[er1lam]
er1lamval

#enet ss2
er2 = cv.glmnet(gmatrix,so2,alpha=0.5,family = "cox", nfolds = 10)
plot(er2)
er2$lambda.min
log(er2$lambda.min)
er2c = coef(er2, s = er2$lambda.min)
er2lam = which(er2c!=0)
er2lam
er2lamval = er2c[er2lam]
er2lamval

#enet ss3
er3 = cv.glmnet(gmatrixe1,so3,alpha=0.5,family = "cox", nfolds = 10)
plot(er3)
er3$lambda.min
log(er3$lambda.min)
er3c = coef(er3, s = er3$lambda.min)
er3lam = which(er3c!=0)
er3lam
er3lamval = er3c[er3lam]
er3lamval

#enet ss4
er4 = cv.glmnet(gmatrixe2,so4,alpha=0.5,family = "cox", nfolds = 10)
plot(er4)
er4$lambda.min
log(er4$lambda.min)
er4c = coef(er4, s = er4$lambda.min)
er4lam = which(er4c!=0)
er4lam
er4lamval = er4c[er4lam]
er4lamval

#extracting the selected genes and their corresponding values
install.packages("qpcR", repos = "http://cran.us.r-project.org", dependencies = T)
library(qpcR)

selectedgenes = as.data.frame(qpcR:::cbind.na(er1lam,er2lam,er3lam,er4lam,lr1lam,lr2lam,lr4lam))
allgene = stack(selectedgenes)
allgene = allgene[,-2]
allgene = as.data.frame(allgene)
uniquegene = unique(allgene)
uniquegene = uniquegene %>% rename(uniquegene = allgene)
uniquegene = sort(uniquegene$uniquegene)
selectedgenes = as.data.frame(qpcR:::cbind.na(er1lam,er2lam,er3lam,er4lam,lr1lam,lr2lam,lr4lam,allgene,uniquegene))

#saving the selected variants onto a csv file
write.csv(selectedgenes, file = "/Users/mdsayeefalam/Documents/Working Papers/TMH/geneset.csv")

#keeping only the required covariates
finalgenedata = genedata_clean %>% select(uniquegene)

finaldata = as.data.frame(c(clinicaldata_clean,finalgenedata))
finaldata = finaldata[,-c(1,16)]
write.csv(finaldata, file = "/Users/mdsayeefalam/Documents/Working Papers/TMH/finaldata.csv")

#applying cox proportional hazards regression to determine
#significant genes in all four survival scenario
genes = finaldata[,15:66]
covariates = colnames(genes)
finaldata$Relapse = as.numeric(c("No" = "0", "Yes" = "1")[finaldata$Relapse])
finaldata$Vital_Status = as.numeric(c("Alive" = "0", "Dead" = "1")[finaldata$Vital_Status])

fcandata = subset(finaldata,finaldata$Event != "Competing Risk")
fcandata$Event = as.numeric(c("Alive" = "0", "Cancer" = "1")[fcandata$Event])
fcmpdata = subset(finaldata,finaldata$Event != "Cancer")
fcmpdata$Event = as.numeric(c("Alive" = "0", "Competing Risk" = "1")[fcmpdata$Event])

s1 = Surv(finaldata$PFS,finaldata$Relapse)
s2 = Surv(finaldata$OS,finaldata$Vital_Status)
s3 = Surv(fcandata$OS,fcandata$Event)
s4 = Surv(fcmpdata$OS,fcmpdata$Event)

#surivival scenario 1
univ_formulas1 = sapply(covariates, function(x) as.formula(paste('s1~', x)))
univ_models1 = lapply( univ_formulas1, function(x){coxph(x, data = finaldata)})
# Extract data 
univ_results1 = lapply(univ_models1,
                       function(x){ 
                         x = summary(x)
                         p.value = signif(x$wald["pvalue"], digits=2)
                         wald.test = signif(x$wald["test"], digits=2)
                         beta = signif(x$coef[1], digits=2);#coeficient beta
                         HR = signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
                         HR = paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
                         res1 = c(beta, HR, wald.test, p.value)
                         names(res1) = c("beta", "HR (95% CI for HR)", "wald.test","p.value")
                         return(res1)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res1 = t(as.data.frame(univ_results1, check.names = FALSE))
as.data.frame(res1)

#surivival scenario 2
univ_formulas2 = sapply(covariates, function(x) as.formula(paste('s2~', x)))
univ_models2 = lapply( univ_formulas2, function(x){coxph(x, data = finaldata)})
# Extract data 
univ_results2 = lapply(univ_models2,
                       function(x){ 
                         x = summary(x)
                         p.value = signif(x$wald["pvalue"], digits=2)
                         wald.test = signif(x$wald["test"], digits=2)
                         beta = signif(x$coef[1], digits=2);#coeficient beta
                         HR = signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
                         HR = paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
                         res2 = c(beta, HR, wald.test, p.value)
                         names(res2) = c("beta", "HR (95% CI for HR)", "wald.test","p.value")
                         return(res2)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res2 = t(as.data.frame(univ_results2, check.names = FALSE))
as.data.frame(res2)

#surivival scenario 3
univ_formulas3 = sapply(covariates, function(x) as.formula(paste('s3~', x)))
univ_models3 = lapply( univ_formulas3, function(x){coxph(x, data = fcandata)})
# Extract data 
univ_results3 = lapply(univ_models3,
                       function(x){ 
                         x = summary(x)
                         p.value = signif(x$wald["pvalue"], digits=2)
                         wald.test = signif(x$wald["test"], digits=2)
                         beta = signif(x$coef[1], digits=2);#coeficient beta
                         HR = signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
                         HR = paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
                         res3 = c(beta, HR, wald.test, p.value)
                         names(res3) = c("beta", "HR (95% CI for HR)", "wald.test","p.value")
                         return(res3)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res3 = t(as.data.frame(univ_results3, check.names = FALSE))
as.data.frame(res3)

#surivival scenario 4
univ_formulas4 = sapply(covariates, function(x) as.formula(paste('s4~', x)))
univ_models4 = lapply( univ_formulas4, function(x){coxph(x, data = fcmpdata)})
# Extract data 
univ_results4 = lapply(univ_models4,
                       function(x){ 
                         x = summary(x)
                         p.value = signif(x$wald["pvalue"], digits=2)
                         wald.test = signif(x$wald["test"], digits=2)
                         beta = signif(x$coef[1], digits=2);#coeficient beta
                         HR = signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
                         HR = paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
                         res4 = c(beta, HR, wald.test, p.value)
                         names(res4) = c("beta", "HR (95% CI for HR)", "wald.test","p.value")
                         return(res4)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res4 = t(as.data.frame(univ_results4, check.names = FALSE))
as.data.frame(res4)

coxresult = as.data.frame(cbind(res1,res2,res3,res4))
write.csv(coxresult,file = "/Users/mdsayeefalam/Documents/Working Papers/TMH/coxresult.csv")

library(dplyr)

d1 = as.data.frame(res1)
d1$p.value = as.numeric(d1$p.value)
sg1 = subset(d1, d1$p.value < 0.05)
sg1 = tibble::rownames_to_column(sg1, "gene")

d2 = as.data.frame(res2)
d2$p.value = as.numeric(d2$p.value)
sg2 = subset(d2, d2$p.value < 0.05)
sg2 = tibble::rownames_to_column(sg2, "gene")

d3 = as.data.frame(res3)
d3$p.value = as.numeric(d3$p.value)
sg3 = subset(d3, d3$p.value < 0.05)
sg3 = tibble::rownames_to_column(sg3, "gene")

d4 = as.data.frame(res4)
d4$p.value = as.numeric(d4$p.value)
sg4 = subset(d4, d4$p.value < 0.05)
sg4 = tibble::rownames_to_column(sg4, "gene")

t1 = intersect(sg1$gene,sg2$gene)
t2 = intersect(sg3$gene,sg4$gene)
t3 = intersect(t1,t2)

signigene = as.data.frame(t3)
write.csv(signigene,file = "/Users/mdsayeefalam/Documents/Working Papers/TMH/parsigene.csv")

library(stringr)
signigene$t3 = str_sub(signigene$t3,2)
finaldata = genedata_clean %>% select(signigene$t3)
finaldata = as.data.frame(cbind(clinicaldata_clean, finaldata))
                          
#name change for probeset genes
names(finaldata)[16] =  'WDHD1'
names(finaldata)[17] =  'MRPL9'
names(finaldata)[18] =  'XPNPEP1'
names(finaldata)[19] =  'GALNT3'
names(finaldata)[20] =  'FAM117A'
names(finaldata)[21] =  'NUS1'
names(finaldata)[22] =  'ZC2HC1A'
names(finaldata)[23] =  'STC2'
names(finaldata)[24] =  'HILPDA'
names(finaldata)[25] =  'PDPK1'
names(finaldata)[26] =  'RPE'
names(finaldata)[27] =  'CCDC93'
names(finaldata)[28] =  'PIAS1'
names(finaldata)[29] =  'NDST1'
##################

parsidf = finaldata %>% select(16:29)

fcandata = subset(finaldata,finaldata$Event != "Competing Risk")
fcandata$Event = as.numeric(c("Alive" = "0", "Cancer" = "1")[fcandata$Event])
fcmpdata = subset(finaldata,finaldata$Event != "Cancer")
fcmpdata$Event = as.numeric(c("Alive" = "0", "Competing Risk" = "1")[fcmpdata$Event])
canparsidf = fcandata %>% select(16:29)
cmpparsidf = fcmpdata %>% select(16:29)

finaldata = cbind(finaldata,gs1 = rowSums(finaldata[16:29]))
fcandata = cbind(fcandata,gs2 = rowSums(fcandata[16:29]))
fcmpdata = cbind(fcmpdata,gs3 = rowSums(fcmpdata[16:29]))

#saving the 
write.csv(parsidf, file = "/Users/mdsayeefalam/Documents/Working Papers/TMH/parsidf.csv")

#heatmap
matparsi = as.matrix(parsidf)
matparsi = scale(matparsi,center = T)
heatmap(matparsi)

#correlation network map
install.packages("qgraph",dependencies = T)
install.packages("psych",dependencies = T)
library(qgraph)
library(psych)
CorMat1 = cor_auto(parsidf)
Graph_pcor = qgraph(CorMat1,graph = "pcor",layout = "spring")

#brierscore and brierplot
library(survival)
library(prodlim)
library(riskRegression)
#f1 <- coxph(Surv(OS,Vital_Status)~gs1,data = finaldata,x=TRUE)
f2 <- coxph(Surv(OS,Event)~gs2,data = fcandata,x=TRUE)
#f3 <- coxph(Surv(OS,Event)~gs3,data = fcmpdata,x=TRUE)

xscore <- Score(list("Gene Signature"=f2),formula=Surv(OS,Event)~1,data=fcandata,times=0:200,metrics="brier",null.model = F)
plotBrier(xscore)

#adjusted  p values
signigene = signigene[order(signigene$p.value),]
genetestdata = signigene[!duplicated(signigene$gene),]

install.packages("stats")
library(stats)
pvalBH = p.adjust(genetestdata$p.value, method = "BH", n = length(genetestdata$p.value))
pvalfdr = p.adjust(genetestdata$p.value, method = "fdr", n = length(genetestdata$p.value))
pvalBonn = p.adjust(genetestdata$p.value, method = "bonferroni", n = length(genetestdata$p.value))
t = as.data.frame(cbind(pvalBH,pvalBonn,pvalfdr))

#plot adjusted pvalues
hist(t$pvalBH)
hist(t$pvalBonn)
hist(t$pvalfdr)

#FDR
install.packages("fuzzySim")
library(fuzzySim)
fdr.fit = FDR(data = mergeddata, sp.cols = 5, var.cols = 16:22230,family = "auto", correction = "fdr", q = 0.05, verbose = TRUE, simplif = FALSE)
hist(fdr.fit$exclude$p.value)
hist(fdr.fit$select$p.value)