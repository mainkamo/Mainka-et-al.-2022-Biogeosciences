#                                                                         #
# Statistical analysis conducted for the                                  #
# Biogeosciences publication                                              #
# Mainka et al. (2022):                                                   #
# Soil geochemistry as a driver of soil organic matter composition        #
# Biogeosciences, 19, 1-15.                                               #
# Author of script:                                                       #
# Moritz Mainka                                                           #
# Date: 17.03.2022                                                        #
###########################################################################
rm(list=ls())

#load dataset and compute ratios ####
setwd("C:/Users/Moritz/Desktop/BachelorStudium/B. Sc. Thesis/Drafting/20200911_Simplified approach")
df <- read.table("20201027_data_exceptDeepSubsoil.txt", header = T, sep="\t")
df <- na.omit(df)

# Linear regression for C/N ratio ####
#Variable importance

# Define Formula
#without SOM variables
form_all2 <- as.formula(C.N~Al+Fe+Mg+Si+
                          Fe.total.DCB+Al.total.DCB+Ti.Zr+Fe.Si+Clay+Silt+Sand+CEC+
                          BSat+pH+TRB+BS.CEC+BS.pH+CEC.Clay)
library(caret)
# Create Resample (80 % of input data and 20 Resamples)
dp_all <- createDataPartition(df$C.N, p=.80, times=20) 
# Define Monte-Carlo cross-validation and model selection 
# LGOCV = Leave-one-group-out-cross-validation
tr_all <- trainControl(method = "LGOCV",p=0.8, number=100,index=dp_all,
                       preProcOptions=list(center=T,scale=T),
                       selectionFunction = "oneSE")
# Define predictors
#without SOM variables
impVars <- df[,c("Al","Fe","Mg","Si",
                 "Fe.total.DCB","Al.total.DCB","Ti.Zr","Fe.Si","Clay",
                 "Silt","Sand","CEC","BSat","pH","TRB","BS.CEC","BS.pH",
                 "CEC.Clay")]

# Variable Inflation Factor
library(car)
all_vifs_test <- try(vif(lm(form_all2, data=df)), silent=TRUE)
if (class(all_vifs_test) == "try-error"){
  lm_alias <- alias(lm(form_all2, data=df))
  broken <- data.frame(lm_alias$Complete)
  broken_var <- row.names(broken)
  nam_var <- names(impVars) %in% broken_var
  impVars <- impVars[!nam_var]
  form_all_new <- as.formula(paste("C.N ~ ", 
                                   paste(names(impVars), collapse=" + "),
                                   sep=""))
  all_vifs <- vif(lm(form_all_new, data=df))
} else {
  all_vifs <- all_vifs_test
}

if(any(all_vifs>5)){
  all_vifs <- as.data.frame(all_vifs)  
  while((nrow(all_vifs) > 2)& (max(all_vifs[, 1]) > 5)  &  
        (class(all_vifs) != "try-error")) {
    remove_var <- rownames(all_vifs)[which(all_vifs[, 1] == max(all_vifs[, 1]))]  
    impVars <- impVars[!names(impVars) %in% remove_var]  
    fullForm <- paste ("C.N ~ ", paste (names(impVars), collapse=" + "), sep="")  
    fullMod <- lm(as.formula(fullForm), data=df)  
    all_vifs <- try(as.data.frame(vif(fullMod)), silent=TRUE) 
  }
  vif_filtered_variables <- names(fullMod$model)[!names(fullMod$model) %in% 
                                                   "C.N"]
} else {
  all_vifs <- as.data.frame(all_vifs)
  vif_filtered_variables <- rownames(all_vifs)[!names(all_vifs) %in% "C.N"]
}

# Filtered Formula
form <- as.formula(paste("C.N~",paste(vif_filtered_variables, collapse='+')))

# Define Resampling method
dp <- createDataPartition(df$C.N, p=.80, times=20)
# Define cross-validation method
tr <- trainControl(method = "LGOCV",p=0.8, number=100,index=dp,
                   preProcOptions=list(center=T,scale=T),
                   selectionFunction = "oneSE")
linregCN <- train(form, data= df, method="lm",metric="RMSE", 
                    trControl=tr)
res_lmCN <- linregCN$results[,-1]
coef_lmCN <- as.data.frame(t(varImp(linregCN$finalModel)))
endtab_lmCN <- merge(res_lmCN, coef_lmCN)
colnames(endtab_lmCN) <- paste("LM", colnames(endtab_lmCN), sep="_")

# Linear regression for d15N value ####

#Variable importance

# Define Formula
#without SOM variables
form_all2 <- as.formula(d15N~Al+Fe+Mg+Si+
                          Fe.total.DCB+Al.total.DCB+Ti.Zr+Fe.Si+Clay+Silt+Sand+CEC+
                          BSat+pH+TRB+BS.CEC+BS.pH+CEC.Clay)
library(caret)
# Create Resample (80 % of input data and 20 Resamples)
dp_all <- createDataPartition(df$d15N, p=.80, times=20) 
# Define Monte-Carlo cross-validation and model selection 
# LGOCV = Leave-one-group-out-cross-validation
tr_all <- trainControl(method = "LGOCV",p=0.8, number=100,index=dp_all,
                       preProcOptions=list(center=T,scale=T),
                       selectionFunction = "oneSE")
# Define predictors
#without SOM variables
impVars <- df[,c("Al","Fe","Mg","Si",
                 "Fe.total.DCB","Al.total.DCB","Ti.Zr","Fe.Si","Clay",
                 "Silt","Sand","CEC","BSat","pH","TRB","BS.CEC","BS.pH",
                 "CEC.Clay")]

# Variable Inflation Factor
library(car)
all_vifs_test <- try(vif(lm(form_all2, data=df)), silent=TRUE)
if (class(all_vifs_test) == "try-error"){
  lm_alias <- alias(lm(form_all2, data=df))
  broken <- data.frame(lm_alias$Complete)
  broken_var <- row.names(broken)
  nam_var <- names(impVars) %in% broken_var
  impVars <- impVars[!nam_var]
  form_all_new <- as.formula(paste("d15N ~ ", 
                                   paste(names(impVars), collapse=" + "),
                                   sep=""))
  all_vifs <- vif(lm(form_all_new, data=df))
} else {
  all_vifs <- all_vifs_test
}

if(any(all_vifs>5)){
  all_vifs <- as.data.frame(all_vifs)  
  while((nrow(all_vifs) > 2)& (max(all_vifs[, 1]) > 5)  &  
        (class(all_vifs) != "try-error")) {
    remove_var <- rownames(all_vifs)[which(all_vifs[, 1] == max(all_vifs[, 1]))]  
    impVars <- impVars[!names(impVars) %in% remove_var]  
    fullForm <- paste ("d15N ~ ", paste (names(impVars), collapse=" + "), sep="")  
    fullMod <- lm(as.formula(fullForm), data=df)  
    all_vifs <- try(as.data.frame(vif(fullMod)), silent=TRUE) 
  }
  vif_filtered_variables <- names(fullMod$model)[!names(fullMod$model) %in% 
                                                   "d15N"]
} else {
  all_vifs <- as.data.frame(all_vifs)
  vif_filtered_variables <- rownames(all_vifs)[!names(all_vifs) %in% "d15N"]
}

# Filtered Formula
form <- as.formula(paste("d15N~",paste(vif_filtered_variables, collapse='+')))

# Define Resampling method
dp <- createDataPartition(df$d15N, p=.80, times=20)
# Define cross-validation method
tr <- trainControl(method = "LGOCV",p=0.8, number=100,index=dp,
                   preProcOptions=list(center=T,scale=T),
                   selectionFunction = "oneSE")
linregd15N <- train(form, data= df, method="lm",metric="RMSE", 
                  trControl=tr)
res_lmd15N <- linregd15N$results[,-1]
coef_lmd15N <- as.data.frame(t(varImp(linregd15N$finalModel)))
endtab_lmd15N <- merge(res_lmd15N, coef_lmd15N)
colnames(endtab_lmd15N) <- paste("LM", colnames(endtab_lmd15N), sep="_")

# Linear regression for d13C value ####
#Variable importance
#without SOM variables
form_all2 <- as.formula(d13C~Al+Fe+Mg+Si+
                          Fe.total.DCB+Al.total.DCB+Ti.Zr+Fe.Si+Clay+Silt+Sand+CEC+
                          BSat+pH+TRB+BS.CEC+BS.pH+CEC.Clay)
library(caret)
# Create Resample (80 % of input data and 20 Resamples)
dp_all <- createDataPartition(df$d13C, p=.80, times=20) 
# Define Monte-Carlo cross-validation and model selection 
# LGOCV = Leave-one-group-out-cross-validation
tr_all <- trainControl(method = "LGOCV",p=0.8, number=100,index=dp_all,
                       preProcOptions=list(center=T,scale=T),
                       selectionFunction = "oneSE")
# Define predictors
#without SOM variables
impVars <- df[,c("Al","Fe","Mg","Si",
                 "Fe.total.DCB","Al.total.DCB","Ti.Zr","Fe.Si","Clay",
                 "Silt","Sand","CEC","BSat","pH","TRB","BS.CEC","BS.pH",
                 "CEC.Clay")]

# Variable Inflation Factor
library(car)
all_vifs_test <- try(vif(lm(form_all2, data=df)), silent=TRUE)
if (class(all_vifs_test) == "try-error"){
  lm_alias <- alias(lm(form_all2, data=df))
  broken <- data.frame(lm_alias$Complete)
  broken_var <- row.names(broken)
  nam_var <- names(impVars) %in% broken_var
  impVars <- impVars[!nam_var]
  form_all_new <- as.formula(paste("d13C ~ ", 
                                   paste(names(impVars), collapse=" + "),
                                   sep=""))
  all_vifs <- vif(lm(form_all_new, data=df))
} else {
  all_vifs <- all_vifs_test
}

if(any(all_vifs>5)){
  all_vifs <- as.data.frame(all_vifs)  
  while((nrow(all_vifs) > 2)& (max(all_vifs[, 1]) > 5)  &  
        (class(all_vifs) != "try-error")) {
    remove_var <- rownames(all_vifs)[which(all_vifs[, 1] == max(all_vifs[, 1]))]  
    impVars <- impVars[!names(impVars) %in% remove_var]  
    fullForm <- paste ("d13C ~ ", paste (names(impVars), collapse=" + "), sep="")  
    fullMod <- lm(as.formula(fullForm), data=df)  
    all_vifs <- try(as.data.frame(vif(fullMod)), silent=TRUE) 
  }
  vif_filtered_variables <- names(fullMod$model)[!names(fullMod$model) %in% 
                                                   "d13C"]
} else {
  all_vifs <- as.data.frame(all_vifs)
  vif_filtered_variables <- rownames(all_vifs)[!names(all_vifs) %in% "d13C"]
}

# Filtered Formula
form <- as.formula(paste("d13C~",paste(vif_filtered_variables, collapse='+')))

# Define Resampling method
dp <- createDataPartition(df$d13C, p=.80, times=20)
# Define cross-validation method
tr <- trainControl(method = "LGOCV",p=0.8, number=100,index=dp,
                   preProcOptions=list(center=T,scale=T),
                   selectionFunction = "oneSE")
linregd13C <- train(form, data= df, method="lm",metric="RMSE", 
                    trControl=tr)
res_lmd13C <- linregd13C$results[,-1]
coef_lmd13C <- as.data.frame(t(varImp(linregd13C$finalModel)))
endtab_lmd13C <- merge(res_lmd13C, coef_lmd13C)
colnames(endtab_lmd13C) <- paste("LM", colnames(endtab_lmd13C), sep="_")

# Linear regression for aliphatic C-H ####
#Variable importance

# Define Formula
#without SOM variables
form_all2 <- as.formula(aliphaticCH~Al+Fe+Mg+Si+
                          Fe.total.DCB+Al.total.DCB+Ti.Zr+Fe.Si+Clay+Silt+Sand+CEC+
                          BSat+pH+TRB+BS.CEC+BS.pH+CEC.Clay)
library(caret)
# Create Resample (80 % of input data and 20 Resamples)
dp_all <- createDataPartition(df$aliphaticCH, p=.80, times=20) 
# Define Monte-Carlo cross-validation and model selection 
# LGOCV = Leave-one-group-out-cross-validation
tr_all <- trainControl(method = "LGOCV",p=0.8, number=100,index=dp_all,
                       preProcOptions=list(center=T,scale=T),
                       selectionFunction = "oneSE")
# Define predictors
#without SOM variables
impVars <- df[,c("Al","Fe","Mg","Si",
                 "Fe.total.DCB","Al.total.DCB","Ti.Zr","Fe.Si","Clay",
                 "Silt","Sand","CEC","BSat","pH","TRB","BS.CEC","BS.pH",
                 "CEC.Clay")]

# Variable Inflation Factor
library(car)
all_vifs_test <- try(vif(lm(form_all2, data=df)), silent=TRUE)
if (class(all_vifs_test) == "try-error"){
  lm_alias <- alias(lm(form_all2, data=df))
  broken <- data.frame(lm_alias$Complete)
  broken_var <- row.names(broken)
  nam_var <- names(impVars) %in% broken_var
  impVars <- impVars[!nam_var]
  form_all_new <- as.formula(paste("aliphaticCH ~ ", 
                                   paste(names(impVars), collapse=" + "),
                                   sep=""))
  all_vifs <- vif(lm(form_all_new, data=df))
} else {
  all_vifs <- all_vifs_test
}

if(any(all_vifs>5)){
  all_vifs <- as.data.frame(all_vifs)  
  while((nrow(all_vifs) > 2)& (max(all_vifs[, 1]) > 5)  &  
        (class(all_vifs) != "try-error")) {
    remove_var <- rownames(all_vifs)[which(all_vifs[, 1] == max(all_vifs[, 1]))]  
    impVars <- impVars[!names(impVars) %in% remove_var]  
    fullForm <- paste ("aliphaticCH ~ ", paste (names(impVars), collapse=" + "), sep="")  
    fullMod <- lm(as.formula(fullForm), data=df)  
    all_vifs <- try(as.data.frame(vif(fullMod)), silent=TRUE) 
  }
  vif_filtered_variables <- names(fullMod$model)[!names(fullMod$model) %in% 
                                                   "aliphaticCH"]
} else {
  all_vifs <- as.data.frame(all_vifs)
  vif_filtered_variables <- rownames(all_vifs)[!names(all_vifs) %in% "aliphaticCH"]
}

# Filtered Formula
form <- as.formula(paste("aliphaticCH~",paste(vif_filtered_variables, collapse='+')))

# Define Resampling method
dp <- createDataPartition(df$aliphaticCH, p=.80, times=20)
# Define cross-validation method
tr <- trainControl(method = "LGOCV",p=0.8, number=100,index=dp,
                   preProcOptions=list(center=T,scale=T),
                   selectionFunction = "oneSE")
linregCH <- train(form, data= df, method="lm",metric="RMSE", 
                    trControl=tr)
res_lmCH <- linregCH$results[,-1]
coef_lmCH <- as.data.frame(t(varImp(linregCH$finalModel)))
endtab_lmCH <- merge(res_lmCH, coef_lmCH)
colnames(endtab_lmCH) <- paste("LM", colnames(endtab_lmCH), sep="_")


# Linear regression for aromatic C=C ####
#Variable importance
#without SOM variables
form_all2 <- as.formula(aromaticCC~Al+Fe+Mg+Si+
                          Fe.total.DCB+Al.total.DCB+Ti.Zr+Fe.Si+Clay+Silt+Sand+CEC+
                          BSat+pH+TRB+BS.CEC+BS.pH+CEC.Clay)
library(caret)
# Create Resample (80 % of input data and 20 Resamples)
dp_all <- createDataPartition(df$aromaticCC, p=.80, times=20) 
# Define Monte-Carlo cross-validation and model selection 
# LGOCV = Leave-one-group-out-cross-validation
tr_all <- trainControl(method = "LGOCV",p=0.8, number=100,index=dp_all,
                       preProcOptions=list(center=T,scale=T),
                       selectionFunction = "oneSE")
# Define predictors
#without SOM variables
impVars <- df[,c("Al","Fe","Mg","Si",
                 "Fe.total.DCB","Al.total.DCB","Ti.Zr","Fe.Si","Clay",
                 "Silt","Sand","CEC","BSat","pH","TRB","BS.CEC","BS.pH",
                 "CEC.Clay")]

# Variable Inflation Factor
library(car)
all_vifs_test <- try(vif(lm(form_all2, data=df)), silent=TRUE)
if (class(all_vifs_test) == "try-error"){
  lm_alias <- alias(lm(form_all2, data=df))
  broken <- data.frame(lm_alias$Complete)
  broken_var <- row.names(broken)
  nam_var <- names(impVars) %in% broken_var
  impVars <- impVars[!nam_var]
  form_all_new <- as.formula(paste("aromaticCC ~ ", 
                                   paste(names(impVars), collapse=" + "),
                                   sep=""))
  all_vifs <- vif(lm(form_all_new, data=df))
} else {
  all_vifs <- all_vifs_test
}

if(any(all_vifs>5)){
  all_vifs <- as.data.frame(all_vifs)  
  while((nrow(all_vifs) > 2)& (max(all_vifs[, 1]) > 5)  &  
        (class(all_vifs) != "try-error")) {
    remove_var <- rownames(all_vifs)[which(all_vifs[, 1] == max(all_vifs[, 1]))]  
    impVars <- impVars[!names(impVars) %in% remove_var]  
    fullForm <- paste ("aromaticCC ~ ", paste (names(impVars), collapse=" + "), sep="")  
    fullMod <- lm(as.formula(fullForm), data=df)  
    all_vifs <- try(as.data.frame(vif(fullMod)), silent=TRUE) 
  }
  vif_filtered_variables <- names(fullMod$model)[!names(fullMod$model) %in% 
                                                   "aromaticCC"]
} else {
  all_vifs <- as.data.frame(all_vifs)
  vif_filtered_variables <- rownames(all_vifs)[!names(all_vifs) %in% "aromaticCC"]
}

# Filtered Formula
form <- as.formula(paste("aromaticCC~",paste(vif_filtered_variables, collapse='+')))

# Linear regression for carboxylate C=O, amide/quetone/quinones ####
#Variable importance

#without SOM variables
form_all2 <- as.formula(amiquiket~Al+Fe+Mg+Si+
                          Fe.total.DCB+Al.total.DCB+Ti.Zr+Fe.Si+Clay+Silt+Sand+CEC+
                          BSat+pH+TRB+BS.CEC+BS.pH+CEC.Clay)
library(caret)
# Create Resample (80 % of input data and 20 Resamples)
dp_all <- createDataPartition(df$amiquiket, p=.80, times=20) 
# Define Monte-Carlo cross-validation and model selection 
# LGOCV = Leave-one-group-out-cross-validation
tr_all <- trainControl(method = "LGOCV",p=0.8, number=100,index=dp_all,
                       preProcOptions=list(center=T,scale=T),
                       selectionFunction = "oneSE")
# Define predictors
#without SOM variables
impVars <- df[,c("Al","Fe","Mg","Si",
                 "Fe.total.DCB","Al.total.DCB","Ti.Zr","Fe.Si","Clay",
                 "Silt","Sand","CEC","BSat","pH","TRB","BS.CEC","BS.pH",
                 "CEC.Clay")]

# Variable Inflation Factor
library(car)
all_vifs_test <- try(vif(lm(form_all2, data=df)), silent=TRUE)
if (class(all_vifs_test) == "try-error"){
  lm_alias <- alias(lm(form_all2, data=df))
  broken <- data.frame(lm_alias$Complete)
  broken_var <- row.names(broken)
  nam_var <- names(impVars) %in% broken_var
  impVars <- impVars[!nam_var]
  form_all_new <- as.formula(paste("amiquiket ~ ", 
                                   paste(names(impVars), collapse=" + "),
                                   sep=""))
  all_vifs <- vif(lm(form_all_new, data=df))
} else {
  all_vifs <- all_vifs_test
}

if(any(all_vifs>5)){
  all_vifs <- as.data.frame(all_vifs)  
  while((nrow(all_vifs) > 2)& (max(all_vifs[, 1]) > 5)  &  
        (class(all_vifs) != "try-error")) {
    remove_var <- rownames(all_vifs)[which(all_vifs[, 1] == max(all_vifs[, 1]))]  
    impVars <- impVars[!names(impVars) %in% remove_var]  
    fullForm <- paste ("amiquiket ~ ", paste (names(impVars), collapse=" + "), sep="")  
    fullMod <- lm(as.formula(fullForm), data=df)  
    all_vifs <- try(as.data.frame(vif(fullMod)), silent=TRUE) 
  }
  vif_filtered_variables <- names(fullMod$model)[!names(fullMod$model) %in% 
                                                   "amiquiket"]
} else {
  all_vifs <- as.data.frame(all_vifs)
  vif_filtered_variables <- rownames(all_vifs)[!names(all_vifs) %in% "amiquiket"]
}

# Filtered Formula
form <- as.formula(paste("amiquiket~",paste(vif_filtered_variables, collapse='+')))

# Define Resampling method
dp <- createDataPartition(df$amiquiket, p=.80, times=20)
# Define cross-validation method
tr <- trainControl(method = "LGOCV",p=0.8, number=100,index=dp,
                   preProcOptions=list(center=T,scale=T),
                   selectionFunction = "oneSE")
linregCO <- train(form, data= df, method="lm",metric="RMSE", 
                  trControl=tr)
res_lmCO <- linregCO$results[,-1]
coef_lmCO <- as.data.frame(t(varImp(linregCO$finalModel)))
endtab_lmCO <- merge(res_lmCO, coef_lmCO)
colnames(endtab_lmCO) <- paste("LM", colnames(endtab_lmCO), sep="_")

#export tables from linregs
library(plyr)
linregs <- rbind.fill(endtab_lmCN, endtab_lmd15N, endtab_lmd13C, endtab_lmCH,endtab_lmCC, endtab_lmCO)
`rownames<-`(linregs, c("C.N","d15N","d13C","C-H/C=O","C=C/C=O"))
library(openxlsx)
write.xlsx(linregs,"lmvar_output.xlsx")

# calculate rRMSE and rMAE ####
2.0/mean(df$C.N)*100
0.6/mean(df$d13C)*100
1.2/mean(df$d15N)*100
2.9/mean(df$aliphaticCH)*100
3.2/mean(df$aromaticCC)*100
20.3/mean(df$amiquiket)*100

#rMAE
1.8/mean(df$C.N)*100
0.6/mean(df$d13C)*100
1.0/mean(df$d15N)*100
2.5/mean(df$aliphaticCH)*100
2.6/mean(df$aromaticCC)*100
18.1/mean(df$amiquiket)*100


# One-way ANOVA #####
library(agricolae)
#dataset preparation
df$Soil.age <- as.character(df$Soil.age)
df$Soil.age <- factor(df$Soil.age)
df$Mean.depth <- as.character(df$Mean.depth)
d$Mean.depth <- factor(df$Mean.depth)

#anova for d15N
anova.15N <- aov(d15N [Mean.depth == 20] ~ Soil.age[Mean.depth == 20], data = df)
anova.15N <- aov(d15N ~ Soil.age + Mean.depth, data = df)
summary(anova.15N)
require("dplyr")
group_by(df, Soil.age) %>%
  summarise(
    count = n(),
    mean = mean(d15N, na.rm = TRUE),
    sd = sd(d15N, na.rm = TRUE))
#Tukey HSD
posthoc.15N <- TukeyHSD(anova.15N)
posthoc.15N
x15n <- HSD.test(anova.15N, "Soil.age", group=TRUE)
x15n
require("car")
leveneTest(df.org$d15N[df.org$Mean.depth == 20], df.org$Soil.age[df.org$Mean.depth == 20], data = df.org)
aov_residuals <- residuals(object = anova.15N)
shapiro.test(x = aov_residuals )


#one-way anova for d13C
anova.13C <- aov(d13C[Mean.depth == 20] ~ Soil.age[Mean.depth == 20], data = df)
anova.13C <- aov(d13C ~ Soil.age+Mean.depth, data = df)
group_by(df.org, Soil.age) %>%
  summarise(
    count = n(),
    mean = mean(d13C, na.rm = TRUE),
    sd = sd(d13C, na.rm = TRUE))
#Tukey HSD
posthoc13C <- TukeyHSD(anova.13C)
posthoc13C
x13c <- HSD.test(anova.13C, "Mean.depth", group=TRUE)
x13c
leveneTest(df.org$d13C[df.org$Mean.depth == 5], df.org$Soil.age[df.org$Mean.depth == 5], data = df.org)
aov_residuals <- residuals(object = anova.13C)
shapiro.test(x = aov_residuals )

# C stocks natgeo 
C <- read.table("Cstocks_Doetterl.txt", header = T, sep = "\t")
str(C)
#C$Soil.age <- factor(C$Soil.age)
#levels(C$Soil.age) <- c("0.1","3","19","295","3000")
C$Mean.depth <- factor(C$Mean.depth)
levels(C$Mean.depth)
aaov.C <- aov(C_perc ~ Soil.age + Mean.depth, data = C)
summary(aaov.C)
HSD.c <- HSD.test(aaov.C, "Soil.age")
HSD.c

#ANOVA for C/N
anova.CN <- aov(C.N[Mean.depth == 20] ~ Soil.age[Mean.depth == 20], data = df)
anova.CN <- aov(C.N ~ Soil.age + Mean.depth, data = df)
group_by(df.org, Soil.age) %>%
  summarise(
    count = n(),
    mean = mean(C.N, na.rm = TRUE),
    sd = sd(C.N, na.rm = TRUE))
#Tukey HSD
posthocCN <- TukeyHSD(anova.CN)
plot(posthocCN)
posthocCN
xCN <- HSD.test(anova.CN, "Mean.depth", group=TRUE)
xCN
leveneTest(df.org$C.N[df.org$Mean.depth == 5], df.org$Soil.age[df.org$Mean.depth == 5], data = df.org)
aov_residuals <- residuals(object = anova.CN)
shapiro.test(x = aov_residuals)

#one-way anova for C-H
anova.CH <- aov(aliphaticCH[Mean.depth == 20] ~ Soil.age[Mean.depth == 20], data = df)
anova.CH <- aov(aliphaticCH ~ Soil.age + Mean.depth, data = df)
group_by(df.org, Soil.age) %>%
  summarise(
    count = n(),
    mean = mean(alipcarbox, na.rm = TRUE),
    sd = sd(alipcarbox, na.rm = TRUE))
#Tukey HSD
posthocCH <- TukeyHSD(anova.CH)
posthocCH
xCH <- HSD.test(anova.CH, "Soil.age", group=TRUE)
xCH

leveneTest(df.org$alipcarbox[df.org.Mean.depth == 20], df.org$Soil.age[df.org.Mean.depth == 20], data = df.org)
aov_residuals <- residuals(object = anova.CH)
shapiro.test(x = aov_residuals)


#one-way anova for C=C
anova.CC <- aov(aromaticCC[Mean.depth == 20] ~ Soil.age[Mean.depth == 20], data = df)
anova.CC <- aov(aromaticCC ~ Soil.age + Mean.depth, data = df)
#Tukey HSD
posthocCC <- TukeyHSD(anova.CC)
posthocCC
xCC <- HSD.test(anova.CC, "Soil.age", group=TRUE)
xCC

leveneTest(aromcarbox, df.org.Soil.age[df.org.Mean.depth == 5], data = dat)
aov_residuals <- residuals(object = anova.CC)
shapiro.test(x = aov_residuals)

#one-way anova for MOM
anova.CO <- aov(amiquiket[Mean.depth == 5] ~ Soil.age[Mean.depth == 5], data = df)
anova.CO <- aov(amiquiket ~ Soil.age + Mean.depth, data = df)
#Tukey HSD
posthocCO <- TukeyHSD(anova.CO)
posthocCO
xCO <- HSD.test(anova.CO, "Soil.age", group=TRUE)