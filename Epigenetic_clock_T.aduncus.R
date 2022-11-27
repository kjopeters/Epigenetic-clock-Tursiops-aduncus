#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
# May 2022
# Epigenetic clock Tursiops aduncus
#
# Katharina Peters & Livia Gerber
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
library(glmnet)
library(coefplot)
library(caret)

setwd("~/Documents/Work/UZH/Manuscripts/Peters et al. EpiClock/Analysis/EpiClock")

set.seed(1234567)

#........................................................
#
#      Data import
#........................................................

#Read in epigenetic data
epi <- readRDS("data/epi_shuff_365.rds")


#delete bottom 62 sites
epi <- epi[1:37492,]

#read in non-randomised file for CGid
epi2 <-read.csv("data/normalized_betasN93_N139_730confidenceNoDups.csv")

#keep CGids and add intercept to merge with non-zero coefficients further below
CGid <- data.frame(epi2$CGid)
CGid <- rbind("intercept",CGid)

#delete bottom 62 sites
CGid <- data.frame(CGid[1:37493,])

#Read in metainfo with known chronological ages
info <- readRDS("data/info_shuff_365.rds")

#........................................................
#
#      LOG LINEAR TRANSFORMATION OF AGE
#
##This is needed because individuals age faster until reaching maturity
# ASM = Age of maturity -> separate males and females and run separate clocks 
# based on the age of maturity of both sexes

## Functions are based on log-linear age for Clock 3 in universal DNA methylation age paper
# Transforms age
#........................................................
matAGEfem <- 7
matAGEmale <- 7
############# sexual maturity

SplitBySex <- split(x = info, f = info$Sex)
### Applies the LLin3 transformation to the input vector x
fun_llin_trans_Female <- Vectorize(function(x, maturity = matAGEfem) {
  if (is.na(x) | is.na(maturity)) {return(NA)}
  k <- 1.5
  y <- 0
  if (x < maturity) {y = log((x+k)/(maturity+k))}
  else {y = (x-maturity)/(maturity+k)}
  return(y)
})

fun_llin_trans_Male <- Vectorize(function(x, maturity = matAGEmale) {
  if (is.na(x) | is.na(maturity)) {return(NA)}
  k <- 1.5
  y <- 0
  if (x < maturity) {y = log((x+k)/(maturity+k))}
  else {y = (x-maturity)/(maturity+k)}
  return(y)
})

SplitBySex$Male$tAge = fun_llin_trans_Male(SplitBySex$Male$Age)
SplitBySex$Female$tAge = fun_llin_trans_Female(SplitBySex$Female$Age)


#inverse sex-specific

fun_llin_inv_Female <- Vectorize(function(y, maturity = matAGEfem) {
  if (is.na(y) | is.na(maturity)) {return(NA)}
  k <- 1.5
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y)-k}
  else {x = (maturity+k)*y+maturity}
  return(x)
})

fun_llin_inv_Male <- Vectorize(function(y, maturity = matAGEmale) {
  if (is.na(y) | is.na(maturity)) {return(NA)}
  k <- 1.5
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y)-k}
  else {x = (maturity+k)*y+maturity}
  return(x)
})

SplitBySex$Female$iAge = fun_llin_inv_Female(SplitBySex$Female$tAge)
SplitBySex$Male$iAge = fun_llin_inv_Male(SplitBySex$Male$tAge)

info <- unsplit(SplitBySex, info$Sex)


################################################################################

#Using alpha=ALPH
ALPH = 0.5

#Create vectors ofor elastic net model construction

epi<-data.matrix(epi)
epi_data<-t(epi[,1:ncol(epi)])

age_data<-info$tAge[(1:nrow(info))]

#create weighting vector
info$AC_rev <- ifelse(info$AgeConfidence <= 365, 365, info$AgeConfidence)
unique(info$AC_rev)
info$weight <- 1/(info$AC_rev/365)
weight <- info$weight


#........................................................
#
# elastic net model including sex-specific maturity
#........................................................
#Before you really start, you will generate a lambda value using the "cv.glmnet()" function (you can use any value of "nfold" here, but we use a value of 10 here, because anything larger won't make that much of a difference). 
glmnet.Training.CV<-cv.glmnet(x=as(epi_data, "dgCMatrix"), y=as.vector(age_data),nfolds=nrow(info),alpha=ALPH,family="gaussian", weights = weight)

# Save this value of lambda to "lambda.training. NOTE: I use "lambda.1se" instead of "lambda.min" for all of my pipelines and clocks, because it reduces overfitting and makes performance more consistent in new data sets
lambda.training=glmnet.Training.CV$lambda.1se


# LOOC 
for (i in 1:nrow(info)){
  epi_data_loo <- epi_data[-i,]
  age_data_loo <-age_data[-i]
  weight_loo <- weight[-i]
  test_id <- as.matrix(epi_data[i,])
  test_id <-t(test_id)
  
  fit1<-glmnet(x=epi_data_loo, y=age_data_loo,lambda = lambda.training,alpha=ALPH,family="gaussian", weights = weight_loo)
  y.fitted <- predict(fit1,test_id, type = "response")
  info$DNAmoClockLOO[i] <- y.fitted
  print(i)
}


#split based on Sex 
SplitBySex <- split(x = info, f = info$Sex)

SplitBySex$Female$invAgePredictedLOOCV = fun_llin_inv_Female(SplitBySex$Female$DNAmoClockLOO)
SplitBySex$Male$invAgePredictedLOOCV = fun_llin_inv_Male(SplitBySex$Male$DNAmoClockLOO)

info <- unsplit(SplitBySex, info$Sex)


#error LOOCV
errorLOOCV <- info$Age - info$invAgePredictedLOOCV
MAE_LOOCV <- sqrt(median(errorLOOCV^2))
MAE_LOOCV

# pearson correlation LOOCV
r_LOOCV <- cor(info$Age, info$invAgePredictedLOOCV, method = "pearson")
r_LOOCV




#plot LOOCV
Plot = ggplot(info, aes(x=Age, y = invAgePredictedLOOCV)) +
  ylim(0,60) + xlim(0,40) +
  geom_point(shape=1, size=2) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, colour="deepskyblue4") +
  geom_abline(intercept = 0, slope = 1, linetype = 3)+
  theme_classic()+
  theme(axis.title=element_blank(), axis.text = element_text(size=14, colour = "black"),);Plot


#run clock with all samples
  
final_clock<-glmnet(x=epi_data, y=age_data,lambda = lambda.training,alpha=ALPH,family="gaussian")



##################################################################################################################
##################################################################################################################
#
#   Sex prediction for samples
#
##################################################################################################################
##################################################################################################################

#read in randomised data
sex_info <- readRDS("data/info_sex_shuff.rds")
epi_all <- readRDS("data/epi_sex_shuff.rds")
epiAll <- epi_all[1:37492,]

fem <-which(sex_info$Sex == "Female")
mal <-which(sex_info$Sex == "Male")

#make binomial glmnet

#Using alpha=ALPH
ALPH = 0.5

epiAll<-data.matrix(epiAll)

#Create vector for elastic net model construction
sex_data<-t(epiAll[,1:ncol(epiAll)])

sex_age<-sex_info$Sex[(1:nrow(sex_info))]

sponge_data <-sex_info$Sponger[(1:nrow(sex_info))]

#10-fold cross validation
glmnet.Training.CV.Sex<-cv.glmnet(x=as(sex_data, "dgCMatrix"), y=as.vector((sex_age)),nfolds=10,alpha=ALPH,family="binomial")
plot(glmnet.Training.CV.Sex)

lambda.training.sex=glmnet.Training.CV.Sex$lambda.1se


# LOOC non-specific
for (i in 1:nrow(sex_info)){
  epi_data_loo <- sex_data[-i,]
  age_data_loo <-sex_age[-i]
  test_id <- as.matrix(sex_data[i,])
  test_id <-t(test_id)
  
  fit1<-glmnet(x=epi_data_loo, y=age_data_loo,lambda = lambda.training.sex,alpha=ALPH,family="binomial")
  y.fitted <- predict(fit1,test_id, type = "response")
  sex_info$SexClockLOO[i] <- y.fitted
  print(i)
}


glmnet.SEX = glmnet(sex_data,(sex_age), family = "binomial", alpha = 0.5, lambda = lambda.training.sex)
plot(glmnet.SEX)

sex_info$DNAmSex = predict(glmnet.SEX,sex_data,type = "response")
