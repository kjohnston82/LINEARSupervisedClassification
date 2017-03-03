#####################################################################
##### run this script to perform AL from start to finish on ASAS data
## and save all relevant output, classifications + performance metrics
####################################################
# load needed packages and files
set.seed(1)

setwd("C:/Users/kyle.johnston/Desktop/Dissertation/starvars/")

source("src/astro/utils/utils_classify.R")
source("src/astro/utils/GenerateLabels.R")

####################################################
####################################################
# load data
library(randomForest)
library(nnet)
library(caTools)
library(RSNNS)
library(caTools)
library(FactoMineR)
library(kernlab)
####################################################
# Load TRAINING data
cat("Loading ASAS + Hipp + OGLE +LINEAR Training data\n")
#####

load("data/Training_Dataset_AohLINEAR.R")
# featureCombinedTrain = trainingList[[1]];
# labelCombinedTrain = trainingList[[2]]


featureCombinedTrain$color_diff_bj <- 
	as.numeric(featureCombinedTrain$color_diff_bj);

featureCombinedTrain$color_diff_hk <- 
	as.numeric(featureCombinedTrain$color_diff_hk);

featureCombinedTrain$color_diff_jh <- 
	as.numeric(featureCombinedTrain$color_diff_jh);

featureCombinedTrain$color_diff_rj <- 
	as.numeric(featureCombinedTrain$color_diff_rj);

featureCombinedTrain$color_diff_vj <- 
	as.numeric(featureCombinedTrain$color_diff_vj);

# Train Classifier on all classes and features
cat("Final training set size: ",dim(featureCombinedTrain )[1],"\n") 
cat("Number of unique labels: ",length(unique(labelCombinedTrain)), "\n")

classType <- levels(labelCombinedTrain);
########################################################
load("data/LINEAR_Dataset.R")
# is in the load, IDLinearReduced 
########################################################
maxFalseAlarm = 0.005;

name <- load("data/ResultsFromLINEAR/PredProbOutput_MLP.R")
load("data/ResultsFromLINEAR/PredProbOutput_MLP.R")
predTestProbMLP <- predTestProb;

name <- load("data/Alpha_MLP.R")
load("data/Alpha_MLP.R")

fr_MLP <- tmpFalseRateArray;
hr_MLP <- tmpHitRateArray;
p_MLP <- tmpPrecisionRateArray;

sizeClasses <- dim(predTestProbMLP)[1];
sizeLINEAR <- dim(predTestProbMLP)[2];

alphaCrit <- sort(seq(from = -0.01, to = 1.01, by = 0.01), decreasing = TRUE);
thresholdArray <- array(dim=c(1, sizeClasses));

for(i in 1:sizeClasses){
	# which alpha 
	thresholdArray [i] <- alphaCrit[which.max(fr_MLP[,i] > maxFalseAlarm)];

	ofClassType <- IDLinearReduced[predTestProbMLP[i,] > thresholdArray[i]];

	classLabel <- strsplit(classType[i], " ")[[1]]
	classLabel <- paste("data/",classLabel[1], "StarType_EstMembers_MLP", sep = "")

	save(ofClassType, file = paste(classLabel, "R", sep = ".")) 
}