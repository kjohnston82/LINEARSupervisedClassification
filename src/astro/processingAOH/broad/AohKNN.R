#####################################################################
##### run this script to perform AL from start to finish on ASAS data
## and save all relevant output, classifications + performance metrics
####################################################
# load needed packages and files
set.seed(1)

setwd("C:/Users/kyle.johnston/Desktop/Dissertation/starvars/")

source("src/astro/utils/utils_classify.R")

####################################################
####################################################
# load data

library(foreign)

####################################################
# Load TRAINING data
cat("Loading ASAS + Hipp + OGLE Training data\n")
#####

feat.train = read.table("data/trainingsets/training_set_features.dat",header=TRUE,sep=',')
ID.train = feat.train$ID.train

## Trim the patterns
exclude = c("ID.train","qso_log_chi2_qsonu","qso_log_chi2nuNULL_chi2nu","freq_n_alias")
feat.train = feat.train[,-which(names(feat.train) %in% exclude)]
ra.dec.train = read.table("data/trainingsets/training_set_ra_dec.dat",header=TRUE,sep=',')
class.train = read.table("data/trainingsets/training_set_id_class.dat",header=FALSE,sep='\t')[,2]
n = length(ID.train)

####################################################
# Train Classifier on all classes and features
cat(paste("Training data loaded: ",n," sources\n"))

cat("Final training set size: ",dim(feat.train)[1],"\n") 
cat("Number of unique sources: ",length(unique(ID.train)),"\n")
cat("Number of unique labels: ",length(unique(class.train)), "\n")

####################################################
## Sub-Function Calls

source("src/astro/utils/GenerateTrainingAndCrossval.R")
source("src/astro/utils/GenerateLabels.R")
source("src/astro/utils/classification/GeneratekNNAndROCCurves.R")

#####################################################################

listLabels = GenerateLabels(class.train);

pulsatingLabels <- listLabels[[1]];
eruptiveLabels <- listLabels[[2]];
multiStarLabels <- listLabels[[3]];
unknownClassLabels <- listLabels[[4]];

####################################################
# Train Pulsating Classifier
listTrainCrossPulse = GenerateTrainingAndCrossval(feat.train, pulsatingLabels);

pulsatingTrainingData = listTrainCrossPulse[[1]]
pulsatingCrossvalData = listTrainCrossPulse[[2]] 
pulsatingTrainingLabel = listTrainCrossPulse[[3]] 
pulsatingCrossvalLabel = listTrainCrossPulse[[4]] 

cat("Fitting kNN Pulsating\n")

performanceCurves = GeneratekNNAndROCCurves(
	pulsatingTrainingData, pulsatingTrainingLabel,
 	pulsatingCrossvalData, pulsatingCrossvalLabel,
	"Pulsating");

pulsatingFalseRate = performanceCurves[[1]];
pulsatingHitRate = performanceCurves[[2]];
pulsatingPrecisionRate = performanceCurves[[3]];

####################################################
# Train Eruptive Classifier
listTrainCrossEruptive = GenerateTrainingAndCrossval(
	feat.train, eruptiveLabels);

eruptiveTrainingData = listTrainCrossEruptive[[1]] 
eruptiveCrossvalData = listTrainCrossEruptive[[2]]
eruptiveTrainingLabel = listTrainCrossEruptive[[3]] 
eruptiveCrossvalLabel = listTrainCrossEruptive[[4]] 

cat("Fitting kNN Eruptive\n")

performanceCurves = GeneratekNNAndROCCurves(
	eruptiveTrainingData , eruptiveTrainingLabel ,
 	eruptiveCrossvalData , eruptiveCrossvalLabel ,
	"Eruptive");

eruptiveFalseRate = performanceCurves[[1]];
eruptiveHitRate = performanceCurves[[2]];
eruptivePrecisionRate = performanceCurves[[3]];

####################################################
# Train MultiStar Classifier
listTrainCrossMultiStar = GenerateTrainingAndCrossval(
	feat.train, multiStarLabels );

multiStarTrainingData = listTrainCrossMultiStar[[1]] 
multiStarCrossvalData = listTrainCrossMultiStar[[2]]
multiStarTrainingLabel = listTrainCrossMultiStar[[3]] 
multiStarCrossvalLabel = listTrainCrossMultiStar[[4]] 

cat("Fitting kNN MultiStar\n")

performanceCurves = GeneratekNNAndROCCurves(
	multiStarTrainingData , multiStarTrainingLabel ,
 	multiStarCrossvalData , multiStarCrossvalLabel ,
	"MultiStar");

multiStarFalseRate = performanceCurves[[1]];
multiStarHitRate = performanceCurves[[2]];
multiStarPrecisionRate = performanceCurves[[3]];


####################################################
# Train Other Classifier
listTrainCrossUnknown = GenerateTrainingAndCrossval(
	feat.train, unknownClassLabels);

unknownTrainingData = listTrainCrossUnknown[[1]] 
unknownCrossvalData = listTrainCrossUnknown[[2]]
unknownTrainingLabel = listTrainCrossUnknown[[3]] 
unknownCrossvalLabel = listTrainCrossUnknown[[4]] 

cat("Fitting kNN Unknown\n")

performanceCurves = GeneratekNNAndROCCurves(
	unknownTrainingData , unknownTrainingLabel ,
 	unknownCrossvalData , unknownCrossvalLabel ,
	"Unknown");

otherFalseRate = performanceCurves[[1]];
otherHitRate = performanceCurves[[2]];
otherPrecisionRate = performanceCurves[[3]];
####################################################
####################################################
