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

library(RSNNS)

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
source("src/astro/utils/classification/GenerateMLPAndROCCurves.R")

#####################################################################

listLabels = GenerateLabels(class.train);

pulsatingLabels <- listLabels[[1]];
eruptiveLabels <- listLabels[[2]];
multiStarLabels <- listLabels[[3]];
unknownClassLabels <- listLabels[[4]];

####################################################
# Train Pulsating Classifier

pulsatingTargets <- decodeClassLabels(pulsatingLabels)
featurePulse <- splitForTrainingAndTest(feat.train, pulsatingTargets , ratio = 0.5)

cat("Fitting MLP Pulsating\n")

performanceCurves = GenerateMLPAndROCCurves(
	featurePulse , pulsatingTargets ,
	"Pulsating");

pulsatingFalseRate = performanceCurves[[1]];
pulsatingHitRate = performanceCurves[[2]];
pulsatingPrecisionRate = performanceCurves[[3]];

####################################################
# Train Eruptive Classifier

eruptiveTargets <- decodeClassLabels(eruptiveLabels)
featureEruptive <- splitForTrainingAndTest(feat.train, eruptiveTargets , ratio = 0.5)

cat("Fitting MLP eruptive\n")

performanceCurves = GenerateMLPAndROCCurves(
	featureEruptive , eruptiveTargets ,
	"Eruptive");

eruptiveFalseRate = performanceCurves[[1]];
eruptiveHitRate = performanceCurves[[2]];
eruptivePrecisionRate = performanceCurves[[3]];

####################################################
# Train MultiStar Classifier
multiStarTargets <- decodeClassLabels(multiStarLabels )
featureMulti <- splitForTrainingAndTest(feat.train, multiStarTargets , ratio = 0.5)

cat("Fitting MLP multiStar\n")

performanceCurves = GenerateMLPAndROCCurves(
	featureMulti , multiStarTargets ,
	"Multi-Star");

multiStarFalseRate = performanceCurves[[1]];
multiStarHitRate = performanceCurves[[2]];
multiStarPrecisionRate = performanceCurves[[3]];


####################################################
# Train Other Classifier
listTrainCrossUnknown = GenerateTrainingAndCrossval(
	feat.train, unknownClassLabels);

unknownTargets <- decodeClassLabels(unknownClassLabels)
featureUnknown <- splitForTrainingAndTest(feat.train, unknownTargets , ratio = 0.5)

cat("Fitting MLP Unknown\n")

performanceCurves = GenerateMLPAndROCCurves(
	featureUnknown , multiStarTargets ,
	"Unknown");

otherFalseRate = performanceCurves[[1]];
otherHitRate = performanceCurves[[2]];
otherPrecisionRate = performanceCurves[[3]];
####################################################
####################################################
