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
cat(paste("Training data loaded: ",n," sources\n"))

cat("Final training set size: ",dim(featureCombinedTrain )[1],"\n") 
cat("Number of unique labels: ",length(unique(labelCombinedTrain)), "\n")

classType <- levels(labelCombinedTrain);

listLabels = GenerateLabels(labelCombinedTrain)
pulsatingLabels <- listLabels[[1]];
eruptiveLabels <- listLabels[[2]];
multiStarLabels <- listLabels[[3]];
unknownClassLabels <- listLabels[[4]];

pca3 = PCA(featureCombinedTrain , graph = FALSE)

plot(pca3$ind$coord[1:length(eruptiveLabels),1], 
pca3$ind$coord[1:length(eruptiveLabels),2],
ylab = "PCA 2", xlab = "PCA 1", pch=1)
points(pca3$ind$coord[pulsatingLabels== 1,1],pca3$ind$coord[pulsatingLabels== 1,2], pch=16, col = 2)
points(pca3$ind$coord[multiStarLabels == 1,1],pca3$ind$coord[multiStarLabels == 1,2], pch=16, col = 3)
points(pca3$ind$coord[unknownClassLabels== 1,1],pca3$ind$coord[unknownClassLabels== 1,2], pch=16, col = 4)
points(pca3$ind$coord[eruptiveLabels== 1,1],pca3$ind$coord[eruptiveLabels== 1,2], pch=16, col = 5)
legend("bottomright", legend = c("pulsating","multiStar", "other",
"eruptive"), pch = c(16,16,16,16), col = c(2,3,4,5), cex = 0.8)
grid()

####################################################
# Load TESTING data
cat("Loading LINEAR TESTING data\n")
#####

load("data/LINEAR_Dataset.R")
# featLINEAR
numOfLINEAR = dim(featLINEAR)[1];
####################################################
####################################################
## Sub-Function Calls

source("src/astro/utils/classification/GenerateMLPAndROCCurvesOptimum.R")
source("src/astro/utils/GenerateTrainingAndCrossval.R")

####################################################
# Train Pulsating Classifier

alphaCrit <- sort(seq(from = -0.01, to = 1.01, by = 0.01), decreasing = TRUE);

rocAucArray <- array(dim=c(1, length(classType)));
prAucArray <- array(dim=c(1, length(classType)));

predTestProb <- array(dim=c(length(classType), numOfLINEAR));

tmpFalseRateArray <- array(dim=c(length(alphaCrit), length(classType)));
tmpHitRateArray <- array(dim=c(length(alphaCrit), length(classType)));
tmpPrecisionRateArray <- array(dim=c(length(alphaCrit), length(classType)));

for(i in 1:length(classType)){

	tmpLabels <- array(dim=c(1,length(labelCombinedTrain)))

	tmpLabels [labelCombinedTrain == classType[i]] <- +1;
	tmpLabels [labelCombinedTrain != classType[i]] <- -1;

	#numClass <- length(tmpLabels [labelCombinedTrain == classType[i]] );
	#cat("Number in class: ", numClass, "\n");

	#tmpTargets <- decodeClassLabels(tmpLabels)
	#featureTmp <- splitForTrainingAndTest(featureCombinedTrain , tmpTargets , ratio = 0.5)

	listTrainCross = GenerateTrainingAndCrossval(featureCombinedTrain , tmpLabels);

	trainingList = data.frame(as.vector(listTrainCross[[3]]), 
		 as.vector(listTrainCross[[3]]));
	
	trainingList [trainingList [1:length(listTrainCross[[3]]),1] == 1,1] = 0
	trainingList [trainingList [1:length(listTrainCross[[3]]),2] == -1,2] = 0 
	trainingList [trainingList [1:length(listTrainCross[[3]]),1] == -1,1] = 1 
	colnames(trainingList) <- c("-1", "+1");

	crossList = data.frame(as.vector(listTrainCross[[4]]), 
		 as.vector(listTrainCross[[4]]));

	crossList [crossList [1:length(listTrainCross[[4]]),1] == 1,1] = 0
	crossList [crossList [1:length(listTrainCross[[4]]),2] == -1,2] = 0 
	crossList [crossList [1:length(listTrainCross[[4]]),1] == -1,1] = 1 
	colnames(crossList ) <- c("-1", "+1");

	featureTmp = list(inputsTrain= listTrainCross[[1]], 
		inputsTest = listTrainCross[[2]], 
		targetsTrain = trainingList,
		targetsTest = crossList)

	performanceCurves = GenerateMLPAndROCCurvesOptimum(
		featureTmp , classType[i]);

	tmpFalseRateArray[,i]= performanceCurves[[1]];
	tmpHitRateArray[,i]= performanceCurves[[2]];
	tmpPrecisionRateArray[,i]= performanceCurves[[3]];
	tmpClassifier = performanceCurves[[4]];
	
	# for(j in 1:numOfLINEAR){ 
	#	sampleLinear <- as.numeric(featLINEAR[j,])
	#	probTemp  <- predict(tmpClassifier , sampleLinear )
	#	predTestProb[i,j]  <- probTemp[1,2];
	#}

	#critIndex <- which.max(tmpFalseRate  > 0.002)
	
	#alphaCrit <- sort(seq(from = -0.01, to = 1.01, by = 0.01), decreasing = TRUE);

	#numInClass = length(predTestProb[2,predTestProb[2,] > alphaCrit[critIndex]])
}

save(tmpFalseRateArray, tmpHitRateArray, tmpPrecisionRateArray,
  file = "data/Alpha_MLP.R")