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
library(e1071)
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

source("src/astro/utils/classification/GenerateSVMAndROCCurvesOptimum.R")
source("src/astro/utils/GenerateTrainingAndCrossval_OC.R")
source("src/astro/utils/GenerateTrainingAndCrossval.R")

####################################################
# Train OC_SVM Classifier

logicalAnomaly <- vector(mode = "logical", 
	length = length(numOfLINEAR))
tmpLabels <- array(+1.0, dim=c(1,length(labelCombinedTrain)))


svmModel<-svm(as.matrix(featureCombinedTrain) ,y=as.vector(tmpLabels),
               type='one-classification',
               nu=(0.03),
               scale=TRUE,
               kernel="radial")

for(j in 1:numOfLINEAR){
	sampleLinear <- t(as.matrix(as.numeric(featLINEAR[j,])))

	logicalAnomaly[j] <- predict(svmModel, newdata = sampleLinear, 
		probability = TRUE);
}

save(logicalAnomaly, file = "data/PredProb_SVM_OC.R")