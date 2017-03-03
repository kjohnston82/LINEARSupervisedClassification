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
library(e1071)
library(caTools)

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

classType <- levels(class.train);
####################################################
####################################################
## Sub-Function Calls

source("src/astro/utils/classification/GenerateSVMAndROCCurves.R")
source("src/astro/utils/GenerateTrainingAndCrossval.R")
 
####################################################
# Train Pulsating Classifier

rocAucArray <- array(dim=c(1, length(classType)));
prAucArray <- array(dim=c(1, length(classType)));

for(i in 1:length(classType)){

	tmpLabels <- array(dim=c(1,length(class.train)))

	tmpLabels [class.train == classType[i]] <- +1;
	tmpLabels [class.train != classType[i]] <- -1;

	numClass <- length(tmpLabels [class.train == classType[i]] );
	cat("Number in class: ", numClass, "\n");

	listTrainCross = GenerateTrainingAndCrossval(feat.train, tmpLabels);

	tmpTrainingData = listTrainCross[[1]]
	tmpCrossvalData = listTrainCross[[2]] 
	tmpTrainingLabel = as.vector(listTrainCross[[3]]) 
	tmpCrossvalLabel = as.vector(listTrainCross[[4]]) 

	cat("Fitting SVM: ", classType[i]  ,"\n")

	performanceCurves = GenerateSVMAndROCCurves(
		tmpTrainingData , tmpTrainingLabel ,
 		tmpCrossvalData , tmpCrossvalLabel ,
		classType[i]);

	tmpFalseRate = performanceCurves[[1]];
	tmpHitRate = performanceCurves[[2]];
	tmpPrecisionRate = performanceCurves[[3]];
	
	alphaCrit <- sort(seq(from = -0.01, to = 1.01, by = 0.01), decreasing = TRUE);

	ROC.AUC = trapz(tmpFalseRate [1,1:length(alphaCrit)], tmpHitRate [1,1:length(alphaCrit)]);
	PR.AUC = trapz(tmpHitRate [1,1:length(alphaCrit)],tmpPrecisionRate [1,1:length(alphaCrit)]);

	rocAucArray[i] = abs(ROC.AUC)
	prAucArray[i] = abs(PR.AUC)

}

rocArraySVM <- rocAucArray;
prArraySVM <- prAucArray;

save(rocArrayKNN , prArrayKNN , file = "data/performance/individual/KNNPerformance.R")

plot(prAucArray, rocAucArray, xlab = "PR-AUC", ylab = "ROC-AUC", 
xlim = c(0,1), ylim = c(0,1), main = "Multi-Layer Perceptron", pch = 1, col = 1);
lines(alphaCrit ,alphaCrit, col=1, lty = 2)
points(prAucArray[c(20,21,22,24,27,28)], rocAucArray[c(20,21,22,24,27,28)], pch = 16, col = 2)
points(prAucArray[c(1,2,6)], rocAucArray[c(1,2,6)], pch = 16, col = 3)
points(prAucArray[7:9], rocAucArray[7:9], pch = 16, col = 4)
points(prAucArray[10:12], rocAucArray[10:12], pch = 16, col = 5)
points(prAucArray[13:19], rocAucArray[13:19], pch = 16, col = 6)
points(prAucArray[29:32], rocAucArray[29:32], pch = 16, col = 7)
points(prAucArray[c(3,4,5,23,25,26)], rocAucArray[c(3,4,5,23,25,26)], pch = 16, col = 8)
grid()
legend("bottomright", legend = c("Eruptive", "Giants", "Cephids", "RRLyr", "Other Pulse", "Multi-Star", "Others"), 
col = c(2:8), pch = 16, cex = 0.8, title = "Class Type")


plot(prAucArray, rocAucArray, xlab = "PR-AUC", ylab = "ROC-AUC", 
xlim = c(0,1), ylim = c(0,1), main = "Multi-Layer Perceptron", pch = 1, col = 1);
lines(alphaCrit ,alphaCrit, col=1, lty = 2)
points(prAucArray[c(20,21,22,24,27,28)], rocAucArray[c(20,21,22,24,27,28)], pch = 16, col = 2)
text(prAucArray[c(20,21,22,24,27,28)], rocAucArray[c(20,21,22,24,27,28)], labels = classType[c(20,21,22,24,27,28)], cex = 0.5, pos = 4)
grid()

plot(prAucArray, rocAucArray, xlab = "PR-AUC", ylab = "ROC-AUC", 
xlim = c(0,1), ylim = c(0,1), main = "Multi-Layer Perceptron", pch = 1, col = 1);
lines(alphaCrit ,alphaCrit, col=1, lty = 2)
points(prAucArray[c(1,2,6)], rocAucArray[c(1,2,6)], pch = 16, col = 3)
text(prAucArray[c(1,2,6)], rocAucArray[c(1,2,6)], labels = classType[c(1,2,6)], cex = 0.5, pos = 4)
grid()

plot(prAucArray, rocAucArray, xlab = "PR-AUC", ylab = "ROC-AUC", 
xlim = c(0,1), ylim = c(0,1), main = "Multi-Layer Perceptron", pch = 1, col = 1);
lines(alphaCrit ,alphaCrit, col=1, lty = 2)
points(prAucArray[7:9], rocAucArray[7:9], pch = 16, col = 4)
text(prAucArray[7:9], rocAucArray[7:9], labels = classType[7:9], cex = 0.5, pos = c(4,3,4))
grid()

plot(prAucArray, rocAucArray, xlab = "PR-AUC", ylab = "ROC-AUC", 
xlim = c(0,1), ylim = c(0,1), main = "Multi-Layer Perceptron", pch = 1, col = 1);
lines(alphaCrit ,alphaCrit, col=1, lty = 2)
points(prAucArray[10:12], rocAucArray[10:12], pch = 16, col = 5)
text(prAucArray[10:12], rocAucArray[10:12], labels = classType[10:12], cex = 0.5, pos = 4)
grid()

plot(prAucArray, rocAucArray, xlab = "PR-AUC", ylab = "ROC-AUC", 
xlim = c(0,1), ylim = c(0,1), main = "Multi-Layer Perceptron", pch = 1, col = 1);
lines(alphaCrit ,alphaCrit, col=1, lty = 2)
points(prAucArray[13:19], rocAucArray[13:19], pch = 16, col = 6)
text(prAucArray[13:19], rocAucArray[13:19], labels = classType[13:19], cex = 0.5, pos = 4)
grid()

plot(prAucArray, rocAucArray, xlab = "PR-AUC", ylab = "ROC-AUC", 
xlim = c(0,1), ylim = c(0,1), main = "Multi-Layer Perceptron", pch = 1, col = 1);
lines(alphaCrit ,alphaCrit, col=1, lty = 2)
points(prAucArray[29:32], rocAucArray[29:32], pch = 16, col = 7)
text(prAucArray[29:32], rocAucArray[29:32], labels = classType[29:32], cex = 0.5, pos = c(4,1,1,1))
grid()

plot(prAucArray, rocAucArray, xlab = "PR-AUC", ylab = "ROC-AUC", 
xlim = c(0,1), ylim = c(0,1), main = "Multi-Layer Perceptron", pch = 1, col = 1);
lines(alphaCrit ,alphaCrit, col=1, lty = 2)
points(prAucArray[c(3,4,5,23,25,26)], rocAucArray[c(3,4,5,23,25,26)], pch = 16, col = 8)
text(prAucArray[c(3,4,5,23,25,26)], rocAucArray[c(3,4,5,23,25,26)], labels = classType[c(3,4,5,23,25,26)], cex = 0.5, pos = c(4,4,4,4,4,4))
grid()

classPopSize <- array(dim=c(1, length(classType)));
for(i in 1:length(classType)){
	classPopSize [i] = length(tmpLabels [class.train == classType[i]]);
}

sizePop = 50;

plot(prAucArray, rocAucArray, xlab = "PR-AUC", ylab = "ROC-AUC", 
xlim = c(0,1), ylim = c(0,1), main = "Multi-Layer Perceptron", pch = 1, col = 1);
lines(alphaCrit ,alphaCrit, col=1, lty = 2)
points(prAucArray[classPopSize > sizePop ], rocAucArray[classPopSize > sizePop ], pch = 16, col = 7)
points(prAucArray[classPopSize <= sizePop ], rocAucArray[classPopSize <= sizePop ], pch = 16, col = 8)
legend("bottomright", legend = c("50<", "50 =>"), 
col = c(7,8), pch = 16, cex = 0.8, title = "# in Pop")
grid()

