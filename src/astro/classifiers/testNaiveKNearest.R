########################################################
setwd("C:/Users/kyle.johnston/Desktop/Dissertation/starvars/")

source("src/astro/utils/GenerateTrainingAndCrossval.R")

######################################

dimI <- dim(iris);

featTrain <- iris[1:dimI[1],1:(dimI[2] - 1)];
classLabels <- iris[1:dimI[1], dimI[2]];
newLabels<- array(dim = c(1,length(classLabels)));

labelTypes <- unique(classLabels);
boolFirstType <- classLabels == labelTypes[1];

newLabels[boolFirstType] <- +1;
newLabels[!boolFirstType] <- -1;

listTrainCross = GenerateTrainingAndCrossval(
	featTrain, newLabels);

trainingData = listTrainCross [[1]]
crossvalData = listTrainCross [[2]] 
trainingLabel = listTrainCross [[3]] 
crossvalLabel = listTrainCross [[4]] 

########################################################

source("src/astro/classifiers/NaiveKNearest.R")


######################################

listKNN <- Naive_K_Nearest(trainingData, as.vector(trainingLabel), crossvalData,
      20, 2.0, 1.0,  "Miss");

fltResponse = listKNN[[1]]
classEstimate = listKNN[[2]] 


plot(crossvalData[1:length(classEstimate),1],
	crossvalData[1:length(classEstimate),2], col = 1, pch = 16,
	xlim = c(4,9), ylim = c(1.5, 5),
	xlab = "Feature 1", ylab = "Feature 2")
points(
	crossvalData[(classEstimate == 1 & crossvalLabel == 1),1],
	crossvalData[(classEstimate == 1 & crossvalLabel == 1),2], col = 1, pch = 1,
	xlab = "Feature 1", ylab = "Feature 2")
points(
	crossvalData[(classEstimate == 1 & crossvalLabel == -1),1],
	crossvalData[(classEstimate == 1 & crossvalLabel == -1),2], col = 2, pch = 1)
points(
	crossvalData[(classEstimate == -1 & crossvalLabel == -1),1],
	crossvalData[(classEstimate == -1 & crossvalLabel == -1),2], col = 3, pch = 1)
points(
	crossvalData[(classEstimate == -1 & crossvalLabel == 1),1],
	crossvalData[(classEstimate == -1 & crossvalLabel == 1),2], col = 4, pch = 1)
grid()
legend("bottomright", legend = c("TP", "FP", "TN", "FN"), col = c(1,2,3,4), pch = c(1,1,1,1))

###