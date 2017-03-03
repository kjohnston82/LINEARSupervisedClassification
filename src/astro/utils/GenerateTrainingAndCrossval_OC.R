GenerateTrainingAndCrossval_OC <- function(featTrain, labels, split){

originalFeatureSpace = as.matrix(featTrain);
nOfFeatures = dim(originalFeatureSpace);
cat("Number of Features:", nOfFeatures[2] ,"\n")

indexTrainingArray <- array(dim=c(1,length(labels)))
indexCrossvalArray <- array(dim=c(1,length(labels)))
counterTrain = 1;
counterCrossval =1;

# 50/50 split of features for OC
for(i in 1:length(labels)){
	x1 <- runif(1, 0.0, 1.0)

	if(x1 > split){
		indexTrainingArray[counterTrain] <- i;
		counterTrain = counterTrain + 1; 
	}else{
		indexCrossvalArray [counterCrossval ] <- i;
		counterCrossval = counterCrossval + 1;
	}
}


counterTrain <- counterTrain - 1;
counterCrossval <- counterCrossval - 1;

trainingData <- array(dim=c(counterTrain, nOfFeatures[2]))
crossvalData <- array(dim=c(counterCrossval, nOfFeatures[2]))
trainingLabel <- vector(length = counterTrain)
crossvalLabel <- vector(length =  counterCrossval)

counter = 1;
for(i in indexTrainingArray[1:counterTrain]){
	trainingData [counter, 1:nOfFeatures[2]] = 
		originalFeatureSpace[i, 1:nOfFeatures[2]];

	# all one label
	trainingLabel [counter] = +1;

	counter = counter + 1;
}

counter = 1;
for(i in indexCrossvalArray [1:counterCrossval]){
	crossvalData [counter, 1:nOfFeatures[2]] = 
		originalFeatureSpace[i, 1:nOfFeatures[2]];

	# all one label
	crossvalLabel [counter] = +1;
	counter = counter + 1;
}

return(list(trainingData, crossvalData , trainingLabel , crossvalLabel ))

}