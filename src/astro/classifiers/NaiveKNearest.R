NaiveKNearest <- function(
	trainingData, grpSource_Training, 
	crossvalData, intKN, p, fltSpread,  strMissValue){

# function [fltResponse, classEstimate] = Naive_K_Nearest(...
#     trainingData, grpSource_Training, crossVal, ...
#     intKN, p, fltSpread,  strMissValue)
#
# Author: Kyle Johnston     22 Jan 2015
#
# Usage: This function generates class estimates based on the K-NN
#   algorithm. Ties are represented by the "strMissValue" string in the
#   output labels. This particular implimentation is a weighted K-NN
#   algorithm based on the distance from the input pattern to other
#   training pattterns.
# 
# Input: 
#       trainingData: the training dataset
#       grpSource_Training: the labels for the training dataset
#       crossVal: the test dataset
#       intKN: value of k
#       p: distance type (2 is quadratic distance)
#       fltSpread: for the distance measurement
#       strMissValue: strnig label for rejected points (ties)
#
# Output:
#       fltResponse: the set of posterior probabilities
#       classEstimate: estimated labels for fltPatternArray_CrossVal

# Initialize Variables and Arrays
dimTraining <- dim(trainingData);
intNumberOfDimensionsTraining <- dimTraining[2];
intNumberOfTrainingPatterns <- dimTraining[1];

dimCrossval <- dim(crossvalData);
intNumberOfDimensionsTest <- dimCrossval[2];
intNumberOfTestPatterns <- dimCrossval[1];

grpUniqueTypes <- sort(unique(grpSource_Training));
intNumberOfUnique <- length(grpUniqueTypes);

fltResponse <- matrix(0, intNumberOfTestPatterns, intNumberOfUnique);

#########################
# Scale the inputs for standardization
scaleMean <- colMeans(trainingData);

sum_pattern2 <- colSums(trainingData^2)/intNumberOfTrainingPatterns;
sumpattern_2 <- (colSums(trainingData)/intNumberOfTrainingPatterns)^2;

scaleStdDev <- sqrt((sum_pattern2 - sumpattern_2));

fltPatternArray_Training <- matrix(0, intNumberOfTrainingPatterns, intNumberOfDimensionsTraining);
fltPatternArray_CrossVal <- matrix(0, intNumberOfTestPatterns , intNumberOfDimensionsTraining);

for(i in 1:intNumberOfDimensionsTraining){
	fltPatternArray_Training[1:intNumberOfTrainingPatterns,i] <- 
			(trainingData[1:intNumberOfTrainingPatterns,i] - scaleMean[i])/scaleStdDev[i] ;
	fltPatternArray_CrossVal[1:intNumberOfTestPatterns,i] <- 
			(crossvalData[1:intNumberOfTestPatterns,i]- scaleMean[i])/scaleStdDev[i] ;
}

#########################
# Check to make sure two datasets operate in the same vector space
if(intNumberOfDimensionsTest == intNumberOfDimensionsTraining){

    # Initialize distance arrays
    fltDistanceSet <- matrix(0,intNumberOfTrainingPatterns,
		intNumberOfDimensionsTest);
    fltDistance <- rep(0, intNumberOfTrainingPatterns);

    # Loop through test patterns provided
    for(i in 1:intNumberOfTestPatterns){

        # Compute distance between training and test
        for(j in 1:intNumberOfTrainingPatterns){

		set <- c(1:intNumberOfDimensionsTest);

            fltDistanceSet[j,set] =(abs(fltPatternArray_CrossVal[i,set] - fltPatternArray_Training[j,set]))^p;
          
            fltDistance[j] = (sum(fltDistanceSet[j,1:intNumberOfDimensionsTest]))^(1/p);

        }

        # Find the K closest points
        for( j in 1:intKN){
            indexCurrent = which.min(fltDistance);

            for( k in 1:intNumberOfUnique){
                if(grpUniqueTypes[k] == grpSource_Training[indexCurrent[1]]){

		        # Unweighted scheme of associated response
		         fltResponse[i,k] = fltResponse[i,k] + 1.0;

			  # Weighted scheme of associated based on distance
                     #fltResponse[i,k] = fltResponse[i,k] + 1.5 - pnorm(fltDistance[indexCurrent], mean = 0 , sd = fltSpread);
                    break;
                }
            }
	
            fltDistance <- fltDistance[min(fltDistance) != fltDistance];
        }

        # Compute posterior probabilities
        fltResponse[i,1:intNumberOfUnique] = 
		fltResponse[i,1:intNumberOfUnique]/sum(fltResponse[i,1:intNumberOfUnique]);
    }

}

# Label the datasets based on posterior probabilities
classEstimate<- array(dim = c(intNumberOfUnique,1));
for( i in 1:intNumberOfTestPatterns){

    fltMaxProb = max(fltResponse[i,1:intNumberOfUnique]);
    intSingle = length(fltResponse[i,fltResponse[i,1:intNumberOfUnique] == fltMaxProb])
    index = match(fltMaxProb, fltResponse[i,1:intNumberOfUnique]);

    if(intSingle == 1){
      classEstimate[i] = grpUniqueTypes[index];
    }else{
	classEstimate[i] = strMissValue;
    }
}

return(list(fltResponse, classEstimate));

}