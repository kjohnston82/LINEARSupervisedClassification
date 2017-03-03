GenerateSVMAndROCCurvesOptimum <- function(
	Training.data, Training.label, Crossval.data, Crossval.label,
	summaryType) {
  
kernelSpread <- 0.1

alphaCrit <- sort(seq(from = -0.01, to = 1.01, by = 0.01), decreasing = TRUE);

hitRate  		= array(dim=c(1,length(alphaCrit)));
falseRate 		= array(dim=c(1,length(alphaCrit)));
precisionRate 	= array(dim=c(1,length(alphaCrit)));

svmTrainingLabel = vector(length = length(Training.label));
svmCrossvalLabel = vector(length = length(Training.label));

svmTrainingLabel <- Training.label;
svmCrossvalLabel <- Crossval.label;


svmStar <- svm(x = as.matrix(Training.data),y = factor(svmTrainingLabel), 
	 kernel = "radial", gamma = 0.1, probability = TRUE)

svmEst <- predict(svmStar, newdata = Crossval.data, probability = TRUE)

pred.cross.prob<-attr(svmEst , "probabilities")
dimenSet <- dim(pred.cross.prob);

counter = 1;
for(i in 1:length(alphaCrit)){

	countTP = 0;
	countFP = 0;
	countTN = 0;
	countFN = 0;

	for(j in 1:dimenSet[1]){
		if(pred.cross.prob[j,2] < alphaCrit[i]){
			if(svmCrossvalLabel[j] == +1){
				countTP = countTP + 1;
			}else{
				countFP  = countFP  + 1;
			} 
		}else{
			if(svmCrossvalLabel[j] == -1){
				countTN = countTN + 1;
			}else{
				countFN = countFN + 1;
			} 
		}
		}

	hitRate[1,i]   		= countTP /(countTP + countFN);
	falseRate[1,i] 		= countFP /(countFP + countTN);
	precisionRate[1,i] 	= countTP /(countTP + countFP);
}

xdivide <- seq(from = 0, to = 1, by = 0.01);

precisionRate[is.nan(precisionRate)] = 0.9999999


####################################################
## Note, these libraries must be downloaded, :then: loaded via "library"
library(caTools)

#########################
##   Non-Uniform Trapezoidal Rule Numerical Integration
return(list(falseRate, hitRate, precisionRate, svmStar));

}