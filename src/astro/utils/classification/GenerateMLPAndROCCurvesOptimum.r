GenerateMLPAndROCCurvesOptimum <- function(
	feature, summaryType) {
  
feature <- normTrainingAndTestSet(feature)

alphaCrit <- sort(seq(from = -0.01, to = 1.01, by = 0.01), decreasing = TRUE);

hitRate  		= array(dim=c(1,length(alphaCrit)));
falseRate 		= array(dim=c(1,length(alphaCrit)));
precisionRate 	= array(dim=c(1,length(alphaCrit)));

mlpModel <- mlp(feature$inputsTrain, feature$targetsTrain, size = 60 ,
	learnFuncParams = 0.5, maxit = 60, inputsTest = feature$inputsTest,
	targetsTest = feature$targetsTest);

pred.cross.prob <- predict(mlpModel , feature$inputsTest)
	
dimenSet <- dim(pred.cross.prob);

## Loop Over to determine ROC curve based on performance (probs
counter = 1;
for(i in 1:length(alphaCrit)){

	countTP = 0;
	countFP = 0;
	countTN = 0;
	countFN = 0;

	for(j in 1:dimenSet[1]){
		if(pred.cross.prob[j,2] > alphaCrit[i]){
			if(feature$targetsTest[j,2] == 1){
				countTP = countTP + 1;
			}else{
				countFP  = countFP  + 1;
			} 
		}else{
			if(feature$targetsTest[j,2] == 0){
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

frame()
par(mfrow = c(1,2))
plot(falseRate[1,1:length(alphaCrit)], hitRate[1,1:length(alphaCrit)],
	xlab = "False Rate", ylab = "Hit Rate",
	xlim = c(0.0,1.0), ylim = c(0.0,1.0),
	col=2, lty = 1, type = "l")
lines(xdivide , xdivide , col=1, lty = 2)
grid()

plot(hitRate[1,1:length(alphaCrit)], precisionRate[1,1:length(alphaCrit)],
	xlab = "Recall", ylab = "Precision", 
	xlim = c(0.0,1.0), ylim = c(0.0,1.0),
	col=2, lty = 2, type = "l")
grid()

####################################################
## Note, these libraries must be downloaded, :then: loaded via "library"
library(caTools)

#########################
##   Non-Uniform Trapezoidal Rule Numerical Integration
ROC.AUC = trapz(falseRate[1,1:length(alphaCrit)], hitRate[1,1:length(alphaCrit)]);
PR.AUC = trapz(hitRate[1,1:length(alphaCrit)],precisionRate[1,1:length(alphaCrit)]);

ROC.AUC = abs(ROC.AUC)
PR.AUC = abs(PR.AUC)

#########################
return(list(falseRate, hitRate, precisionRate, mlpModel));

}