GenerateKNNAndROCCurvesOptimum <- function(
	Training.data, Training.label, Crossval.data, Crossval.label,
	summaryType) {
  
#######################################################
setwd("C:/Users/kyle.johnston/Desktop/Dissertation/starvars/")

source("src/astro/classifiers/NaiveKNearest.R")

#######################################################

alphaCrit <- sort(seq(from = -0.01, to = 1.01, by = 0.01), decreasing = TRUE);

hitRate  		= array(dim=c(1,length(alphaCrit)));
falseRate 		= array(dim=c(1,length(alphaCrit)));
precisionRate 	= array(dim=c(1,length(alphaCrit)));

listKNN <- NaiveKNearest(Training.data, Training.label, 
	Crossval.data, 10, 1, 1.0,  "Miss");

fltResponse= listKNN[[1]]
classEstimate= listKNN[[2]] 

pred.cross.prob <- fltResponse
dimenSet <- dim(pred.cross.prob);

counter = 1;
for(i in 1:length(alphaCrit)){

	countTP = 0;
	countFP = 0;
	countTN = 0;
	countFN = 0;

	for(j in 1:dimenSet[1]){
		if(pred.cross.prob[j,2] > alphaCrit[i]){
			if(Crossval.label[j] == +1){
				countTP = countTP + 1;
			}else{
				countFP  = countFP  + 1;
			} 
		}else{
			if(Crossval.label[j] == -1){
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

return(list(falseRate, hitRate, precisionRate));

}