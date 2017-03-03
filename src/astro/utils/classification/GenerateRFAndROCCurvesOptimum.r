GenerateRFAndROCCurvesOptimum <- function(
	Training.data, Training.label, Crossval.data, Crossval.label,
	summaryType) {
  
alphaCrit <- sort(seq(from = -0.01, to = 1.01, by = 0.01), decreasing = TRUE);

hitRate  		= array(dim=c(1,length(alphaCrit)));
falseRate 		= array(dim=c(1,length(alphaCrit)));
precisionRate 	= array(dim=c(1,length(alphaCrit)));



rf.tr = randomForest(
	x = as.matrix(Training.data),
	y = factor(Training.label),
	mtry = 8, ntree = 200)

pred.cross.prob = predict(rf.tr,
	newdata = Crossval.data,type='prob')
	
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

return(list(falseRate, hitRate, precisionRate, rf.tr));

}