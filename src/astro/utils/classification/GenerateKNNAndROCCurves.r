GeneratekNNAndROCCurves <- function(
	Training.data, Training.label, Crossval.data, Crossval.label,
	summaryType) {
  
#######################################################
setwd("C:/Users/kyle.johnston/Desktop/Dissertation/starvars/")

source("src/astro/classifiers/NaiveKNearest.R")

#######################################################

polySet <- c(1,2,3);
intKNArray <- c(1,2,3);

alphaCrit <- sort(seq(from = -0.01, to = 1.01, by = 0.01), decreasing = TRUE);

hitRate  		= array(dim=c(length(polySet ),length(alphaCrit)));
falseRate 		= array(dim=c(length(polySet ),length(alphaCrit)));
precisionRate 	= array(dim=c(length(polySet ),length(alphaCrit)));


for(k in 1:length(polySet )){

	listKNN <- Naive_K_Nearest(Training.data, Training.label, 
		Crossval.data, 10, polySet[k], 1.0,  "Miss");

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

		hitRate[k,i]   		= countTP /(countTP + countFN);
		falseRate[k,i] 		= countFP /(countFP + countTN);
		precisionRate[k,i] 	= countTP /(countTP + countFP);
	}
}

xdivide <- seq(from = 0, to = 1, by = 0.01);

precisionRate[is.nan(precisionRate)] = 0.9999999

frame()
par(mfrow = c(1,2))
plot(falseRate[1,1:length(alphaCrit)], hitRate[1,1:length(alphaCrit)],
	xlab = "False Rate", ylab = "Hit Rate",
	xlim = c(0.0,1.0), ylim = c(0.0,1.0),
	col=2, lty = 1, type = "l")
lines(falseRate[2,1:length(alphaCrit)], hitRate[2,1:length(alphaCrit)], col=3, lty = 3)
lines(falseRate[3,1:length(alphaCrit)], hitRate[3,1:length(alphaCrit)], col=4, lty = 4)
lines(xdivide , xdivide , col=1, lty = 2)
legend("bottomright", legend = polySet , col = c(2,3,4), lty = c(1,3,4), cex = 0.8, title = "poly-n")
grid()

plot(hitRate[1,1:length(alphaCrit)], precisionRate[1,1:length(alphaCrit)],
	xlab = "Recall", ylab = "Precision", 
	xlim = c(0.0,1.0), ylim = c(0.0,1.0),
	col=2, lty = 2, type = "l")
lines(hitRate[2,1:length(alphaCrit)], precisionRate[2,1:length(alphaCrit)], col=3, lty = 3)
lines(hitRate[3,1:length(alphaCrit)], precisionRate[3,1:length(alphaCrit)], col=4, lty = 4)
legend("bottomright", legend = polySet , col = c(2,3,4), lty = c(1,3,4), cex = 0.8, title = "poly-n")
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

cat("Poly 1: ROC AUC: \n", ROC.AUC,"\n")
cat("Poly 1: PR AUC: \n", PR.AUC,"\n")

ROC.AUC = trapz(falseRate[2,1:length(alphaCrit)], hitRate[2,1:length(alphaCrit)]);
PR.AUC = trapz(hitRate[2,1:length(alphaCrit)],precisionRate[2,1:length(alphaCrit)]);

ROC.AUC = abs(ROC.AUC)
PR.AUC = abs(PR.AUC)

cat("Poly 2: ROC AUC: \n", ROC.AUC,"\n")
cat("Poly 1: PR AUC: \n", PR.AUC,"\n")

return(list(falseRate, hitRate, precisionRate));

}