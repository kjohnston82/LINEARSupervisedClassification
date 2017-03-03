#####################################################################
##### run this script to perform AL from start to finish on ASAS data
## and save all relevant output, classifications + performance metrics
####################################################
# load needed packages and files
set.seed(1)

setwd("C:/Users/kyle.johnston/Desktop/Dissertation/starvars/")

source("src/astro/utils/utils_classify.R")

##########################
# http://cran.r-project.org/web/packages/e1071/e1071.pdf
##########################
library(e1071)
library(RSNNS)
library(FactoMineR)
library(kernlab)
####################################################
####################################################
# load data

library(foreign)

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

####################################################
## Sub-Function Calls

source("src/astro/utils/GenerateTrainingAndCrossval_OC.R")
source("src/astro/utils/GenerateLabels.R")

####################################################
#######################################
# Self Test
sequenceNu = seq(from = -5, to = -0.5, by = 0.1);
sequenceNu = 10^sequenceNu;
lengthArray <- array(dim=c(1, length(sequenceNu )));

for(i in 1:length(sequenceNu)){
svmModel<-svm(trainingData ,y=trainingLabel ,
               type='one-classification',
               nu=sequenceNu[i],
               scale=TRUE,
               kernel="radial")

probTemp  <- predict(svmModel, 
		newdata = trainingData , probability = TRUE)

lengthArray[i] <- length(trainingLabel [!probTemp])/length(trainingLabel )

}
plot(sequenceNu, lengthArray, log ='x', ylab = "Fraction of Anomalous Points", xlab = "Nu")
grid()

####################################################
## Cross-val
lengthArray <- array(dim=c(50, length(sequenceNu)));

for(k in 1:50){

split<- runif(1, 0.4, 0.7)

listOC = GenerateTrainingAndCrossval_OC(feat.train, class.train, split);

trainingData <- listOC[[1]];
crossvalData <- listOC[[2]];
trainingLabel <- listOC[[3]];
crossvalLabel <- listOC[[4]];

####################################################

for(i in 1:length(sequenceNu)){
	svmModel<-svm(trainingData ,y=trainingLabel ,
               type='one-classification',
               nu=sequenceNu[i],
               scale=TRUE,
               kernel="radial")

	probTemp  <- predict(svmModel, 
		newdata = crossvalData , probability = TRUE)

	lengthArray[k,i] =  length(crossvalLabel[!probTemp])/length(crossvalLabel);
}

}

avg = colMeans(lengthArray);
sdev = apply(lengthArray,2,sd);

plot(sequenceNu, avg,
	ylim=range(c(avg-sdev, avg+sdev)),
	log = "x", pch=19, xlab="Nu", ylab="Fraction Anomolous Mean +/- SD",
	main="OC-SVM Applied to Cross-Val Data"
)
# hack: we draw arrows but with very special "arrowheads"
arrows(sequenceNu, avg-sdev, sequenceNu, avg+sdev, length=0.05, angle=90, code=3)

grid()


#############################
####################################################
# Run PCA and see how things break out

svmModel<-svm(trainingData ,y=trainingLabel ,
               type='one-classification',
               nu=0.001,
               scale=TRUE,
               kernel="radial")

	probTemp  <- predict(svmModel, 
		newdata = trainingData , probability = TRUE)


pca3Train = PCA(trainingData , graph = FALSE)
pca3CrossVald = PCA(trainingData , graph = FALSE)

plot(pca3Train $ind$coord[,1], 
pca3Train $ind$coord[,2],
ylab = "PCA 2", xlab = "PCA 1", pch=16)

points(pca3CrossVald$ind$coord[probTemp  ,1],
pca3CrossVald$ind$coord[probTemp ,2], pch=1, col = 2)

points(pca3CrossVald$ind$coord[!probTemp  ,1],
pca3CrossVald$ind$coord[!probTemp ,2], pch=1, col = 3)

legend("bottomright", legend = c("OK", "Anom"), pch = c(1,1), col = c(2,3), cex = 0.8)
grid()

