#####################################################################
##### run this script to perform AL from start to finish on ASAS data
## and save all relevant output, classifications + performance metrics
####################################################
# load needed packages and files
set.seed(1)

setwd("C:/Users/kyle.johnston/Desktop/Dissertation/starvars/")

source("src/astro/utils/utils_classify.R")

##########################
# http://gastonsanchez.com/blog/how-to/2012/06/17/PCA-in-R.html
##########################
library(FactoMineR)
library(kernlab)
library(superpc)
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

source("src/astro/utils/GenerateTrainingAndCrossval.R")
source("src/astro/utils/GenerateLabels.R")

#####################################################################

listLabels = GenerateLabels(class.train);

pulsatingLabels <- listLabels[[1]];
eruptiveLabels <- listLabels[[2]];
multiStarLabels <- listLabels[[3]];
unknownClassLabels <- listLabels[[4]];

####################################################
# Run PCA and see how things break out
pca3 = PCA(feat.train, graph = FALSE)

####################################################
plot(pca3$ind$coord[1:length(eruptiveLabels),1], 
pca3$ind$coord[1:length(eruptiveLabels),2],
ylab = "PCA 2", xlab = "PCA 1", pch=16)
points(pca3$ind$coord[pulsatingLabels== "True",1],pca3$ind$coord[pulsatingLabels== "True",2], pch=1, col = 2)
points(pca3$ind$coord[multiStarLabels == "True",1],pca3$ind$coord[multiStarLabels == "True",2], pch=1, col = 3)
points(pca3$ind$coord[unknownClassLabels== "True",1],pca3$ind$coord[unknownClassLabels== "True",2], pch=1, col = 4)
points(pca3$ind$coord[eruptiveLabels== "True",1],pca3$ind$coord[eruptiveLabels== "True",2], pch=1, col = 5)
legend("bottomright", legend = c("pulsating","multiStar", "other",
"eruptive"), pch = c(1,1,1,1), col = c(2,3,4,5), cex = 0.8)
grid()


plot(pca3$ind$coord[1:length(eruptiveLabels),1], 
pca3$ind$coord[1:length(eruptiveLabels),2],
ylab = "PCA 2", xlab = "PCA 1", pch=16)
points(pca3$ind$coord[pulsatingLabels== "True",1],pca3$ind$coord[pulsatingLabels== "True",2], pch=1, col = 2)
legend("bottomright", legend = c("pulsating"), pch = c(1), col = c(2), cex = 0.8)
grid()

plot(pca3$ind$coord[1:length(eruptiveLabels),1], 
pca3$ind$coord[1:length(eruptiveLabels),2],
ylab = "PCA 2", xlab = "PCA 1", pch=16)
points(pca3$ind$coord[multiStarLabels == "True",1],pca3$ind$coord[multiStarLabels == "True",2], pch=1, col = 3)
legend("bottomright", legend = c("multiStar"), pch = c(1), col = c(3), cex = 0.8)
grid()

plot(pca3$ind$coord[1:length(eruptiveLabels),1], 
pca3$ind$coord[1:length(eruptiveLabels),2],
ylab = "PCA 2", xlab = "PCA 1", pch=16)
points(pca3$ind$coord[unknownClassLabels== "True",1],pca3$ind$coord[unknownClassLabels== "True",2], pch=1, col = 4)
legend("bottomright", legend = c("other"), pch = c(1), col = c(4), cex = 0.8)
grid()

plot(pca3$ind$coord[1:length(eruptiveLabels),1], 
pca3$ind$coord[1:length(eruptiveLabels),2],
ylab = "PCA 2", xlab = "PCA 1", pch=16)
points(pca3$ind$coord[eruptiveLabels== "True",1],pca3$ind$coord[eruptiveLabels== "True",2], pch=1, col = 5)
legend("bottomright", legend = c("eruptive"), pch = c(1), col = c(5), cex = 0.8)
grid()
####################################################

data(iris)
kpc <- kpca(~.,data=iris,kernel="rbfdot",kpar=list(sigma=0.2),features=2)
plot(pcv(kpc)[1:150,1], pcv(kpc)[1:150,2], pch=16)
points(pcv(kpc)[iris$Species == "setosa",1],pcv(kpc)[iris$Species == "setosa",2], pch=1, col = 2)
points(pcv(kpc)[iris$Species == "virginica",1],pcv(kpc)[iris$Species == "virginica",2], pch=1, col = 3)
points(pcv(kpc)[iris$Species == "versicolor",1],pcv(kpc)[iris$Species == "versicolor",2], pch=1, col = 4)


kpc <- kpca(~.,data=feat.train,kernel="laplacedot",kpar=list(sigma=0.001),features=2)
plot(pcv(kpc)[1:length(eruptiveLabels),1], 
pcv(kpc)[1:length(eruptiveLabels),2], 
ylab = "KPCA 2", xlab = "KPCA 1",xlim = c(-0.1, 0.24), pch=16)
points(pcv(kpc)[pulsatingLabels== "True",1],pcv(kpc)[pulsatingLabels== "True",2], pch=1, col = 2)
points(pcv(kpc)[multiStarLabels == "True",1],pcv(kpc)[multiStarLabels == "True",2], pch=1, col = 3)
points(pcv(kpc)[unknownClassLabels== "True",1],pcv(kpc)[unknownClassLabels== "True",2], pch=1, col = 4)
points(pcv(kpc)[eruptiveLabels== "True",1],pcv(kpc)[eruptiveLabels== "True",2], pch=1, col = 5)
legend("bottomright", legend = c("pulsating","multiStar", "other",
"eruptive"), pch = c(1,1,1,1), col = c(2,3,4,5), cex = 0.8)
grid()



####################################################
####################################################
### Supervised PCA
####################################################

## Testing
x<-matrix(rnorm(1000*100),ncol=100)
v1<- svd(x[1:80,])$v[,1]

y<-2+5*v1+ .05*rnorm(100)

xtest<-x
ytest<-2+5*v1+ .05*rnorm(100)
censoring.status<- sample(c(rep(1,80),rep(0,20)))
censoring.status.test<- sample(c(rep(1,80),rep(0,20)))
featurenames <- paste("feature",as.character(1:1000),sep="")

data<-list(x=x,y=y, censoring.status= censoring.status, 
featurenames=featurenames)
data.test<-list(x=xtest,y=ytest, censoring.status=censoring.status.test, 
featurenames= featurenames)

train.obj<- superpc.train(data, type="survival")

## For Real
mixedLabels <- array(dim=c(1,length(class.train)))

mixedLabels[pulsatingLabels == +1] <- 1;
mixedLabels[multiStarLabels == +1] <- 2;
mixedLabels[unknownClassLabels == +1] <- 3;
mixedLabels[eruptiveLabels == +1] <- 4;

#source("src/astro/utils/GenerateTrainingAndCrossval_OC.R")
#listTrainingCrossval <- GenerateTrainingAndCrossval_OC(
#	feat.train, mixedLabels, 0.5);

#trainingData<- listLabels[[1]];
#crossvalData <- listLabels[[2]];
#trainingLabel <- listLabels[[3]];
#crossvalLabel<- listLabels[[4]];


x <- as.matrix(feat.train)
y <- as.vector(mixedLabels)

featurenames <- colnames(feat.train);
data <- list(x=t(x), y=y, featurenames= featurenames);
trainObj<- superpc.train(data, type="regression")

dataOut <- matrix(ncol = 72, nrow = length(pulsatingLabels))
dataOut [,1:68] <- x;
dataOut [,69] <- pulsatingLabels;
dataOut [,70] <- multiStarLabels;
dataOut [,71] <- unknownClassLabels;
dataOut [,72] <- eruptiveLabels;

write.csv(dataOut, "dataOut.csv", row.names=FALSE)


plot(1:length(featurenames), trainObj$feature.scores, 
	xlab = "Features", ylab = "SPCA - Score")

x11()
cvObj<-superpc.cv(trainObj, data)
superpc.plotcv(cvObj)


fitCts<- superpc.predict(trainObj, data, data, threshold=6.0, n.components=3, prediction.type="continuous")

plot(fitCts$v.pred[1:length(eruptiveLabels),1], 
fitCts$v.pred[1:length(eruptiveLabels),2], 
ylab = "KPCA 2", xlab = "KPCA 1", pch=16)
points(fitCts$v.pred[pulsatingLabels== +1,1],fitCts$v.pred[pulsatingLabels== +1,2], pch=1, col = 2)
points(fitCts$v.pred[multiStarLabels == +1,1],fitCts$v.pred[multiStarLabels == +1,2], pch=1, col = 3)
points(fitCts$v.pred[unknownClassLabels== +1,1],fitCts$v.pred[unknownClassLabels== +1,2], pch=1, col = 4)
points(fitCts$v.pred[eruptiveLabels== +1,1],fitCts$v.pred[eruptiveLabels== +1,2], pch=1, col = 5)
legend("bottomright", legend = c("pulsating","multiStar", "other",
"eruptive"), pch = c(1,1,1,1), col = c(2,3,4,5), cex = 0.8)
grid()

plot(fitCts$v.pred[1:length(eruptiveLabels),1], 
fitCts$v.pred[1:length(eruptiveLabels),3], 
ylab = "KPCA 2", xlab = "KPCA 1", pch=16)
points(fitCts$v.pred[pulsatingLabels== +1,1],fitCts$v.pred[pulsatingLabels== +1,3], pch=1, col = 2)
points(fitCts$v.pred[multiStarLabels == +1,1],fitCts$v.pred[multiStarLabels == +1,3], pch=1, col = 3)
points(fitCts$v.pred[unknownClassLabels== +1,1],fitCts$v.pred[unknownClassLabels== +1,3], pch=1, col = 4)
points(fitCts$v.pred[eruptiveLabels== +1,1],fitCts$v.pred[eruptiveLabels== +1,3], pch=1, col = 5)
legend("bottomright", legend = c("pulsating","multiStar", "other",
"eruptive"), pch = c(1,1,1,1), col = c(2,3,4,5), cex = 0.8)
grid()

dataSPCA <- matrix(ncol = 7, nrow = length(pulsatingLabels))
dataSPCA[,1:3] <- fitCts$v.pred[,1:3];
dataSPCA[,4] <- pulsatingLabels;
dataSPCA[,5] <- multiStarLabels;
dataSPCA[,6] <- unknownClassLabels;
dataSPCA[,7] <- eruptiveLabels;

write.csv(dataSPCA , "dataSPCA .csv", row.names=FALSE)