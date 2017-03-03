#####################################################################
##### run this script to perform AL from start to finish on ASAS data
## and save all relevant output, classifications + performance metrics
######################################################################

####################################################
# load needed packages and files

set.seed(1)

setwd("C:/Users/kyle.johnston/Desktop/Dissertation/starvars/")

source("src/astro/utils/utils_classify.R")

####################################################
####################################################
# load data

library(foreign)

####################################################
# Load TRAINING data
cat("Loading ASAS + Hipp + OGLE Training data\n")
feat.train = read.table("data/trainingsets/training_set_features.dat",header=TRUE,sep=',')
ID.train = feat.train$ID.train
exclude = c("ID.train","qso_log_chi2_qsonu","qso_log_chi2nuNULL_chi2nu","freq_n_alias")
feat.train = feat.train[,-which(names(feat.train) %in% exclude)]
ra.dec.train = read.table("data/trainingsets/training_set_ra_dec.dat",header=TRUE,sep=',')
class.train = read.table("data/trainingsets/training_set_id_class.dat",header=FALSE,sep='\t')[,2]
n = length(ID.train)

cat(paste("Training data loaded: ",n," sources\n"))


####################################################
# load LINEAR data
cat("Loading LINEAR data\n")

lineardat=read.arff(file="data/newdata/200k_combo.arff")
#lineardat = read.table("data/newdata/200k_final_combo.csv",sep=",",header=TRUE)
lineardat$freq_rrd = rep(0,dim(lineardat)[1])
feat.tmp = data.frame(lineardat)
feat.linear = data.frame(matrix(0,dim(feat.tmp)[1],dim(feat.train)[2]))
  for(ii in 1:dim(feat.train)[2]){
    feat.linear[,ii]= feat.tmp[,names(feat.tmp)==names(feat.train)[ii]]
  }
colnames(feat.linear) = names(feat.train)
feat.linear[is.na(feat.linear)] = 0

# fix the ? problem
fixfeat = function(feature){
  feature = paste(feature)
  feature[feature=="?"] = 0
  feature = as.numeric(feature)
  return(feature)
}

feat.linear$amplitude = fixfeat(feat.linear$amplitude)
feat.linear$freq_frequency_ratio_21 = fixfeat(feat.linear$freq_frequency_ratio_21)
feat.linear$freq_frequency_ratio_31 = fixfeat(feat.linear$freq_frequency_ratio_31)
feat.linear$freq_signif_ratio_21 = fixfeat(feat.linear$freq_signif_ratio_21)
feat.linear$freq_signif_ratio_31 = fixfeat(feat.linear$freq_signif_ratio_31)
feat.linear$freq1_lambda = fixfeat(feat.linear$freq1_lambda)
feat.linear$max_slope = fixfeat(feat.linear$max_slope)
feat.linear$median_absolute_deviation = fixfeat(feat.linear$median_absolute_deviation)
feat.linear$p2p_scatter_2praw = fixfeat(feat.linear$p2p_scatter_2praw)
feat.linear$p2p_scatter_over_mad = fixfeat(feat.linear$p2p_scatter_over_mad)
feat.linear$p2p_scatter_pfold_over_mad = fixfeat(feat.linear$p2p_scatter_pfold_over_mad)
feat.linear[is.na(feat.linear)] = 0

###################################################
ID.linear = lineardat$source_id
N = dim(feat.linear)[1]

## add RRLd feature
feat.linear$freq_rrd = ifelse(abs(feat.linear$freq_frequency_ratio_21 - 0.746) < 0.0035 | abs(feat.linear$freq_frequency_ratio_31 - 0.746) < 0.0035, 1,0)

##############################################
# read in mastermain file (ra/dec plus colors)
ra.dec.linear = read.table("data/newdata/masterMain.dat.txt",header=TRUE)

## Pick the ones that are similar between the two
infeature = which(ra.dec.linear$objectID %in% ID.linear)
raDecLINEAR <- ra.dec.linear[infeature,]

## Common ID between the two arrays (
commonID = raDecLINEAR$objectID;

switch <- logical(length(commonID));
nonUniqueLogical <- logical(length(ID.linear))

for(i in 1:length(ID.linear)){

	# which in common id is it?
	indexMatch = which(ID.linear[i] == commonID)

	if(length(indexMatch) == 0){
		next	
	}

	# have I seen it before?
	if(!switch[indexMatch]){

		nonUniqueLogical[i] = TRUE;
		switch[indexMatch] = TRUE;
	}
}

featLINEAR <- feat.linear[nonUniqueLogical,]

IDLinearReduced <- ID.linear[nonUniqueLogical]length(

# align the mastermain table to the feature table
orderfeat = order(IDLinearReduced)
ordermain = order(raDecLINEAR$objectID)

# sort features and main table by object ID
featLINEAR  = featLINEAR[orderfeat,]
IDLinearReduced = IDLinearReduced[orderfeat]
raDecLINEAR  = raDecLINEAR [ordermain,]


##############################################
# compute colors for the LINEAR data
source("src/astro/utils/ColorTransformUgriz2BVR.R")

linear.bvr = ColorTransformUgriz2BVR(
	raDecLINEAR$uMod,
	raDecLINEAR$gMod, raDecLINEAR$rMod,
	raDecLINEAR$iMod, raDecLINEAR$zMod)

featLINEAR$color_diff_bj = linear.bvr$B - raDecLINEAR$J
featLINEAR$color_diff_rj = linear.bvr$R - raDecLINEAR$J
featLINEAR$color_diff_vj = linear.bvr$V - raDecLINEAR$J


cat(paste("LINEAR data loaded: ",n," sources\n"))

save(featLINEAR, IDLinearReduced ,raDecLINEAR ,
	 file = "data/LINEARDataset.R")

####################################################
# load LINEAR labels on the 7k
## cat("Loading LINEAR Labels for the 7k data set\n")
## linear.eb = read.table("data/EBs_IDs.txt",header=FALSE)[,1]
## linear.rrc = read.table("data/data_list_RRc.txt",header=TRUE)[,1]
## linear.rrab = read.table("data/RRab.dat",header=FALSE)[,1]
## linear.sxe = read.table("data/2013_02_04_SXPhe.dat",header=TRUE)[,1]

## linear.class = rbind(cbind(linear.eb, rep("EB",length(linear.eb))),
##   cbind(linear.rrab, rep("RRab",length(linear.rrab))),
##   cbind(linear.rrc, rep("RRc",length(linear.rrc))),
##   cbind(linear.sxe, rep("SXPhe",length(linear.sxe))))

### 2013-06-17 USING THE LINEARattributesFinalApr2013.dat file
# 1 = ab RR Lyr (2612/2923); 2 = c RR Lyr
# (864/990); 4 = Algol-like with 2 minima (342/357); 5 = contact binary (2246/2385);
# 6 = delta Scu/SX Phe (82/112);
linear.class = read.table("data/trainingsets/LINEARattributesFinalApr2013.dat")[,c(14,13)]
linear.class[,2] = factor(linear.class[,2],labels=c("RRab","RRc","Algol","Contact","DScu/SXPhe"))

cat(paste("Number of manually labeled LINEAR objects: ",dim(linear.class)[1]," sources\n"))


####################################################
## add LINEAR - trainingset overlap data to training set
cat("Adding LINEAR training objects that overlap with Existing Training Set (ASAS + Debosscher)\n")
counter = 0
for(ii in 1:dim(ra.dec.train)[1]){
  ind = which(round(ra.dec.linear$ra,2)==round(ra.dec.train[ii,2],2) & round(ra.dec.linear$dec,2)==round(ra.dec.train[ii,3],2))
  if(length(ind)>0){
    ind1 = which(ID.linear == ra.dec.linear$objectID[ind])
    ind2 = which(ID.train == ra.dec.train[ii,1])
    if(length(ind2)==1){
      counter = counter+1
      cl.tmp = paste(class.train[ind2])
      cat("Exchanging ID: ",ID.train[ind2]," Class: ",cl.tmp,"\n")
      ## remove old Training (ASAS / Deb) data
      feat.train = feat.train[-ind2,]
      class.train = class.train[-ind2]
      ID.train = ID.train[-ind2]
      ## add new LINEAR data
      feat.train = rbind(feat.train,feat.linear[ind1,])
      class.train = c(paste(class.train),paste(cl.tmp))
      ID.train = c(ID.train,ID.linear[ind1])
    }
  }
}
class.train = factor(class.train)
nn = length(class.train)
# total size: 
print(table(class.train[(nn-counter+1):nn])[table(class.train[(nn-counter+1):nn])>0])
cat("number added:",counter,"\n")
# 1 RRL FM, LINEAR ID 6196048 / 
# Exchanging ID:  230655  Class:  g. RR Lyrae FM 
featureCombinedTrain <- feat.train;
labelCombinedTrain <- class.train;

save(featureCombinedTrain , labelCombinedTrain ,
	 file = "data/Training_Dataset_AohLINEAR.R")
