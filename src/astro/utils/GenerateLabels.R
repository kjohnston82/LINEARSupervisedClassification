GenerateLabels <- function(
	class.train) {

####################################################
# Reduce classes to "pulsating, eruptive, multi-star, 
# and unknown" in "two class" system

pulsatingLabels <- array(dim=c(1,length(class.train)))
eruptiveLabels <- array(dim=c(1,length(class.train)))
multiStarLabels <- array(dim=c(1,length(class.train)))
unknownClassLabels <- array(dim=c(1,length(class.train)))

for(i in 1:length(class.train)){

	isOther <- +1;

	##########################
	# Eruptive Type Variables 
	if(
		length(grep(class.train[i],"r. Wolf-Rayet "))>0 || 
		length(grep(class.train[i],"q. Chem. Peculiar"))>0 || 
		length(grep(class.train[i],"p. Per. Var. SG"))>0 ||
		length(grep(class.train[i],"t. Herbig AE/BE"))>0 ||
		length(grep(class.train[i],"u. S Doradus"))>0 ||
		length(grep(class.train[i],"s1. Class. T Tauri"))>0){

		eruptiveLabels[i] <- +1;
		isOther <- -1;
	}else{
		eruptiveLabels[i] <- -1;
	}

	############################
	## Pulsating Type Variables
	if(
		length(grep(class.train[i],"a. Mira"))>0 ||
		length(grep(class.train[i],"b1. Semireg PV"))>0 ||
		length(grep(class.train[i],"c. RV Tauri"))>0){
	## Giants
		pulsatingLabels[i] <- +1
		isOther <- -1;
	}else if(
		length(grep(class.train[i],"d. Classical Cepheid"))>0 ||
		length(grep(class.train[i],"e. Pop. II Cepheid"))>0 ||
		length(grep(class.train[i],"f. Multi. Mode Cepheid"))>0){
	## Cepheids
		pulsatingLabels[i] <- +1
		isOther <- -1;
	}else if(
		length(grep(class.train[i],"h. RR Lyrae FO"))>0 ||
		length(grep(class.train[i],"g. RR Lyrae FM"))>0 ||
		length(grep(class.train[i],"i. RR Lyrae DM"))>0){
	## RR Lyrae
		pulsatingLabels[i] <- +1
		isOther <- -1;
	}else if(
		length(grep(class.train[i],"j. Delta Scuti"))>0 ||
		length(grep(class.train[i],"k. Lambda Bootis"))>0 ||
		length(grep(class.train[i],"l. Beta Cephei"))>0 ||
		length(grep(class.train[i],"m. Slowly Puls. B"))>0 ||
		length(grep(class.train[i],"o. Pulsating Be"))>0 ||
		length(grep(class.train[i],"j1. SX Phe"))>0 ||
		length(grep(class.train[i],"n. Gamma Doradus"))>0){
	## "Other"
		pulsatingLabels[i] <- +1
		isOther <- -1;
	}else{
		pulsatingLabels[i] <- -1
	}

	##############
	## Multi-Star
	if(
		length(grep(class.train[i],"v. Ellipsoidal"))>0 ||
		length(grep(class.train[i],"w. Beta Persei"))>0 ||
		length(grep(class.train[i],"x. Beta Lyrae"))>0 ||
		length(grep(class.train[i],"y. W Ursae Maj."))>0){
	## "Multi-Star"
		multiStarLabels [i] <- +1
		isOther <- -1;
	}else{
		multiStarLabels [i] <- -1
	}

	##############
	## Other
	if(isOther == +1){
		unknownClassLabels[i] <- +1;
	}else{
 		unknownClassLabels[i] <- -1;
	}

}

cat("Number of Eruptive Type Variables:",sum(eruptiveLabels == +1),"\n")
cat("Number of Pulsating Type Variables:",sum(pulsatingLabels== +1),"\n")
cat("Number of Multi-Star Variables:",sum(multiStarLabels == +1),"\n")
cat("Number of Other Variables:",sum(unknownClassLabels== +1),"\n")

return(list(pulsatingLabels, eruptiveLabels , multiStarLabels , unknownClassLabels))

}
