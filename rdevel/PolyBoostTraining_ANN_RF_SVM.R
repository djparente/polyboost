# Daniel J. Parente, MD PhD
# University of Kansas Medical Center
#
# Training functions for RF, SVM, and ANN
# XGBoost has too many variables to run hyperparameter grid search
# so instead approximate this by coordinate descent (see separate file)

# Requires that humvar-combined-na.txt and humdiv-combined-na.txt are in the working directory

# This code is patterned after the machine learning examples found in the supplement to
# Taylor et al. Predicting urinary tract infections in the emergency department with 
# machine learning. PLoS One. 2018 Mar 7;13(3):e0194085

# Import the required libraries
require(caret)
library(RANN)
library(dplyr)

# Read in the data
humvar.full <- read.csv('humvar-combined-na.txt', header=T, sep='\t')
humvar.train <- filter(humvar.full, split == "training")

humdiv.full <- read.csv('humdiv-combined-na.txt', header=T, sep='\t')
humdiv.train <- filter(humdiv.full, split == "training")

# Declare the model function
fx <- truth ~ Score1 + dScore + Nobs + dVol + PFamHitBool + IdPmax + IdQmin + CpG + NormASA + Bfact + dProp

# Train artificial neural network
trainANN <- function(dataset) {
	# 10-fold cross validation with binary classification
	tc <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary, classProbs = T)
		
	# Declare the hyperparameter grid for ANN
	hyperpara.grid = expand.grid( size = seq(20, 100, 20), decay = seq(0.1, 0.8, .1) )

	# Set a consistent seed, so that if we run it multiple times all the random variables will be the same with each fun
	set.seed(45)
	
	# Perform training
	ann.tune <- train(
		fx,
		data=dataset,
		method="nnet",
		trControl=tc,
		tuneGrid=hyperpara.grid,
		metric="ROC",
		MaxNWts = 10000,
		preProcess="medianImpute",
		na.action=na.pass 
	)
	
	# Return the result
	return(ann.tune)
}

# Train random forest
trainRF <- function(dataset) {
	# 10-fold cross validation with binary classification
	tc <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary, classProbs = T)
		
	# Declare the hyperparameter grid for RF
	hyperpara.grid = expand.grid(.mtry=c(10:20))
	
	# Set consistent seed
	set.seed(45)
	
	# Run random forests
	rf.tune <-train(
		fx,
		data=dataset,
		method="rf",
		trControl=tc,
		tuneGrid=hyperpara.grid,
		verbose=T,
		metric="ROC",
		nTree = 100, #TODO: Set 1000
		preProcess="medianImpute",
		na.action=na.pass 
	)
	
	# Return the result
	return(rf.tune)
}

# Train SVM
trainSVM <- function(dataset) {
	# 10-fold cross validation with binary classification
	tc <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary, classProbs = T)
	
	# Set consistent seed
	set.seed(45)
	
	# Run the tune
	svm.tune <-train(
		fx,
		data=dataset,
		method="svmLinear",
		trControl=tc,
		tuneLength = 5,	# Tune specification
		metric="ROC",
		preProcess= c("medianImpute", "center", "scale"),
		na.action=na.pass 
	)
	
	return(svm.tune)
}

# Train ANN on HumVar and HumDiv
#tunes.ann.hv <- trainANN(humvar.train)
#tunes.ann.hd <- trainANN(humdiv.train)

# Train RF on HumVar and HumDiv
#tunes.rf.hv <- trainRF(humvar.train)
#tunes.rf.hd <- trainRF(humdiv.train)

# Train SVM on HumVar and HumDiv
#tunes.svm.hv <- trainSVM(humvar.train)
#tunes.svm.hd <- trainSVM(humdiv.train)

# Cleanup the workspace
#rm(humvar.full)
#rm(humvar.train)
#rm(humdiv.full)
#rm(humdiv.train)
#rm(fx)
