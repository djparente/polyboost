# Daniel J. Parente, MD PhD
# University of Kansas Medical Center
#
# Training functions for XGBoost
# XGBoost has too many hyperparameters to run grid search over the whole search space
# so instead approximate this by coordinate descent

# Requires that humvar-combined-na.txt and humdiv-combined-na.txt are in the working directory

# This code is heavily adapted from the machine learning examples found in the 
# supplement to Taylor et al. Predicting urinary tract infections in the emergency
# department with machine learning. PLoS One. 2018 Mar 7;13(3):e0194085

# Load required libraries
library(caret)
library(doParallel)
library(dplyr)

# Set up parallelism
cl <- makePSOCKcluster(detectCores()-1) # number of cores to use
registerDoParallel(cl)

xgBoostGradDescent <- function(training) {
	# Declare the prediction formula
	fx <- truth ~ Score1 + dScore + Nobs + dVol + PFamHitBool + IdPmax + IdQmin + CpG + NormASA + Bfact + dProp

	# Perform 10-fold cross validation separating into two classes
	tc <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary, classProbs = TRUE, allowParallel=T)

	# Set random seed so separate runs should be identical
	set.seed(1138)

	# Declare the initial (abbreviated) search grid
	initialGrid <- expand.grid(
	  nrounds=c(1000),
	  max_depth = c(3:10),  
	  eta = c(0.05, 0.1, 0.2), 
	  gamma = c(0.01),
	  colsample_bytree = c(0.75),
	  subsample = c(0.50),
	  min_child_weight = c(0)
	)

	# Perform the initial hyperparametric tuning of XGBoost
	initial_tune <- train(fx,
					 data=training,
					 method="xgbTree",
					 trControl=tc,
					 tuneGrid=initialGrid,
					 verbose=T,
					 metric="ROC",
					 nthread =3,
					 na.action=na.pass
	)

	# Save the initial ROC result
	initialROC <- max(initial_tune$results$ROC)

	# Initialize the best ROC and bestTune parameters from the initial grid search tuning
	bestROC <- initialROC
	bestTune <- initial_tune$bestTune

	abortCounter = 0 # Keeps track of rounds completed without change in parameters

	# Perform coordinate descent
	for(dRound in 1:10) {
		for(i in 1:6) {

			# Extract the best parameters
			Nmax_depth <- bestTune$max_depth
			Neta <- bestTune$eta
			Ncolsample_bytree <- bestTune$colsample_bytree
			Ngamma <- bestTune$gamma
			Nsubsample <- bestTune$subsample
			Nmin_child_weight <- bestTune$min_child_weight
				
			# Give the selected parameter this round a range of values to optimize
			if (i==1) {
				Nmax_depth <- c(3:10)
			}

			if(i == 2) {
				Neta <- seq(0.05, 0.2, 0.05)
			}
			
			if(i == 3) {
				Ngamma <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
			}
			
			if(i == 4) {
				Ncolsample_bytree <- seq(0.5, 1, 0.05)
			}
			
			if(i == 5) {
				Nsubsample <- seq(0.5, 1, 0.05)
			}
			
			if(i == 6) {
				Nmin_child_weight <- seq(0, 1, .1)
			}
			
			# Build the coordinate descent grid
			descentGrid <- expand.grid(
				nrounds = c(1000),
				max_depth = Nmax_depth,  
				eta = Neta,
				gamma = Ngamma,
				colsample_bytree = Ncolsample_bytree,
				subsample = Nsubsample,
				min_child_weight = Nmin_child_weight
			)
			
			# Status report
			print(paste('Running tune',toString(i)))
			
			# Stop tuning after we've gone a full round of coordinate descent without change
			if( abortCounter < 7 ) {
				# Tune the selected coordinate
				descentTune <- train(fx,
					data=training,
					method="xgbTree",
					trControl=tc,
					tuneGrid=descentGrid,
					verbose=T,
					metric="ROC",
					nthread =3,
					na.action=na.pass
				)
			}
			
			# If the trained model achieves a higher ROC than previously seen, save this as the new best tune
			curROC <- max(descentTune$results$ROC)
			if( curROC > bestROC ) {
				bestTune <- descentTune$bestTune
				bestROC <- max(descentTune$results$ROC)
				abortCounter <- 0
			} else {
				abortCounter <- abortCounter + 1
			}
				
			# Status report
			print(paste("ROC: ",toString(bestROC)))
			print(bestTune)
			print(paste("Abort Counter: ",toString(abortCounter)))
		}
	}

	# Extract the best tuning parameters
	Nmax_depth <- bestTune$max_depth
	Neta <- bestTune$eta
	Ncolsample_bytree <- bestTune$colsample_bytree
	Ngamma <- bestTune$gamma
	Nsubsample <- bestTune$subsample
	Nmin_child_weight <- bestTune$min_child_weight

	# Run a final tuning run using the best parameters
	finalGrid <- expand.grid(
		nrounds = c(1000),
		max_depth = Nmax_depth,  
		eta = Neta,
		gamma = Ngamma,
		colsample_bytree = Ncolsample_bytree,
		subsample = Nsubsample,
		min_child_weight = Nmin_child_weight
	)

	final_tune <- train(fx,
					 data=training,
					 method="xgbTree",
					 trControl=tc,
					 tuneGrid=finalGrid,
					 verbose=T,
					 metric="ROC",
					 nthread =3,
					 na.action=na.pass
	)
	
	return(final_tune)
}

# Read in the data
humvar.full <- read.csv('humvar-combined-na.txt', header=T, sep='\t')
humvar.train <- filter(humvar.full, split == "training")

humdiv.full <- read.csv('humdiv-combined-na.txt', header=T, sep='\t')
humdiv.train <- filter(humdiv.full, split == "training")

# Actually run the tuning
tunes.polyboost.hv <- xgBoostGradDescent(humvar.train)
tunes.polyboost.hd <- xgBoostGradDescent(humdiv.train)

# Clean up the workspace
rm(humvar.full)
rm(humvar.train)
rm(humdiv.full)
rm(humdiv.train)
