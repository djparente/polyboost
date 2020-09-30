# Daniel J. Parente, MD PhD
# University of Kansas Medical Center

# This file analyzes the performance of various classifiers. It assumes that these have already been loaded into
# the workspace named like: tunes.polyboost.hv, tunes.polyboost.hd, tunes.rf.hv, tunes.rf.hd, tunes.svm.hv, tunes.svm.hd, tunes.ann.hv, and tunes.ann.hd for XGBoost (polyboost), random forest (RF), support vector machine (SVM) and artificial neural network (ANN) models, respectively.

# It also assumes that there are files humvar-combined-na.txt and humdiv-combined-na.txt that are pre-processed (preprocess.py) versions of PolyPhen-2 output on the humvar and humdiv files. Finally, it assumes that there are two files clinvar.small.hv.txt and clinvar.small.hd.txt that are pre-processed (preprocess.py) versions of PolyPhen-2 predictions for ClinVar that have had 11-feature vectors that are identical to HumVar and HumDiv features removed.

# We also want to compare our classifier to the performance of PolyPhen-2 on the HumVar and HumDiv training sets. However, since these were TRAINED on HumVar and HumDiv, we can't just run the training sets through PolyPhen-2 and report out the statistics (otherwise it would be validating on the training data, so we will instead compare it to the 5-fold cross validation results published in the supplement to Adzheubi et al.)

# Load all the libraries we will need
require(ROCR)
require(plotROC)
library(psych)
library(dplyr)
require(pROC)
library(Hmisc)
library(ggplot2)
library(cowplot)

# Data from the Adzhubei et al. supplement on PolyPhen-2 HumVar and HumDiv performance during cross-validation
pph.fpr <- seq(from=0.1, to=0.8, by=.05)
pph.tpr.hv <- c(0.552,0.659,0.734,0.787,0.836,0.867,0.897,0.921,0.940,0.954,0.965,0.973,0.979,0.984,0.988)
pph.tpr.hd <- c(0.767,0.866,0.918,0.947,0.963,0.973,0.976,0.978,0.980,0.982,0.984,0.986,0.988,0.990,0.992)

# Load the test datasets
# HumVar
humvar.full <- read.csv('humvar-combined-na.txt', header=T, sep='\t')
humvar.train <- filter(humvar.full, split == "training")
humvar.val <- filter(humvar.full, split == "validation")

#HumDiv
humdiv.full <- read.csv('humdiv-combined-na.txt', header=T, sep='\t')
humdiv.train <- filter(humdiv.full, split == "training")
humdiv.val <- filter(humdiv.full, split == "validation")

# ClinVar
clinvar.hv <- read.csv('clinvar.small.hv.txt', header=T, sep='\t')
clinvar.hd <- read.csv('clinvar.small.hd.txt', header=T, sep='\t')

# Calculated parameters from the sensitivity, specificity and prevalence
deriveParameters <- function(sen, spec, prev) {
	youden=sen+spec-1	# Youden index/
	ppv = (sen * prev) / ( (sen * prev) + (1-spec) * (1-prev)) # Pos predictive value`
	npv = (spec * (1-prev)) / ((1-sen)*prev + spec * (1-prev)) # Neg predictive value
	lrp = sen / (1-spec) # Positive likelihood ratio
	lrn = (1-sen)/spec # Negative likelihood ratio
	dor = lrp/lrn # Diagnostic odds ratio
	f1 = harmonic.mean(c(ppv, sen)) # F1 score
	
	# Return the parameters in a list
	return(list(youden=youden, ppv=ppv, npv=npv, lrp=lrp, lrn=lrn, dor=dor, f1=f1))
}

# For a model and set of training data, get the coordinates of the best threshold
getTrainThreshold <- function(model, training) {
	# Get deleterious probabilities under the model
	trainPreds <- predict(model, training, type="prob", na.action=na.pass)[,1]
	
	# Construct an ROC curve for these predictions
	trainRoc <- roc(training$truth, trainPreds, levels=c("n", "d"), direction="<")
	
	# Get the best coordinates
	trainCoord <- coords(trainRoc, "best", transpose=F)
	
	# Return the threshold
	return(trainCoord$threshold)
}

# Evaluation subroutine
doEval <- function(model, valid, modelname, trainingname, validname,  preds=NULL) {
	# Obtain deleterious probabilities, if not specified (required for PPH-2 evaluation)
	if(missing('preds')) {
		preds <- predict(model, valid, type="prob", na.action=na.pass)[,1] # Get deleterious
	}
	# Get an ROC curve and the coordinates threshold/sensitivity/specificity at the best threshold
	curRoc <- roc(valid$truth, preds, levels=c("n", "d"), direction="<") 
	bestCoord <- coords(curRoc, "best", transpose=F)
	
	# Assume a pre-test probability ("prevalence") of 50% -- i.e. a situation of clinical equipoise
	prev = 0.5
	
	# Extract the sensitivity, specificity, threshold, youden index, AUC and other derived parameters from the sen/spec/prev (i.e the LR+, LR-, DOR, PPV, NPV, F1 score)
	sen = bestCoord$sensitivity
	spec = bestCoord$specificity
	thresh = bestCoord$threshold
	youden = bestCoords = sen + spec - 1
	auc = curRoc$auc
	dParam <- deriveParameters(sen, spec, prev)
	
	# Build a row of data frame specifying the model, training set, validation set, and all of the performance statistics of interest
	res <- data.frame(model=modelname, train=trainingname, valid=validname, auc=auc, youden=youden, sen=sen, spec=spec, lrp=dParam$lrp, lrn=dParam$lrn, dor=dParam$dor, ppv=dParam$ppv, npv=dParam$npv, f1=dParam$f1, stringsAsFactors=F)
	
	ret <- list(roc=curRoc, best=bestCoord, thresh=thresh, preds=preds, resultframe=res)
	
	return(ret)
}

# Use the evaluation method above to calculate internal validations
intEvaluations <- list(	
	# Internal
	polyboost.hv.i=doEval(tunes.polyboost.hv, humvar.val, "XGBoost", "HumVar", "Internal"),
	polyboost.hd.i=doEval(tunes.polyboost.hd, humdiv.val, "XGBoost", "HumDiv", "Internal"),
	
	rf.hv.i=doEval(tunes.rf.hv, humvar.val, "Random forest", "HumVar", "Internal"),
	rf.hd.i=doEval(tunes.rf.hd, humdiv.val, "Random forest", "HumDiv", "Internal"),
	
	ann.hv.i=doEval(tunes.ann.hv, humvar.val, "Neural network", "HumVar", "Internal"),
	ann.hd.i=doEval(tunes.ann.hd, humdiv.val, "Neural network", "HumDiv", "Internal"),
	
	svm.hv.i=doEval(tunes.svm.hv, humvar.val, "SVM", "HumVar", "Internal"),	
	svm.hd.i=doEval(tunes.svm.hd, humdiv.val, "SVM", "HumDiv", "Internal")
)

# Use the evaluation method above to calculate external validations
extEvaluations <- list(
	# External
	polyphen2.hv.e=doEval(NULL, clinvar.hv, "PolyPhen-2", "HumVar", "External", preds=clinvar.hv$pph2_prob),
	polyphen2.hd.e=doEval(NULL, clinvar.hd, "PolyPhen-2", "HumDiv", "External", preds=clinvar.hd$pph2_prob),
	
	polyboost.hv.e=doEval(tunes.polyboost.hv, clinvar.hv, "XGBoost", "HumVar", "External"),
	polyboost.hd.e=doEval(tunes.polyboost.hd, clinvar.hd, "XGBoost", "HumDiv", "External"),
	
	rf.hv.e=doEval(tunes.rf.hv, clinvar.hv, "Random forest", "HumVar", "External"),
	rf.hd.e=doEval(tunes.rf.hd, clinvar.hd, "Random forest", "HumDiv", "External"),
	
	ann.hv.e=doEval(tunes.ann.hv, clinvar.hv, "Neural network", "HumVar", "External"),
	ann.hd.e=doEval(tunes.ann.hd, clinvar.hd, "Neural network", "HumDiv", "External"),
	
	svm.hv.e=doEval(tunes.svm.hv, clinvar.hv, "SVM", "HumVar", "External"),
	svm.hd.e=doEval(tunes.svm.hd, clinvar.hd, "SVM", "HumDiv", "External")
)

# Accumulate the tables for internal and external validation, and limit output to 3 decimals
accumulateFrameTable <- function(frameList, element) {
	res <- data.frame()
	for(i in frameList) {
		res <- rbind(res, i[[element]]	)
	}	
		
	is.num <- sapply(res, is.numeric)
	res[is.num] <- lapply(res[is.num], round, 3)
	
	return(res)
}

# Internal validation on PolyPhen-2 relies on the published tables in the Adzhubei et al. supplement. Sensitivity and specificity at
# the optimal threshold were determined by a 5th order polynomial fit (and specified in the call to this function, below). We have to manually add PolyPhen-2 to the evaluation table using this sensitivity and specificity. Note that AUC cannot be determined because the tables do not specify the full range of false positive rates (so we only have a partial ROC curve); set AUC to zero as a reminder that this value is not meaningful in the table (and remove from the final data table).
pphFrameTable <- function(sen, spec, modelname, trainingname, validname) {
	# Calculate performance statistics from the sensitivity and specificity at pretest probability of 0.5 (clinical equipoise)
	dParam <- deriveParameters(sen, spec, 0.5)
	
	# Calculate the youden index
	youden <- sen+spec - 1
	
	# Return a row for the table with the performance statistics
	res <- data.frame(model=modelname, train=trainingname, valid=validname, auc=0, youden=youden, sen=sen, spec=spec, lrp=dParam$lrp, lrn=dParam$lrn, dor=dParam$dor, ppv=dParam$ppv, npv=dParam$npv, f1=dParam$f1, stringsAsFactors=F)
	
	return (res)
}

# Accumulate a table summarizing the external validation performance
extdf <- accumulateFrameTable(extEvaluations, 'resultframe')
extdf[order(-extdf$auc),]	# Display the table ordered by AUC in the console

# Accumulate a table summarizing the internal validation performance
intdf <- accumulateFrameTable(intEvaluations, 'resultframe')
# Add the manual PolyPhen-2 calculations to the internal evaluation frame -- sensitivities and specificities are determined by
# 5th order polynomial fit
intdf <- rbind(intdf, pphFrameTable(0.787, 0.750, "PolyPhen-2", "HumVar", "Internal"))
intdf <- rbind(intdf, pphFrameTable(0.900, 0.820, "PolyPhen-2", "HumDiv", "Internal"))
intdf # Display the internal validation table on the console


# Begin code for the external validation figure

# Helper function to extract a dataframe from an ROC curve for plotting purposes
extractROCFrame <- function(roc, name) {
	data.frame(tpr=roc$sensitivities, fpr=1-roc$specificities, method=name)
}

# Function to create the external validation figure
createExternalValidationFigure <- function(extEvaluations) {
	# Extract the data frames form the ROC curves for the external validation figure (XGBoost, RF and PolyPhen-2 all on HumVar, plus PolyPhen-2 on HumDiv)
	extROC <- data.frame()
	extROC <- rbind(extROC, extractROCFrame(extEvaluations$polyboost.hv.e$roc, "XGBoost/HumVar"))
	extROC <- rbind(extROC, extractROCFrame(extEvaluations$rf.hv.e$roc, "Random Forest/HumVar"))
	extROC <- rbind(extROC, extractROCFrame(extEvaluations$polyphen2.hv.e$roc, "PolyPhen-2/HumVar"))
	extROC <- rbind(extROC, extractROCFrame(extEvaluations$polyphen2.hd.e$roc, "PolyPhen-2/HumDiv"))

	# Create the ROC curve figure
	res <- ggplot(extROC, aes(x=fpr, y=tpr, group=method)) +
		geom_line(aes(linetype=method, color=method), size=1) + 
		scale_linetype_manual(values=c("solid", "longdash", "dotted", "dotdash")) + 
		scale_color_manual(values=c("#000000", "#FF0000", "#00FF00", "#0000FF")) + 
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
		theme(legend.position = c(0.7, 0.2)) +
		xlab('False Positive Rate') +
		ylab('True Positive Rate') +
		geom_abline(intercept = 0, slope = 1, size=1.2) +
		labs(color="External: ClinVar", linetype="External: ClinVar")
	
	return(res)
}

ext.val.fig <- createExternalValidationFigure(extEvaluations)

# Internal validation figures
# Figure type is either validation on HumVar or HumDiv; partial amount limits the x and y axis to show the most
# important part of the partial ROC curve (limited space in the main multipanel figure; when partialAmount=0 shows the full ROC curve, which is included in the supplement)
intValidFigure <- function(figType, partialAmount=0) {
	# Extract data for internal validation on HumVar or HumDiv
	intROC <- data.frame()
	if( figType == "HumVar" ) {
		intROC <- rbind(intROC, extractROCFrame(intEvaluations$polyboost.hv.i$roc, "XGBoost"))
		intROC <- rbind(intROC, extractROCFrame(intEvaluations$rf.hv.i$roc, "Random Forest"))
		intROC <- rbind(intROC, extractROCFrame(intEvaluations$svm.hv.i$roc, "SVM"))
		intROC <- rbind(intROC, extractROCFrame(intEvaluations$ann.hv.i$roc, "Neural network"))
		intROC <- rbind(intROC, data.frame(tpr=pph.tpr.hv, fpr=pph.fpr, method="PolyPhen-2 (Pub.)")) # PolyPhen-2 is based on the tabular supplement dat
	} else if (figType == "HumDiv") {
		intROC <- rbind(intROC, extractROCFrame(intEvaluations$polyboost.hd.i$roc, "XGBoost"))
		intROC <- rbind(intROC, extractROCFrame(intEvaluations$rf.hd.i$roc, "Random Forest"))
		intROC <- rbind(intROC, extractROCFrame(intEvaluations$svm.hd.i$roc, "SVM"))
		intROC <- rbind(intROC, extractROCFrame(intEvaluations$ann.hd.i$roc, "Neural network"))
		intROC <- rbind(intROC, data.frame(tpr=pph.tpr.hd, fpr=pph.fpr, method="PolyPhen-2 (Pub.)")) # PolyPhen-2 is based on the tabular supplement dat
	} else {
		stop("Invalid internal validation figure type")
	}
	
	# Set up the legend title to specify the internal validation dataset
	legendtitle <- paste("Internal: ", figType)
	
	# Plot the figure
	res <- ggplot(intROC, aes(x=fpr, y=tpr, group=method)) +
		geom_line(aes(linetype=method, color=method), size=0.5) + 
		geom_point(aes(shape=method, color=method), size=2) + 
		scale_linetype_manual(values=c("solid", "longdash", "dotted", "dotdash", "solid")) + 
		scale_color_manual(values=c("#000000", "#FF0000", "#00FF00", "#0000FF", "333333")) + 
		scale_shape_manual(values=c(NA, NA, NA, NA, 15)) + # Only show data points for PolyPhen-2 (since it is interpolating these tabular data points, whereas the rest of the ROC curves are essentially continuous
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
		theme(legend.position = c(0.7, 0.35)) +
		xlab('False Positive Rate') +
		ylab('True Positive Rate') +
		xlim(0, 1-partialAmount) +
		ylim(partialAmount, 1) + 
		geom_abline(intercept = 0, slope = 1, size=1.2) +
		labs(color=legendtitle, linetype=legendtitle, shape=legendtitle)
	
	# Return the figure
	return(res)
}

# Extract a true positive rate nomogram table -- contemplates the possibility that you might want multiple ROCs represented
extractTPRTable <- function(fprs, rocs, rocnames) {
	ret <- data.frame()

	# Set up the threshold names by adding 'thresh' to the end of the ROC name
	threshnames <- rocnames
	threshnames <- sapply(rocnames, paste0, 'thresh')

	# For each FPR...
	for(fpr in fprs) {
		spec = 1-fpr # Calculate specificity
		newRow <- c(fpr) # Start building a new row
		# Extract the sensitivity from each ROC at the given specificity and add to the row
		for(roc in rocs) {
			newRow <- c(newRow, coords(roc, x=spec, input="specificity", transpose=F)$sensitivity)
		}
		# Extract the threshold from each ROC at the given specificity and add to the row
		for(roc in rocs) {
			# You actually cannot use the threshold as an input to get this is a reasonable way, so instead extract
			# all coordinates and use dplyr to filter them such that the first row is the first time the specificity is above
			# a given threshold; extract the first row ([1,]) and then the threshold ($threshold)
			thresh <- filter(coords(roc, x="all", transpose=F), specificity > spec)[1,]$threshold
			newRow <- c(newRow, thresh) # Add this to the growing row
		}
		# Add this row to the growing table
		ret <- rbind(ret, newRow)
	}
	
	# Set the column names: FPR, <Name for each ROC curve>, <Name for ROC curve + plus 'thresh'>
	colnames(ret) <- cbind(c('FPR', rocnames, threshnames))
	
	# Format the FPR to 2 decimals
	ret$FPR <- sapply(sapply(ret$FPR, round, 2), format, nsmall=2)
	
	# Format the TPRs to 3 digits
	for(rocname in rocnames) {
		ret[[rocname]] <- sapply(sapply(ret[[rocname]], round, 3), format, nsmall=3)
	}
	
	# Format the thresholds to 3 digits
	for(tname in threshnames) {
		ret[[tname]] <- sapply(sapply(ret[[tname]], round, 3), format, nsmall=3)
	}
	
	# Return the data table
	return(data.frame(ret))
}

# Extract a false positive rate nomogram table
extractFPRTable <- function(tprs, rocs, rocnames) {
	# Empty data frame
	ret <- data.frame()

	# Set up the threshold names by adding 'thresh' to the end of the ROC name
	threshnames <- rocnames
	threshnames <- sapply(rocnames, paste0, 'thresh')

	# For each TPR...
	for(tpr in tprs) {
		# Start a new table row and add the TPR
		newRow <- c(tpr)
		# For each ROC, extract the FPR (1-specificity) at the given sensitivity (TPR)
		for(roc in rocs) {
			newRow <- c(newRow, 1-coords(roc, x=tpr, input="sensitivity", transpose=F)$specificity)
		}
		# For each ROC, extract the threshold at a given sensitivity (see comment in similar loop above for
		# why this requires using dplyr instead of using coords directly to query threshold)
		for(roc in rocs) {
			thresh <- filter(coords(roc, x="all", transpose=F), sensitivity <= tpr)[1,]$threshold
			newRow <- c(newRow, thresh)
		}
		# Add this row to the growing table
		ret <- rbind(ret, newRow)
	}
	
	# Set the column names: FPR, <Name for each ROC curve>, <Name for ROC curve + plus 'thresh'>
	colnames(ret) <- cbind(c('TPR', rocnames, threshnames))
	
	# Format the TPR to 2 digits
	ret$TPR <- sapply(sapply(ret$TPR, round, 2), format, nsmall=2)
	
	# Format the FPR to 3 digits
	for(rocname in rocnames) {
		ret[[rocname]] <- sapply(sapply(ret[[rocname]], round, 3), format, nsmall=3)
	}
	
	# Format the threshold to 3 digits
	for(tname in threshnames) {
		ret[[tname]] <- sapply(sapply(ret[[tname]], round, 3), format, nsmall=3)
	}
	
	# Return
	return(data.frame(ret))
}

# We care more about low-ish FPRs and high-ish TPRs, so calculate a table sampling these values more finely (every 1%), and then broaden out to sampling every 5%, up to a maximum of 99%
nomogramRates <- c(seq(0.01, 0.25, 0.01), seq(0.3, 0.95, 0.05), c(.99))

# Extract the tables
extTPR <- extractTPRTable(nomogramRates, list(extEvaluations$polyboost.hv.e$roc), c("PolyBoost"))
extTPR

extFPR <- extractFPRTable(1-nomogramRates, list(extEvaluations$polyboost.hv.e$roc), c("PolyBoost"))
extFPR

# Write some data output
write.table(extTPR, "External TPR Nomogram.txt", sep="\t", quote=F)
write.table(extFPR, "External FPR Nomogram.txt", sep="\t", quote=F)

write.table(extdf, "External Validation Table.txt", sep="\t", quote=F)
write.table(intdf, "Internal Validation Table.txt", sep="\t", quote=F)

# Construct internal validation figures
int.val.fig.hv.partial <- intValidFigure("HumVar", partialAmount=0.3)
int.val.fig.hv.full <- intValidFigure("HumVar", partialAmount=0)
int.val.fig.hd.partial <- intValidFigure("HumDiv", partialAmount=0.3)
int.val.fig.hd.full<- intValidFigure("HumDiv", partialAmount=0)

# Create the combination internal and external validation figure for the main article

# Upper two panels of the combined figure
val.fig.comb.up <- plot_grid(
	int.val.fig.hv.partial,
	int.val.fig.hd.partial,
	labels=c("A", "B"),
	ncol=2,
	nrow=1
)

# Add the external validation (lower panel) to the combined figure
val.fig.combined <- plot_grid(
	val.fig.comb.up,
	ext.val.fig,
	labels=c("", "C"),
	ncol = 1,
	nrow = 2,
	rel_heights=c(1,2)
)

# Display the combined figure
val.fig.combined

# Save a high-resolution summary validation figure and the full internal validation figures for the supplement
ggsave('Combined-Validation-Figure.png', plot=val.fig.combined, dpi=1200, height=9, width=6, unit="in")
ggsave('HumVar-Internal-Full.png', int.val.fig.hv.full, dpi=1200, height=6, width=6, unit="in")
ggsave('HumDiv-Internal-Full.png', int.val.fig.hd.full, dpi=1200, height=6, width=6, unit="in")

# Finally get a table showing data set statistics
getDatasetTableRow <- function(dataset, setname, splitname) {
	# For each dataset, get the num and percent
	counts <- table(dataset$truth)	# Get the counts of deleterious and benign variants
	total <- sum(counts) # Get the sum of all variants
	pct <- counts / total * 100 # Calculate the percentage of benign and deleterious variants
	
	# Format this nicely like, for example: 20,000 (51.2) to indicate 20,000 variants are 51.2% of the total dataset; this is the typesetting format for the data table in the article. Note that counts and pct is a vector showing the deleterious and benign totals at index 1 and index 2, respectively; so it is really doing two formats here
	formatted <- paste0(counts, " (", format(round(pct, digits=1), nsmall=1), ")")
	
	# Return a row with the statistics
	return(data.frame(dataset=setname, split=splitname, total=total, deleterious=formatted[1], neutral=formatted[2]))
	
	
}

# Build the dataset summary table
buildDatasetTable <- function() {
	df <- data.frame()
	df <- rbind(df, getDatasetTableRow(humvar.train, "HumVar", "Training"))
	df <- rbind(df, getDatasetTableRow(humvar.val, "HumVar", "Validation"))
	df <- rbind(df, getDatasetTableRow(humdiv.train, "HumDiv", "Training"))
	df <- rbind(df, getDatasetTableRow(humdiv.val, "HumDiv", "Validation"))
	df <- rbind(df, getDatasetTableRow(clinvar.hv, "ClinVar", "Validation"))
	
	return(df)
}

# Actually build the dataset table
datasetTable <- buildDatasetTable()
datasetTable # Display to console

# Write dataset table to a file
write.table(datasetTable, "Table-Datasets.txt", sep="\t", quote=F)