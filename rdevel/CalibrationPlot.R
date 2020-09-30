# Daniel Parente, MD PhD
# University of Kansas Medical Center

# Performs calibration analysis and generates calibration plots
# Can only be run after global_evaluation.R has already been completed


getCalibrationData <- function(dat, cuts) {
	dat <- data.frame(dat)
	
	# Quantize the predictions into the cut point bins
	dat$quantPred <- cut(dat$pred, seq(0, 1, 0.1), include.lowest=T, right=F)
	
	# Quantize the truth into 0 or 1. This has the useful property that mean(quantTruth) over some range
	# of interest (i.e. datapoints falling into a quantPred bins) is the percentage of deleterious variants
	# in the bin
	dat$quantTruth <- ifelse(dat$truth == 'd', 1, 0)
	
	# Aggregate means and medians, stratified by quantPred
	ag <- aggregate(. ~ quantPred, dat, FUN=function(x) c(mean=mean(x), median=median(x)))
	
	# quantTruth[,1] is the mean of quantTruth within each bin (i.e. the fraction deleterious) while pred[,1] is the mean and pred[,2] the median of the predictions within the bin
	return(data.frame(calib=ag$quantTruth[,1], means=ag$pred[,1], medians=ag$pred[,2]))
}

# Get a data frame returning the % actually deleterious, mean score and median score for various 
# bins ('cuts') of deleterious probabilities.
# N.B.: This is actually a pretty ugly solution that I developed first before realizing that
# this could be done in fewer lines (see above, based on rmseModelEvaluation.R's rmseDF function).
getCalibrationData_old <- function(dat, cuts) {
	res <- c()
	means <- c()
	medians <- c()
	

	# For each cut
	for(i in 1:(length(cuts) - 1)) {
		lLim <- cuts[i]
		uLim <- cuts[i+1]
		
		if(i == 1) {
			lLim = lLim - 0.0000001 # Subtract a small epsilon to allow greater than or equal to for the lowest cut
		}

		# Extract the part of the data table with the relevant probability range
		datInCut <- filter(dat, pred > lLim & pred <= uLim)
		
		# Determine the fraction that are deleterious within this probability range
		tableInCut <- table(datInCut$truth)
		dProb = tableInCut[1] / sum(tableInCut)
		
		# Add this to the growing calibration result
		res <- c(res, dProb)
		
		# Add the means and medians
		means <- c(means, mean(datInCut$pred))
		medians <- c(medians, median(datInCut$pred))
		
	}
	
	return(data.frame(calib = res, means=means, medians=medians))
}

# Get a calibration plot for a single model on one graph (not actually used to generate final figures)
calibrationFigSingle <- function(model, valid, preds=NULL) {
	# Use the model to obtain deleterious probabilities, if predictions are not explicitly passed (req'd for PolyPhen-2)
	if(missing('preds')) {
		preds <- predict(model, valid, type="prob", na.action=na.pass)[,1]
	}
	
	# Extract the calibration data, using the above method
	predframe <- data.frame(truth=valid$truth, pred=preds)
	calibDat <- getCalibrationData(predframe, seq(0, 1, 0.1))
	
	# Construct a plot
	res <- ggplot(calibDat, aes(x=medians, y=calib)) +
		geom_line() +
		geom_point() +		
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
		xlab('Deleterious probability (predicted)') +
		ylab('Deleterious rate (observed)') +
		geom_abline(intercept = 0, slope = 1, size=0.5)
	
	# Return the plot
	return(res)
}
# Get calibration data for a model in a data frame in preparation for plotting multiple lines on same plot
getCalibrationDataMeta <- function(model, valid, method, preds=NULL) {
	# Get predictions from the model and validation set, if not explicitly specified in the function call (needed for PolyPhen-2)
	if(missing('preds')) {
		preds <- predict(model, valid, type="prob", na.action=na.pass)[,1] # Get deleterious probability
	}
	
	# Set up a data frame relating ground truth to predictions
	predframe <- data.frame(truth=valid$truth, pred=preds)
	
	# Get the calibration data for cut points of 0-10%, 10-20%, etc.
	calibDat <- getCalibrationData(predframe, seq(0, 1, 0.1))
	
	# Append the name of the method (as a factor) to the data frame as a column (so we can rbind these into one larger frame with many methods)
	calibDat$method = method
	calibDat$method = as.factor(calibDat$method)
	
	return(calibDat)
}

# Construct a figure that has multiple calibration plots on it for the three possible validation datasets (HumVar, HumDiv, ClinVar)
calibrationFigMultiple <- function(figType) {
	# Set up a data frame for the three possible plots, HumVar, HumDiv or ClinVar
	calibDF <- data.frame()
	if( figType == "HumVar" ) {
		# Get calibration data for each model and rbind them to the global dataframe (calibDF)
		calibDF <- rbind(calibDF, getCalibrationDataMeta(tunes.polyboost.hv, humvar.val, "XGBoost"))
		calibDF <- rbind(calibDF, getCalibrationDataMeta(tunes.rf.hv, humvar.val, "Random Forest"))
		calibDF <- rbind(calibDF, getCalibrationDataMeta(tunes.svm.hv, humvar.val, "SVM"))
		calibDF <- rbind(calibDF, getCalibrationDataMeta(tunes.ann.hv, humvar.val, "Neural network"))
		panelannotation <- "Internal calibration: HumVar"
	} else if (figType == "HumDiv") {
		calibDF <- rbind(calibDF, getCalibrationDataMeta(tunes.polyboost.hd, humdiv.val, "XGBoost"))
		calibDF <- rbind(calibDF, getCalibrationDataMeta(tunes.rf.hd, humdiv.val, "Random Forest"))
		calibDF <- rbind(calibDF, getCalibrationDataMeta(tunes.svm.hd, humdiv.val, "SVM"))
		calibDF <- rbind(calibDF, getCalibrationDataMeta(tunes.ann.hd, humdiv.val, "Neural network"))
		panelannotation <- "Internal calibration: HumDiv"
	} else if (figType == "ClinVar") {
		calibDF <- rbind(calibDF, getCalibrationDataMeta(tunes.polyboost.hv, clinvar.hv, "XGBoost"))
		calibDF <- rbind(calibDF, getCalibrationDataMeta(tunes.rf.hv, clinvar.hv, "Random Forest"))
		calibDF <- rbind(calibDF, getCalibrationDataMeta(NULL, clinvar.hv, preds=clinvar.hv$pph2_prob, "PolyPhen-2/HumVar"))
		calibDF <- rbind(calibDF, getCalibrationDataMeta(NULL, clinvar.hd, preds=clinvar.hd$pph2_prob, "PolyPhen-2/HumDiv"))
		panelannotation <- "External calibration: ClinVar"
	} else {
		stop("Invalid internal validation figure type")
	}
	
	# Set the legend title
	legendtitle <- panelannotation
	
	# Create the figure
	res <- ggplot(calibDF, aes(x=medians, y=calib, group=method)) +
		xlim(0, 1) +
		ylim(0, 1) +
		geom_line(aes(linetype=method, color=method), size=1) + 
		geom_point(aes(shape=method, color=method), size=2) + 
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
		xlab('Deleterious probability (predicted)') +
		ylab('Deleterious rate (observed)') +
		geom_abline(intercept = 0, slope = 1, size=0.5) +
		scale_linetype_manual(values=c("dashed", "dashed", "dashed", "dashed")) + 
		scale_color_manual(values=c("#000000", "#FF0000", "#00FF00", "#0000FF")) + 
		scale_shape_manual(values=c(15, 16, 17, 18)) + 
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
		theme(legend.position = c(0.3, .8)) +
		

		labs(color=legendtitle, linetype=legendtitle, shape=legendtitle)
		
	# Return the figure and the calibration data frame
	return(list(fig=res, df=calibDF))
}

# Get the calibration figures and data frames
calibration.clinvar <- calibrationFigMultiple("ClinVar")
calibration.hv <- calibrationFigMultiple("HumVar")
calibration.hd <- calibrationFigMultiple("HumDiv")

# Save high resolution figures
ggsave("Calib-HV.png", calibration.hv$fig, height=6.5, width=6.5, unit="in", dpi=1200)
ggsave("Calib-HD.png", calibration.hd$fig, height=6.5, width=6.5, unit="in", dpi=1200)
ggsave("Calib-Clinvar.png", calibration.clinvar$fig, height=6.5, width=6.5, unit="in", dpi=1200)