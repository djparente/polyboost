# Daniel J. Parente, MD PhD
# University of Kansas Medical Center
# 2020-09-28

# Calculates the RMSE on 10-bin calibration plots, including calculation of 95% CIs by bootstrap analysis

# Declare libraries
library(boot)

# Calculates the RMSE of a data frame given cutpoints
rmsdDF <- function(dat, cuts) {
	dat <- data.frame(dat)
	
	# Quantize the predictions into the cut point bins
	dat$quantPred <- cut(dat$pred, seq(0, 1, 0.1), include.lowest=T, right=F)
	
	# Quantize the truth into 0 or 1. This has the useful property that mean(quantTruth) over some range
	# of interest (i.e. datapoints falling into a quantPred bins) is the percentage of deleterious variants
	# in the bin
	dat$quantTruth <- ifelse(dat$truth == 'd', 1, 0)
	
	# Aggregate means and medians, stratified by quantPred
	ag <- aggregate(. ~ quantPred, dat, FUN=function(x) c(mean=mean(x), median=median(x)))
	
	# quantTruth[,1] is the mean of quantTruth within a bin (i.e. the percent deleterious) while pred[,2] is the 
	# median of the predictions (i.e. the characteristic prediction of the bin). The sqrt of the mean of the square
	# difference is the RMSD, the quantity we want
	rmsd <- sqrt(mean((ag$quantTruth[,1] - ag$pred[,2])^2))
	return(rmsd)
}

# Calculates the RMSE of an individual model
rmsdModel <- function(model, valid, method, cuts=seq(0, 1, 0.1), preds=NULL) {
	# Get predictions from the model and validation set, if not explicitly specified in the function call (needed for PolyPhen-2)
	if(missing('preds')) {
		preds <- predict(model, valid, type="prob", na.action=na.pass)[,1] # Get deleterious probability
	}
	
	# Set up a data frame relating ground truth to predictions
	predframe <- data.frame(truth=valid$truth, pred=preds)
	
	# Calculate the RMSD of this data frame
	return(rmsdDF(predframe, cuts))
}

# Helper function for bootstrap analysis; calculates RMSE on indices from data passed by boot
rmsdModelBootstrap <- function(cuts, data, indices) {
	d <- data[indices,]
	return(rmsdDF(d, cuts))
}

# Calculates the RMSE of a model including error from boostrap analysis (95% CI; lci = 2.5% bound, uci = 97.5% bound)
rmsdModelWithError <- function(model, valid, method, cuts=seq(0, 1, 0.1), preds=NULL) {
	# Get predictions from the model and validation set, if not explicitly specified in the function call (needed for PolyPhen-2)
	if(missing('preds')) {
		preds <- predict(model, valid, type="prob", na.action=na.pass)[,1] # Get deleterious probability
	}
	
	# Set up a data frame relating ground truth to predictions
	predframe <- data.frame(truth=valid$truth, pred=preds)
	
	# Calculate bootstrap object
	bootobj <- boot(predframe, statistic=rmsdModelBootstrap, R=10000, cuts=cuts)
	# Use the bootstrap object to get a 95% CI
	ci <- boot.ci(bootobj, type="norm")
	
	# Return a data frame showing the method, point estimate and 95% CI, for inclusion in a larger table
	return(data.frame(method=method, rmsd=ci$t0, lci=ci$normal[1,2],  uci=ci$normal[1,3]))
}

# Builds the RMSD tables
getRMSDModelTable <- function(tableType) {
	# Declare an empty data frame to store results in
	res <- data.frame()

	# Switch on the type of table desired (HumVar, HumDiv, ClinVar) and build the table
	# by calling rmsdModelWithError for each calibration RMSD we want to include in the table
	if( tableType == "HumVar" ) {
		res <- rbind(res, rmsdModelWithError(tunes.polyboost.hv, humvar.val, "XGBoost"))
		res <- rbind(res, rmsdModelWithError(tunes.rf.hv, humvar.val, "Random Forest"))
		res <- rbind(res, rmsdModelWithError(tunes.svm.hv, humvar.val, "SVM"))
		res <- rbind(res, rmsdModelWithError(tunes.ann.hv, humvar.val, "Neural network"))
	} else if (tableType == "HumDiv") {
		res <- rbind(res, rmsdModelWithError(tunes.polyboost.hd, humdiv.val, "XGBoost"))
		res <- rbind(res, rmsdModelWithError(tunes.rf.hd, humdiv.val, "Random Forest"))
		res <- rbind(res, rmsdModelWithError(tunes.svm.hd, humdiv.val, "SVM"))
		res <- rbind(res, rmsdModelWithError(tunes.ann.hd, humdiv.val, "Neural network"))
	} else if (tableType == "ClinVar") {
		res <- rbind(res, rmsdModelWithError(tunes.polyboost.hv, clinvar.hv, "XGBoost"))
		res <- rbind(res, rmsdModelWithError(tunes.rf.hv, clinvar.hv, "Random Forest"))
		res <- rbind(res, rmsdModelWithError(NULL, clinvar.hv, preds=clinvar.hv$pph2_prob, "PolyPhen-2/HumVar"))
		res <- rbind(res, rmsdModelWithError(NULL, clinvar.hd, preds=clinvar.hd$pph2_prob, "PolyPhen-2/HumDiv"))
	}

	return(res)
}

# Get the tables
rmsd.clinvar <- getRMSDModelTable("ClinVar")
rmsd.hv <- getRMSDModelTable("HumVar")
rmsd.hd <- getRMSDModelTable("HumDiv")

# Write the tables out to the console
print("ClinVar RMSD Table")
rmsd.clinvar

print("HumVar RMSD Table")
rmsd.hv

print("HumDiv RMSD Table")
rmsd.hd

# Write the tables out to the disk
write.table(rmsd.clinvar, "rmsd-clinvar.txt", sep="\t", quote=F)
write.table(rmsd.hv, "rmsd-humvar.txt", sep="\t", quote=F)
write.table(rmsd.hd, "rmsd-humdiv.txt", sep="\t", quote=F)

