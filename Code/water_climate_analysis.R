# water_climate_analysis.R analyzes data for the snow persistence and watershed signatures project
# Authors: Kristo Elijah Krugger (Architecture/base code), Edward Le (Customization/application to SP project)
# Statistical/Code Validation By: Joseph Janssen          

# Note: By default, only saves Johnson-Neyman Plots for p < 0.01. Please see summarize_results() for commentted out 
#       code for alternative plotting methods (e.g., also including p < 0.05, simple slopes analysis and plot
#       /model diagnostics in the same plot)
#       Model diagnostics are also produced by default (but saved to seperate folder)
#		The call to print JN objects (e.g., "print(jnPlotsP001)") is commented out but can be uncommented to see numerical intervals

# *** Change log ***
# July 26th, 2021 Changes (+/- a week): Refactored code to create reusable function for
# collating statistical analyses, provides optimizations due to uniform transformation
# of climate variables (aridity index only)
# Purpose of refactoring was to find ways to efficiently extract statistical model
# results and report them on for the paper
# Further experimentation may require older versions of the code or manual manipulation without use of the functions
#       Note: Used <<- in collate_analyses function (will assign globally regardless of
#       parameter variable names, an unfortunate side-effect of using the R language)
# July 28th, 2021
# Added comments to transformation section
# Jan 9th, 2022
# Changed filepaths and change pdf output to png output
# Jan 11th, 2022
# Added standardization of y-axis to Johnson-Neyan plots, per type of plot (e.g., p < 0.01; p < 0.05)
# Jan 16th, 2022
# Added double bracketing of some variables to reduce silent errors (refactoring)
# Jan 20th, 2022
# Added creation of plotting folder to make code run standalone
# Jan 27th, 2022
# Making code more standalone (running based on relative directories)
# Jan 30th, 2022
# Cleanup of code, renaming j/jn/j2/jn2 variables to more expressive variables
# Feb 2nd, 2022
# Added functionality to plot all figures at once
# Feb 17th-19th, 2022
# Cleanup of signatures

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(readxl)
library(GGally)
library(repr)
library(cowplot)
library(broom)
library(car)
require(ggiraph)
require(ggiraphExtra)
library(ggResidpanel)
library(interactions)
library(effects)
library(stringr)
library(psych)
library(extrafont)
library(sandwich)
library(ggstance)
library(broom.mixed)
library(iml)
library(gridExtra)
library(collections)
# Customize function meanings for code (suggestions from Kristo Ellijah Krugger)
filter <- dplyr::filter
select <- dplyr::select
summarise_all <- dplyr::summarise_all
funs <- dplyr::funs


## Stage 1: Preparation and Staging of Data

# Code Variables
# Plot titles 
titleDict <- dict()
titleDict$set('mean_si', "SI")

titleDict$set("bfi", dict(list(manuscript = "(a) BFI",
	supplementary =  "BFI Model Diagnostics",
 	spYAxis = "Effect of SP on BFI",
	jnOrder = "a")))
titleDict$set("cbrt_q5Frac", dict(list(manuscript = expression("(b) Normalized"~Q[5]*"\n"), 
	supplementary = expression("Normalized"~Q[5]~"Model Diagnostics"),
	spYAxis = expression("Effect of SP on Normalized"~Q[5]),
	jnOrder = "b")))
titleDict$set("log_low_fdc", dict(list(manuscript = "(c) Low-FDC",
	supplementary = "Low-FDC Model Diagnostics",
	spYAxis = "Effect of SP on Low-FDC",
	jnOrder = "c")))
titleDict$set("log_q95Frac", dict(list(manuscript = expression("(d) Normalized"~Q[95]*"\n"),
	supplementary = expression("Normalized"~Q[95]~"Model Diagnostics"),
	spYAxis = expression("Effect of SP on Normalized"~Q[95]),
	jnOrder = "d")))
titleDict$set("sqrt_high_fdc", dict(list(manuscript = "(e) High-FDC",
	supplementary = "High-FDC Model Diagnostics",
	spYAxis = "Effect of SP on High-FDC",
	jnOrder = "e")))

# Ordering utility dictionary
figureOrder <- "fig03"
orderDict <- dict()
orderDict$set("bfi", dict(list(manuscript = "a", supplementary = "1")))
orderDict$set("cbrt_q5Frac", dict(list(manuscript = "b", supplementary = "2")))
orderDict$set("log_low_fdc", dict(list(manuscript = "c", supplementary = "3")))
orderDict$set("log_q95Frac", dict(list(manuscript = "d", supplementary = "4")))
orderDict$set("sqrt_high_fdc",dict(list(manuscript = "e", supplementary = "5")))

# Create dataframes for exported statistics
rSquares <- data.frame()
hStats <- data.frame()
modelCoefficients <- data.frame()
combinedStatsDf <- data.frame()

# Directory management
# Setwd (if necessary; depends on IDE)
# setwd(<fill_in>)
masterDir <- '../'
dataDir <- paste(masterDir, 'Data/', sep ="")
plotDir <- paste(masterDir, 'Plots/', sep ="")
diagPlotDir <- paste(plotDir, "Model Diagnostics/", sep = "")
jnPlotDir <- paste(plotDir, "JN Plots/", sep = "")
statisticsDir <- paste(masterDir, "Statistics/", sep="")

# Data name (presumes data folder already exists)
combinedDataName <- 'water_climate_data.csv'

# Read data
raw.data <- read_csv(paste(dataDir, combinedDataName, sep=""))

# Create plots folder (if it does not already exist)
dir.create(plotDir, showWarnings = FALSE)
dir.create(diagPlotDir, showWarnings = FALSE)
dir.create(jnPlotDir, showWarnings = FALSE)
dir.create(statisticsDir, showWarnings = FALSE)

# Get names of predictor columns
# First three cols catchment number, longitude, & latitude
predictors <- names(raw.data)[c(4:6)]

# Get names of outcome columns
outcomes <- names(raw.data)[c(7:length(names(raw.data)))]


# Select the predictor and outcome columns
raw.data <- raw.data %>% select(all_of(c(predictors, outcomes)))
data.filtered <- raw.data

# Drop NA data
data.filtered <- drop_na(data.filtered)

# Transformations for some predictors and responses
data.trans <- data.filtered %>% mutate(
						# Snow persistence
						log_mean_sp = log(mean_sp),
						sqrt_mean_sp = sqrt(mean_sp),
						cbrt_mean_sp = mean_sp^(1/3),
						# Aridity Index
						log_mean_ai = log(mean_ai),
						sqrt_mean_ai = sqrt(mean_ai),
						cbrt_mean_ai = mean_ai^(1/3),
						# Baseflow Index
						log_bfi = log(bfi),
						sqrt_bfi = sqrt(bfi),
						cbrt_bfi = bfi^(1/3),
						# q5Frac
						log_q5Frac = log(q5Frac),
						sqrt_q5Frac = sqrt(q5Frac),
						cbrt_q5Frac = (q5Frac)^(1/3),
						# Low-FDC
						log_low_fdc = log(low_fdc),
						sqrt_low_fdc = sqrt(low_fdc),
						cbrt_low_fdc = low_fdc^(1/3),
						# q95 Frac
						log_q95Frac = log(q95Frac),
						sqrt_q95Frac = sqrt(q95Frac),
						cbrt_q95Frac = (q95Frac)^(1/3),
						# High-FDC
						log_high_fdc = log(high_fdc),
						sqrt_high_fdc = sqrt(high_fdc),
						cbrt_high_fdc = high_fdc^(1/3),
						)

# Get lower quartile (25%) of mean_ai
aridIdx <- "mean_ai"
AI_25 <- quantile(data.trans[[aridIdx]], na.rm = TRUE)[['25%']]
# Get median of mean_ai
AI_50 <- quantile(data.trans[[aridIdx]], na.rm = TRUE)[['50%']]
# Get upper quartile (75%) of mean_ai
AI_75 <- quantile(data.trans[[aridIdx]], na.rm = TRUE)[['75%']]
AI_vals <- c(AI_25,AI_50,AI_75)

# Get lower quartile (25%) of mean_si
sI <- "mean_si"
SI_25 <- quantile(data.trans[[sI]], na.rm = TRUE)[['25%']]
# Get median of mean_si
SI_50 <- quantile(data.trans[[sI]], na.rm = TRUE)[['50%']]
# Get upper quartile (75%) of mean_si
SI_75 <- quantile(data.trans[[sI]], na.rm = TRUE)[['75%']]
SI_vals <- c(SI_25,SI_50,SI_75)

# Get maximum and minimum AI values
min_AI <- min(data.trans[[aridIdx]])
max_AI <- max(data.trans[[aridIdx]])

# Get values for log_ai and log_si
logAridIdx <- "log_mean_ai"
logAI_25 <- quantile(data.trans[[logAridIdx]], na.rm = TRUE)[['25%']]
# Get median of mean_ai
logAI_50 <- quantile(data.trans[[logAridIdx]], na.rm = TRUE)[['50%']]
# Get upper quartile (75%) of mean_ai
logAI_75 <- quantile(data.trans[[logAridIdx]], na.rm = TRUE)[['75%']]
logAI_vals <- c(logAI_25,logAI_50,logAI_75)
min_logAI <- min(data.trans[[logAridIdx]])
max_logAI <- max(data.trans[[logAridIdx]])


# Get values for sqrt_ai and sqrt_si
sqrtAridIdx <- "sqrt_mean_ai"
sqrtAI_25 <- quantile(data.trans[[sqrtAridIdx]], na.rm = TRUE)[['25%']]
# Get median of mean_ai
sqrtAI_50 <- quantile(data.trans[[sqrtAridIdx]], na.rm = TRUE)[['50%']]
# Get upper quartile (75%) of mean_ai
sqrtAI_75 <- quantile(data.trans[[sqrtAridIdx]], na.rm = TRUE)[['75%']]
sqrtAI_vals <- c(sqrtAI_25,sqrtAI_50,sqrtAI_75)
min_sqrtAI <- min(data.trans[[sqrtAridIdx]])
max_sqrtAI <- max(data.trans[[sqrtAridIdx]])

# Get values for cbrt_ai and cbrt_si
cbrtAridIdx <- "cbrt_mean_ai"
cbrtAI_25 <- quantile(data.trans[[cbrtAridIdx]], na.rm = TRUE)[['25%']]
# Get median of mean_ai
cbrtAI_50 <- quantile(data.trans[[cbrtAridIdx]], na.rm = TRUE)[['50%']]
# Get upper quartile (75%) of mean_ai
cbrtAI_75 <- quantile(data.trans[[cbrtAridIdx]], na.rm = TRUE)[['75%']]
cbrtAI_vals <- c(cbrtAI_25,cbrtAI_50,cbrtAI_75)
min_cbrtAI <- min(data.trans[[cbrtAridIdx]])
max_cbrtAI <- max(data.trans[[cbrtAridIdx]])

## Stage 2: Outlining of methods and functions used

# Plot slopes at different levels of interaction
# Helper function to help see relationships (not required)
plot_interaction <- function(model, var, 
							 inter, inter2,
							 AIvals, SIvals,
							 minAI, maxAI){
	# to adjust plots size
	options(repr.plot.width=20, repr.plot.height=20)
	var <- rlang::enquo(var)
	inter <- rlang::enquo(inter)
	inter2 <- rlang::enquo(inter2)
	int_plot <- interact_plot(model, pred = !!var, 
							  modx = !!inter, mod2=!!inter2,
							  modx.values=AIvals, mod2.values=SIvals,
							  facet.modx = FALSE, control.fdr = TRUE,
							  plot.points = TRUE)
	int_plot + theme_gray() + theme(axis.title = element_text(family = "sans", ),
		   legend.text = element_text(family = "sans", ),
		   legend.title = element_text(family = "sans", ),
		   strip.text = element_text(family = "sans", ),
								   axis.text.x = element_text(), 
											 axis.text.y = element_text())
}

# Helper to create plot title objects (for all plots)
# Dependencies: titleDict, orderDict properly encoded with titles/order, respectively
# Order is one of NULL, "manuscript", or "supplementary"
create_plot_super_titles <- function(newTitle){
		superTitle <- ggdraw() + 
        	draw_label(newTitle,
                   x = 0,
                   hjust = 0) +
        	theme(plot.margin = margin(0,0,0,7))
	return(superTitle)
}

# Helper to create plot titles for variable title, generally rounded to two decimal places
# Dependencies: Global dict for variable names, varaibleName is assumed to be quosure
create_plot_variable_title_string <- function(varVal, variableName, roundedDigits = 2){
	variableNameString <- quo_name(variableName)
	variableValRounded <- format(round(varVal, roundedDigits), nsmall = roundedDigits)
	titleName <- paste(titleDict$get(variableNameString), " = ", variableValRounded)
	return(titleName)
}

# Helper to create figure labelling
# Dependencies: orderDict properly encoded with order
# Order is one of NULL, "manuscript", or "supplementary"
create_plot_filename <- function(plotDir, sigName, prefixString = "", suffixString = "", order = NULL){
	plotFileName <- paste(plotDir, prefixString, as.character(sigName), suffixString, sep ="")
	if (!is.null(order)){
		plotFileName <- paste(plotDir, prefixString,
		 orderDict$get(sigName)$get(order), "_",
		 as.character(sigName),
		 suffixString, sep ="")
	}
	return (plotFileName)
}

# Plot confidence interval of slopes (simple slopes analysis)
plot_sim_slopes_analysis <- function(model, var,
									 inter, inter2,
									 inter1_vals, inter2_vals){
	
	var <- rlang::enquo(var)
	inter <- rlang::enquo(inter)
	inter2 <- rlang::enquo(inter2)
	sp <- plot(sim_slopes(model, pred = !!var,
						  modx = !!inter, mod2 = !!inter2, 
						  modx.values = inter1_vals, mod2.values = inter2_vals,
						  johnson_neyman = FALSE, 
						  control.fdr = TRUE))
	sp + theme_gray()  + theme(axis.title = element_text(family = "sans", ),
		   legend.text = element_text(family = "sans", ),
		   legend.title = element_text(family = "sans", ),
		   strip.text = element_text(family = "sans", ),
		   axis.text.x = element_text(), 
		   axis.text.y = element_text(), 
		   axis.title.y = element_blank())
	
}

# Create and Format Johnson-Neyman Plots, save diagnostic plots
# Returns JN Plots
# var, inter, inter2 are treated as quosures
summarize_results <- function(model, var, plotName,
							  inter, inter1_vals, 
							  inter2, inter2_vals,
							  minAI, maxAI){
	# Adjust plot sizes
	options(repr.plot.width=20, repr.plot.height=23)
	
	var <- rlang::enquo(var)
	inter <- rlang::enquo(inter)
	inter2 <- rlang::enquo(inter2)
	
	# Create Johnson-Neyman analysis and plots (p < 0.05)
	jnP005 <- sim_slopes(model, pred = !!var, modx = !!inter, mod2 = !!inter2, 
					modx.values = inter1_vals, mod2.values = inter2_vals,
					johnson_neyman = TRUE, jnplot = TRUE, interval = T, control.fdr = TRUE,
					mod.range= c(minAI, maxAI))
	# Plots for p < 0.05
	jnPlotsP005 <- jnP005$jn
	# Get max and min ranges for y-axis; p < 0.05
	slopeMinP005 <- min(c(min(jnPlotsP005[[1]]$cbands$Lower),
						  min(jnPlotsP005[[2]]$cbands$Lower),
						  min(jnPlotsP005[[3]]$cbands$Lower)))
	slopeMaxP005 <- max(c(max(jnPlotsP005[[1]]$cbands$Upper),
						  max(jnPlotsP005[[2]]$cbands$Upper),
						  max(jnPlotsP005[[3]]$cbands$Upper)))
	# Ensure zero line is shown when slope max is below zero
	if (slopeMaxP005 < 0){
		slopeMaxP005 = -1 * slopeMaxP005
	}
	# Ensure zero line is shown when slope min is below zero
	if (slopeMinP005 > 0){
		slopeMinP005 = -1 * slopeMinP005
	}
	# Third element is the smallest in the list of moderators
	p1 <- jnPlotsP005[[3]]$plot  +
		# Assume first elem is the smallest
		ggtitle(create_plot_variable_title_string(inter2_vals[1], inter2)) +
		scale_y_continuous(
			name = titleDict$get(plotName)$get("spYAxis"),
			limits = c(slopeMinP005, slopeMaxP005)) +
		xlab("Ln(AI)") +
		geom_rect(aes(xmin = 0, xmax= Inf, ymin= -Inf, ymax= Inf), fill = 'beige', alpha=0.1) +
		geom_rect(aes(xmin = -Inf, xmax= 0, ymin= -Inf, ymax= Inf), fill = 'darkslategray1', alpha=0.1) +
		annotate(geom = 'text', label = 'Wet Climate', x = minAI/2, y = Inf, hjust = 0.5, vjust = 1, size = 3) +
		annotate(geom = 'text', label = 'Dry Climate', x = maxAI/2, y = Inf, hjust = 0.5, vjust = 1, size = 3) +
		geom_vline(xintercept = 0) +
		theme_gray() +
		theme(axis.title = element_text(family = "sans", ),
											 strip.text = element_text(family = "sans", ),
											 axis.text.x = element_text(),
											 axis.text.y = element_text(),
											 plot.title = element_text(family = "sans",
																	   hjust = 0.5),
											 plot.margin=margin(l=0.1, unit = 'cm'),
											 legend.position="none")
	p2 <- jnPlotsP005[[2]]$plot +
		ggtitle(create_plot_variable_title_string(inter2_vals[2], inter2)) +
		scale_y_continuous(
			name = titleDict$get(plotName)$get("spYAxis"),
			limits = c(slopeMinP005, slopeMaxP005)) +
		xlab("Ln(AI)") +
		geom_rect(aes(xmin = 0, xmax= Inf, ymin= -Inf, ymax= Inf), fill = 'beige', alpha=0.1) +
		geom_rect(aes(xmin = -Inf, xmax= 0, ymin= -Inf, ymax= Inf), fill = 'darkslategray1', alpha=0.1) +
		annotate(geom = 'text', label = 'Wet Climate', x = minAI/2, y = Inf, hjust = 0.5, vjust = 1, size = 3) +
		annotate(geom = 'text', label = 'Dry Climate', x = maxAI/2, y = Inf, hjust = 0.5, vjust = 1, size = 3) +
		geom_vline(xintercept = 0) +
		theme_gray() + 
		theme(axis.title = element_text(family = "sans", ),
											 strip.text = element_text(family = "sans", ),
											 axis.text.x = element_text(), 
											 axis.text.y = element_blank(),
											 plot.title = element_text(family = "sans",
																	   hjust = 0.5 ),
											 axis.title.y = element_blank(),
											 plot.margin=margin(l=0.1, unit = 'cm'),
											 legend.position="none")
	p3 <- jnPlotsP005[[1]]$plot +
		ggtitle(create_plot_variable_title_string(inter2_vals[3], inter2)) +
		scale_y_continuous(
			name = titleDict$get(plotName)$get("spYAxis"),
			limits = c(slopeMinP005, slopeMaxP005)) +
		xlab("Ln(AI)") +
		geom_rect(aes(xmin = 0, xmax= Inf, ymin= -Inf, ymax= Inf), fill = 'beige', alpha=0.1) +
		geom_rect(aes(xmin = -Inf, xmax= 0, ymin= -Inf, ymax= Inf), fill = 'darkslategray1', alpha=0.1) +
		annotate(geom = 'text', label = 'Wet Climate', x = minAI/2, y = Inf, hjust = 0.5, vjust = 1, size = 3) +
		annotate(geom = 'text', label = 'Dry Climate', x = maxAI/2, y = Inf, hjust = 0.5, vjust = 1, size = 3) +
		geom_vline(xintercept = 0) +
		theme_gray() + 
		theme(axis.title = element_text(family = "sans", ),
											 strip.text = element_text(family = "sans", ),
											 axis.text.x = element_text(), 
											 axis.text.y = element_blank(),
											 plot.title = element_text(family = "sans",
																	   hjust = 0.5),
											 plot.margin=margin(l=0.1, unit = 'cm'),
											 axis.title.y = element_blank(),
											 legend.text = element_text(family = "sans",
																		size = 7))  +
		labs(fill="")
	# NOTE: For publication, we are using p < 0.01 only
	# Create Johnson-Neyman analysis and plots (p < 0.01)
	jnP001 <- sim_slopes(model, pred = !!var, modx = !!inter, mod2 = !!inter2, 
					modx.values = inter1_vals, mod2.values = inter2_vals,
					johnson_neyman = TRUE, jnplot = TRUE, interval = T, control.fdr = TRUE, jnalpha=0.01,
					mod.range= c(minAI, maxAI))
	jnPlotsP001 <- jnP001$jn
	# Get max and min ranges for y-axis; p < 0.01
	slopeMinP001 <- min(c(min(jnPlotsP001[[1]]$cbands$Lower),
					  min(jnPlotsP001[[2]]$cbands$Lower),
					  min(jnPlotsP001[[3]]$cbands$Lower)))
	slopeMaxP001 <- max(c(max(jnPlotsP001[[1]]$cbands$Upper),
					  max(jnPlotsP001[[2]]$cbands$Upper),
					  max(jnPlotsP001[[3]]$cbands$Upper)))
	# Ensure zero line is shown when slope max is below zero
	if (slopeMaxP001 < 0){
		slopeMaxP001 = -1 * slopeMaxP001
	}
	# Ensure zero line is shown when slope min is below zero
	if (slopeMinP001 > 0){
		slopeMinP001 = -1 * slopeMinP001
	}
	# Format plotting
	p1b <- jnPlotsP001[[3]]$plot  +
		ggtitle(create_plot_variable_title_string(inter2_vals[1], inter2)) +
		scale_y_continuous(
			name = titleDict$get(plotName)$get("spYAxis"),
			limits = c(slopeMinP001, slopeMaxP001)) +
		xlab("Ln(AI)") +
		geom_rect(aes(xmin = 0, xmax= Inf, ymin= -Inf, ymax= Inf), fill = 'beige', alpha=0.1) +
		geom_rect(aes(xmin = -Inf, xmax= 0, ymin= -Inf, ymax= Inf), fill = 'darkslategray1', alpha=0.1) +
		annotate(geom = 'text', label = 'Wet Climate', x = minAI/2, y = Inf, hjust = 0.5, vjust = 1, size = 3) +
		annotate(geom = 'text', label = 'Dry Climate', x = maxAI/2, y = Inf, hjust = 0.5, vjust = 1, size = 3) +
		geom_vline(xintercept = 0) +
		theme_gray() + 
		theme(axis.title = element_text(family = "sans", ),
											 strip.text = element_text(family = "sans",
																	   ),
											 axis.text.x = element_text(), 
											 axis.text.y = element_text(),
											 plot.title = element_text(family = "sans",
																	   hjust = 0.5),
											 plot.margin=margin(l=0.1, unit = 'cm'),
											 legend.position="none")
	p2b <- jnPlotsP001[[2]]$plot  +
		ggtitle(create_plot_variable_title_string(inter2_vals[2], inter2)) +
		scale_y_continuous(
			name = titleDict$get(plotName)$get("spYAxis"),
			limits = c(slopeMinP001, slopeMaxP001)) +
		xlab("Ln(AI)") +
		geom_rect(aes(xmin = 0, xmax= Inf, ymin= -Inf, ymax= Inf), fill = 'beige', alpha=0.1) +
		geom_rect(aes(xmin = -Inf, xmax= 0, ymin= -Inf, ymax= Inf), fill = 'darkslategray1', alpha=0.1) +
		annotate(geom = 'text', label = 'Wet Climate', x = minAI/2, y = Inf, hjust = 0.5, vjust = 1, size = 3) +
		annotate(geom = 'text', label = 'Dry Climate', x = maxAI/2, y = Inf, hjust = 0.5, vjust = 1, size = 3) +
		geom_vline(xintercept = 0) +
		theme_gray() +  
		theme(axis.title = element_text(family = "sans", ),
											 strip.text = element_text(family = "sans", ),
											 axis.text.x = element_text(), 
											 axis.text.y = element_blank(),
											 plot.title = element_text(family = "sans",
																		hjust = 0.5),
											 axis.title.y = element_blank(),
											 plot.margin=margin(l=0.1, unit = 'cm'),
											 legend.position="none")
	p3b <- jnPlotsP001[[1]]$plot  +
		ggtitle(create_plot_variable_title_string(inter2_vals[3], inter2)) +
		scale_y_continuous(
			name = titleDict$get(plotName)$get("spYAxis"),
			limits = c(slopeMinP001, slopeMaxP001)) +
		xlab("Ln(AI)") +
		geom_rect(aes(xmin = 0, xmax= Inf, ymin= -Inf, ymax= Inf), fill = 'beige', alpha=0.1) +
		geom_rect(aes(xmin = -Inf, xmax= 0, ymin= -Inf, ymax= Inf), fill = 'darkslategray1', alpha=0.1) +
		annotate(geom = 'text', label = 'Wet Climate', x = minAI/2, y = Inf, hjust = 0.5, vjust = 1, size = 3) +
		annotate(geom = 'text', label = 'Dry Climate', x = maxAI/2, y = Inf, hjust = 0.5, vjust = 1, size = 3) +
		geom_vline(xintercept = 0) +
		theme_gray() +
		theme(axis.title = element_text(family = "sans",),
											 strip.text = element_text(family = "sans", ),
											 axis.text.x = element_text(), 
											 axis.text.y = element_blank(),
											 plot.title = element_text(family = "sans",
																	   hjust = 0.5),
											 plot.margin=margin(l=0.1, unit = 'cm'),
											 axis.title.y = element_blank(),
											 legend.text = element_text(family = "sans",
																		size = 7)) +
		labs(fill="")
	# Print jnPlotsP001 object (to see values for aridity/see JN intervals)    
	# print(jnPlotsP001)
	# Create plot diagnostic plots
	pDiagName <- create_plot_filename(plotDir = diagPlotDir,
										sigName = plotName,
										prefixString = "figS",
										suffixString =  "_model_diagnostics.png",
										order = "supplementary")
	residualPanels <- resid_panel(model, plots = c("resid", "qq", "hist"), theme="grey", nrow=1, title.opt = FALSE)
	# Residual title
	residualsTitle <- ggdraw() + 
        draw_label("(a) Residual Plot",
                   x = 0,
                   hjust = 0) +
        theme(plot.margin = margin(0,0,0,7))
	# QQ Title
	qqTitle <- ggdraw() + 
        draw_label("(b) QQ Plot",
                   x = 0,
                   hjust = 0) +
        theme(plot.margin = margin(0,0,0,7))
	# Histogram Title
	histTitle <- ggdraw() + 
        draw_label("(c) Histogram",
                   x = 0,
                   hjust = 0) +
        theme(plot.margin = margin(0,0,0,7))
	residualPanelTitles <- plot_grid(residualsTitle, qqTitle, histTitle, rel_widths = c(0.33, 0.33, 0.33), nrow = 1)

	pDiag <- plot_grid(create_plot_super_titles(titleDict$get(plotName)$get("supplementary")),
						residualPanelTitles,
						residualPanels,
						nrow = 3,
						rel_heights = c(0.1, 0.1, 1))
	save_plot(pDiagName,
			  pDiag,
			  base_asp = 1.85,
			  dpi = 320)
	# Create JN plots
	# Create name of JN plot:
	pJnName <- create_plot_filename(plotDir = jnPlotDir,
										sigName = plotName,
										prefixString = figureOrder,
										suffixString =  "_model_plot_three_interact.png",
										order = "manuscript")
	pJnAllTitle <- create_plot_super_titles(titleDict$get(plotName)$get("manuscript"))
	# JN plots without diagnostics (p < 0.01 only)
	pJnAll <- plot_grid(pJnAllTitle,
					plot_grid(p1b, p2b, p3b,
					  rel_widths=c(1.1,0.80,1.23),
					  nrow=1,
					  scale = 0.95),
					nrow=2,
					rel_heights = c(0.1,1))
	# COMMENT OUT IF NOT p < 0.01 Custom plotting aspect ratio for p < 0.01
	save_plot(pJnName,
			  pJnAll,
			  base_asp = 1.85,
			  dpi = 320)
	# dev.off()
	# JN plots without diagnostics (with sim slopes, p < 0.01, p < 0.05)
	# pJnAll <- plot_grid(pJnAllTitle),
	#                   plot_interaction(model, !!var,
	#                                    inter=!!inter, inter2=!!inter2,
	#                                    AIvals = inter1_vals, SIvals = inter2_vals,
	#                                    minAI = minAI, maxAI = minAI),
	#                   plot_grid(p1, p2, p3, rel_widths=c(rel_widths=c(1,0.88,1.23)),nrow=1),
	#                   plot_grid(p1b, p2b, p3b, rel_widths=c(rel_widths=c(1,0.88,1.23)),nrow=1),
	#                   nrow= 4,
	#                   rel_heights = c(0.1,1,1,1),
	#                   scale = 0.95)
	# # JN plots with full diagnostics, p < 0.05, p < 0.01 (if required)
	# pJnAll <- plot_grid(pJnAllTitle),
	#                     plot_grid(
	#                             plot_interaction(model, !!var,
	#                                              inter=!!inter, inter2=!!inter2,
	#                                              AIvals = inter1_vals, SIvals = inter2_vals,
	#                                             minAI = minAI, maxAI = minAI),
	#                     plot_grid(p1, p2, p3, rel_widths=c(rel_widths=c(1,0.88,1.23)),nrow=1),
	#                     plot_grid(p1b, p2b, p3b, rel_widths=c(rel_widths=c(1,0.88,1.23)),nrow=1),
	#                     nrow=3, scale = 0.95),
	#         plot_grid(plot_sim_slopes_analysis(model, !!var,
	#                                            inter=!!inter, inter2=!!inter2,
	#                                            inter1_vals = inter1_vals,
	#                                            inter2_vals = inter2_vals),
	#                     resid_panel(model, theme="grey", nrow=2),
	#                     rel_widths=c(0.7,1), ncol=2),
	#         rel_heights= c(0.1,1,0.8),nrow=3)
	# # Common plot save (p > 0.01 + p > 0.05, p > 0.01 + p > 0.05 + diagnostics)
	# ggsave(filename = pJnName,
	#        plot = pJnAll,
	#        device= 'png',
	#        dpi = 320)
	return(pJnAll)
}

# Coordination function to extract R^2s, p-values, perform Johnson-Neymann analysis, perform h-stat analysis
# Returns JN plot
# Hard coded to consider log_mean_ai, mean_sp, and mean_si here. Many bugs introduced if used otherwise
# Pre-assumed to have min_logAI and max_logAI calculated as a global parameter
# <<- with rbind to work around variable mutability issues
collate_analyses <- function(currModel, currName, rSquares, hStats, modelCoefficients,
							 currMinAI = min_logAI,
							 currMaxAI = max_logAI){
	# Summary of model
	currModelSum <- summary(currModel)
	# Extract R^2 values and update global variable
	rSquares <<- rbind(rSquares, data.frame(signature = currName,
											   R2 = currModelSum$r.squared,
											   R2Adj = currModelSum$adj.r.squared))
	# Extract coefficients, reset indices, and update global variable
	tempCoeffs <- data.frame(currModelSum$coefficients)
	tempCoeffs <- cbind(coeff = rownames(tempCoeffs), tempCoeffs)
	rownames(tempCoeffs) <- 1:nrow(tempCoeffs)
	# Add signature name
	tempCoeffs <- data.frame(signature = currName, tempCoeffs)
	# Update global variable
	modelCoefficients <<- rbind(modelCoefficients, tempCoeffs)
	
	# Perform Johnson-Neymann analysis 
	pJnAllForCurrSig <- summarize_results(currModel, mean_sp, currName,                  
					  inter = log_mean_ai,
					  inter1_vals = logAI_vals,
					  inter2 = mean_si,
					  inter2_vals = SI_vals,
					  minAI = currMinAI,
					  maxAI = currMaxAI)
	# H-stat analysis (create model and extract h-values)
	currModelData <- currModel$model
	currModelImlObj <- Predictor$new(currModel, data = currModelData)
	currModelInteractions <- Interaction$new(currModelImlObj, feature = 'mean_sp')
	currHValues <- data.frame(currModelInteractions$results)
	currHValues <- cbind(signature = currName, currHValues)
	# Extract h-values and update global variable
	hStats <<- rbind(hStats, currHValues)
	return(pJnAllForCurrSig)
}

# Stage 3: Data Analysis and Modelling
jnPlotList <- list()

# BFI as untransformed
modelBFI <- lm(bfi ~
				   mean_sp*log_mean_ai*mean_si,
			   data=data.trans)
modelBFIName <- 'bfi'
modelBFIJnPlotOrder <- which(letters == titleDict$get(modelBFIName)$get("jnOrder"))
jnPlotList[[modelBFIJnPlotOrder]] <- collate_analyses(modelBFI, modelBFIName, rSquares, hStats, modelCoefficients)

# cbrt q5Frac
modelQ5Frac <- lm(cbrt_q5Frac ~
					   mean_sp*log_mean_ai*mean_si,
				  data=data.trans)
modelQ5FracName <- 'cbrt_q5Frac'
modelQ5FracPlotOrder <- which(letters == titleDict$get(modelQ5FracName)$get("jnOrder"))
jnPlotList[[modelQ5FracPlotOrder]] <- collate_analyses(modelQ5Frac, modelQ5FracName, rSquares, hStats, modelCoefficients)

#log low_fdc
modelLowFDC <- lm(log_low_fdc ~
					  mean_sp*log_mean_ai*mean_si,
			 data=data.trans)
modelLowFDCName <- 'log_low_fdc'
modelLowFDCPlotOrder <- which(letters == titleDict$get(modelLowFDCName)$get("jnOrder"))
jnPlotList[[modelLowFDCPlotOrder]]  <- collate_analyses(modelLowFDC, modelLowFDCName, rSquares, hStats, modelCoefficients)  

# log q95Frac
modelQ95Frac <- lm(log_q95Frac ~
					   mean_sp*log_mean_ai*mean_si,
				  data=data.trans)
modelQ95FracName <- 'log_q95Frac'
modelQ95FracPlotOrder <- which(letters == titleDict$get(modelQ95FracName)$get("jnOrder"))
jnPlotList[[modelQ95FracPlotOrder]] <- collate_analyses(modelQ95Frac, modelQ95FracName, rSquares, hStats, modelCoefficients) 

#sqrt high_fdc
modelHighFDC <- lm(sqrt_high_fdc ~
					   mean_sp*log_mean_ai*mean_si,
			 data=data.trans)
modelHighFDCName <- 'sqrt_high_fdc'
modelHighFDCPlotOrder <- which(letters == titleDict$get(modelHighFDCName)$get("jnOrder"))
jnPlotList[[modelHighFDCPlotOrder]] <- collate_analyses(modelHighFDC, modelHighFDCName, rSquares, hStats, modelCoefficients)

# Plot all panels as one 
allJnPlotsName <- paste(jnPlotDir, figureOrder, ".png", sep = "")
allJnPlots <- grid.arrange(grobs=jnPlotList, ncols = 2)
ggsave(allJnPlotsName,
          allJnPlots,
          height = 12,
          width = 15,
          dpi = 320)

# write csvs of analytical indices (R^2), H interaction statstics, and model coefficients
write.csv(rSquares, paste(statisticsDir, "signature_correlations_three_interact.csv", sep=""), row.names = FALSE)
write.csv(hStats, paste(statisticsDir, "h-statistic_interactions_strengths_three_interact.csv", sep=""), row.names = FALSE)
write.csv(modelCoefficients, paste(statisticsDir, "model_pValues_three_interact.csv", sep=""), row.names = FALSE)

