# Purpose:  plotter_sites_signatures_climatic_metrics.R plots overall signature/climatic metric data for all available watersheds in the study
# Note: this version drops any watershed with na values
# Input:    - Overall signature/climatic metric data (assumed to come pre-loaded with latitudes and longitudes)
# Output:   - Plots for all signatures found within the trend and signature/climatic metric data
# Authors: Edward Le, Joe Janssen

# *** Change log ***
# August 3rd 2021
# Original version
# Jan 9th, 2022
# Changed output from pdf to png, changed file path locations
# Jan 10th, 2022
# Added plotting for combined plot
# Jan 13th, 2022
# Added code to handle change from "catchment" to "gridcode" marker
# January 26th, 2022
# Reintroduced log plotting (see code for handling/transformations)
# January 30th, 2022
# Added plain map (no colour bar for plotting)
# Added some functionality to auto label; similar to analysis code, with simpler code as plotting is simpler
# in this context (design decision)
# Feb 2nd, 2022
# Added functionality to plot all figures at once
# Feb 16-19th, 2022
# Clean up of signatures

library(ggplot2)
library(sf)
library(stringr)
library(rnaturalearth)
library(rnaturalearthdata)
library(readr)
library(sp)
library(maps)
library(ggspatial)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(gsubfn)
library(hash)
library(ggeasy)
library(patchwork)
library(ggsn)
library(BAMMtools)
library(ggnewscale)
library(scales)
library(viridis)
library(ggnewscale)
library(cowplot)
library(egg)

# Set working directory
setwd("./")

# Create plot folders
plotDir <- "../Plots/Geo/"
dir.create("../Plots", showWarnings = FALSE)
dir.create("../Plots/Geo", showWarnings = FALSE)

# Read attribute data
overallAttributes <- read_csv("../Data/water_climate_data.csv")
# Drop na
overallAttributes <- drop_na(overallAttributes)

northAmerica <- ne_states(c("united states of america", "canada"), returnclass = "sf")

panelOrderHash <- hash()
panelOrderHash[['mean_sp']] <- "a"
panelOrderHash[['mean_ai']] <- "b"
panelOrderHash[['mean_si']] <-"c"
panelOrderHash[["bfi"]] <-  "d"
panelOrderHash[["q5Frac"]] <- "e"
panelOrderHash[["low_fdc"]] <-  "f"
panelOrderHash[["q95Frac"]] <-  "g"
panelOrderHash[["high_fdc"]] <- "h"

# Titles List:
# Unfortunately, expression() is buggy with string/expression concatenation for Q5/Q95 subscripting; duplicate
# data on ordering here
plotTitlesHash <- hash()
plotTitlesHash[['mean_sp']] <- "(a) Snow persistence\n(% of Jan-July)"
plotTitlesHash[['mean_ai']] <- "(b) Aridity index"
plotTitlesHash[['mean_si']] <-"(c) Seasonality index"
plotTitlesHash[['bfi']] <- "(d) Baseflow index"
plotTitlesHash[['q5Frac']] <- expression("(e) Normalized"~Q[5],"")
plotTitlesHash[["low_fdc"]]<- "(f) Low flow duration curve slope\n(0.05-0.30)"
plotTitlesHash[['q95Frac']] <- expression("(g) Normalized"~Q[95],"")
plotTitlesHash[["high_fdc"]]<- "(h) High flow duration curve slope\n(0.70-0.95)"
plotTitlesHash[['sites']] <- "Map of study sites in North America"

# Figure prefix ordering
sitesPlotFigure <- "fig01"
attributeGeoPlotFigures <- "fig02"

# Helper functions 

# Create attribute plot filenames
# Dependencies: panelOrderHash encoded with panel orders
create_sig_filenames_string <- function(destPlotDir, sigName, prefixString, suffixString){
  paste(destPlotDir, prefixString, panelOrderHash[[sigName]] , "_", sigName, suffixString, sep="")
}

# Find minimum, midpoint, and maximum number from a distribution
# Will omit min and 25% to middle if they are the same (when rounded to two decimal places untransformed or log transform)
get_min_mids_max <- function(nums){
  min <- min(nums)
  quarterToMiddle <- min(nums) + ((max(nums)- min(nums))/4)
  middle <- min(nums) + ((max(nums)- min(nums))/2)
  threeQuartersToMiddle <- max(nums) - ((max(nums)- min(nums))/4)
  max <- max(nums)
  minMidsMax <- c(
              # Min
              min,
              # 25% to middle
              quarterToMiddle,
              # Middle
              middle,
              # 75%
              threeQuartersToMiddle,
              # Max
              max)
  if (round(min, 2) == round(quarterToMiddle, 2) | (round(exp(min), 2) == round(exp(quarterToMiddle), 2))) {
        minMidsMax <- c(
              # Min
              min,
              # Middle
              middle,
              # 75%
              threeQuartersToMiddle,
              # Max
              max)
  }
  return(minMidsMax)
}
# Formats min, mid, max with two decmial places
format_list_decimal_places <- function(nums, decPlaces = 2){
    nums <- lapply(nums, round , decPlaces)
    nums <- lapply(nums, format, nsmall = decPlaces)
    return(nums)   
}

# Procedural Code:

# Plotting sites only
# "cat" refers to the "category" of data (e.g., sites, signature, climatic metrics)
currCat <- "sites"
sitesPlot <- ggplot(data=northAmerica) +
  theme(plot.margin=grid::unit(c(-0.30,0,0,0), "null")) +
     geom_sf()+ xlab("Longitude") + ylab("Latitude") +
      theme(legend.text =  element_text(size = 11),
            legend.key.height = unit(0.5, "cm"),
            legend.spacing.y = unit(0.25, "cm")
            ) +
      coord_sf(xlim = c(-145, -50), ylim = c(25, 72), expand = FALSE) +
      annotation_scale(location = "bl", width_hint = 0.15) +
      annotation_north_arrow(location = "bl", which_north = "true",
                             pad_x = unit(0.08, "in"), pad_y = unit(0.25, "in"),
                             style = north_arrow_fancy_orienteering,
                             height = unit(1, "cm"),
                             width = unit(1, "cm")) +
      geom_point(data = overallAttributes,
                 na.rm = TRUE,
                 shape = 21,
                 fill = "black",
                 aes_string(x= "longitude",
                           y = "latitude"),
                inherit.aes = FALSE) +
      labs(fill = NULL) +
      # No ordering required for standalone plot
      plot_annotation(title = paste(plotTitlesHash[[currCat]], sep=""),
                      theme = theme(plot.title = element_text(hjust = 0.5,
                                                              size = 22,
                                                              margin = margin(b = 15, unit = "pt"))))
ggsave(filename = paste(plotDir, sitesPlotFigure, ".png", sep=""),
       plot = sitesPlot,
       device= "png",
       dpi = 320)

# Plot Characteristics
catPlotList <- list()

# First col of attribute data is catchment number, 2nd col is longitude, 3rd is latitude
categories <- colnames(overallAttributes)[4: length(colnames(overallAttributes))]
catLen <- length(categories)
for (i in 1:catLen){
  currCat <- categories[i]
  if (currCat %in% keys(plotTitlesHash)){
    catDataFrame <- data.frame(overallAttributes[c("catchment", "latitude", "longitude", currCat)])
    # Min/Midpoints/Max breaks; midpoint defined as middle between min and max
    minMidsMax <-  get_min_mids_max(catDataFrame[,currCat])
    # Round and then pad 0s, turn into strings
    minMidsMaxUntransformedLabels <- format_list_decimal_places(minMidsMax)          
     # Apply log transformation for visual clarity (if not si; si has negative values that do not log transform)
    if (currCat != "mean_sp" && currCat != "mean_si" && currCat != "bfi"){
        catDataFrame[, currCat] <- log(catDataFrame[, currCat])
        # Min/Midpoints/Max breaks; midpoint defined as middle between min and max
        minMidsMax <-  get_min_mids_max(catDataFrame[,currCat])
        # Untransform min, mid, max for labels
        minMidsMaxUntransformedLabels <- exp(minMidsMax)
        # Round and then pad 0s, turn into strings
        minMidsMaxUntransformedLabels <- format_list_decimal_places(minMidsMaxUntransformedLabels)             
    }
    
    # Long-Term Value Plot
    longTermPlot <- ggplot(data=northAmerica) +
     geom_sf()+ xlab("Longitude") + ylab("Latitude") +
      theme(legend.text =  element_text(size = 11),
            legend.key.height = unit(0.5, "cm"),
            legend.spacing.y = unit(0.25, "cm")
            ) +
      coord_sf(xlim = c(-145, -50), ylim = c(25, 72), expand = FALSE) +
      annotation_scale(location = "bl", width_hint = 0.15) +
      annotation_north_arrow(location = "bl", which_north = "true",
                             pad_x = unit(0.08, "in"), pad_y = unit(0.25, "in"),
                             style = north_arrow_fancy_orienteering,
                             height = unit(1, "cm"),
                             width = unit(1, "cm")) +
      geom_point(data = catDataFrame,
                 na.rm = TRUE,
                 shape = 21,
                 aes_string(x= "longitude",
                           y = "latitude",
                           fill = currCat
                           ),
                inherit.aes = FALSE) +
      scale_fill_gradient2(
                low = "#d7191c", 
                mid = "#ffffbf", 
                high = "#2c7bb6",
                midpoint = as.numeric(minMidsMax[ceiling(length(minMidsMax)/2)]),
                breaks = minMidsMax,
                labels = minMidsMaxUntransformedLabels,
      ) +         
      labs(title = plotTitlesHash[[currCat]],
           fill = NULL) +                  
      theme(legend.key.height = unit(0.85, 'cm'),
            plot.title = element_text(hjust = 0.5,
                                      ),
            plot.margin=grid::unit(c(0.01,0.01,0.01,0.01), "null"))
    longTermPlotFileName <- create_sig_filenames_string(destPlotDir = plotDir,
                                                      sigName = currCat,
                                                      prefixString = attributeGeoPlotFigures,
                                                      suffixString = "_long_term_plot.png")
    ggsave(filename = longTermPlotFileName,
      plot = longTermPlot,
      device= "png",
      dpi = 320)

    # Convert letter order to number 
    catPlotOrder <- which(letters == panelOrderHash[[currCat]])
    catPlotList[[catPlotOrder]] <- longTermPlot
  }
}

# Plot all panels as one 
allCatPlotsName <- paste(plotDir, attributeGeoPlotFigures, ".png", sep = "")
allCatPlots <- grid.arrange(grobs=catPlotList, ncol=2)
ggsave(allCatPlotsName,
          allCatPlots,
          height = 11,
          width = 8.5,
          dpi = 320)


