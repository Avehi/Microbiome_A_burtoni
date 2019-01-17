#############################################################################
#       Microbiome analysis for dominant and subordinate A. burtoni         #
#     Last updated 12/12/2018 by Josh Faber-Hammond and Avehi Singh         #
#             Original script written by Josh Faber-Hammond                 #
# See associated manuscript for more details about the analysis and results #
#############################################################################

# Load libraries
library(MASS)
library(plyr)
library(tidyr)
#install.packages("tidyr")
library(dplyr)
#install.packages("dplyr")
library(phangorn)
#install.packages("phangorn")
library(vegan)
#biocLite("vegan")
library("DESeq2")
#packageVersion("DESeq2")
#biocLite("DESeq2")
library(phyloseq)
#packageVersion("phyloseq")
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite(c("phyloseq"))
library(ggplot2)
library("scales")
library("grid")
library(tibble)
#install.packages('tibble', dependencies=TRUE, type="source")
library(FSA)
# install.packages("FSA")
library (labdsv)
# install.packages("labdsv")
library(magrittr)
library (rcompanion)
# install.packages("rcompanion")
library  (multcomp)
# install.packages("multcomp")
require(reshape2)
library(ape)
library(devtools)
#install.packages("devtools")
#source("http://bioconductor.org/biocLite.R")
#biocLite("ggtree")
library("ggtree")
#biocLite("ggtree")
library(plotrix)
#install.packages("plotrix")
#devtools::install_github('reptalex/phylofactor')
library(phylofactor)
require(phylofactor)
library(phytools)
#install.packages("phytools")
??PhyloFactor ### for help topics
#install.packages("inflection")
library(inflection)
#install.packages("inflection")
require(inflection)
#install.packages("chngpt")
library(chngpt)
library(Hmisc)
#install.packages("Hmisc")
library(mgcv)
#install.packages("picante")
library(picante)
#install.packages("abind")
library(abind)
#install.packages("gridExtra")
library(gridExtra)
#install.packages('lme4')
library(lme4)
#install.packages('lmerTest')
library(lmerTest)
#install.packages('MuMIn')
library(MuMIn)

theme_set(theme_bw())
set.seed(42)

#Load helper functions
ste <- function(x) sd(x)/sqrt(length(x))

calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

#' @title
#' Find the Elbow in a Curve
#'
#' @description
#' This utility function finds the elbow in a curve which is concave
#' relative to a line drawn between the first and last points.
#' The elbow is defined as the point with the greatest
#' orthogonal distance from that line.
#'
#' @param y Numeric vector of y values for the curve.
#' 
#' @param plot Logical. Should a plot be made?
#' 
#' @param returnIndex Logical. Should the return value
#' be the index of the elbow point? 
#' 
#' @return If \code{returnIndex = TRUE}, the index of
#' the elbow point.  If \code{returnIndex = FALSE},
#' a data frame containing an index values (x),
#' the y values passed to the function, and the
#' the orthogonal distances of the y values from
#' the line connecting the first and last points.
#' \code{which.max(data_frame_name$dist)} will give the index of 
#' the elbow point.  
#'
#' @references The concept of this function is based on the
#' clever idea in the
#' Stackoverflow post at stackoverflow.com/a/2022348/633251
#' and relies on code posted at
#' paulbourke.net/geometry/pointlineplane/pointline.r
#' (referenced in the SO post).  Minor modifications
#' to the code were made to that code in order to vectorize it.
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @export
#'
#' @importFrom stats lm coef
#'
#' @section Warning:
#' This function makes some simple checks that the data is concave
#' as defined above.  Even so, it may give
#' answers in some cases that are not valid.  Please check
#' on typical data that you encounter to verify that it works
#' in your cases.
#' 
#' @examples
#' tmp <- findElbow(c(0.9, 1.1, 1.1, 1.9, 2.5, 2.8, 4.9, 8.5),
#' 	plot = TRUE) # wandering
#' tmp <- findElbow(c(0.9, 1.0, 1.2, 1.3, 1.5, 1.5, 10.0, 22.0),
#' 	plot = TRUE) # late rise
#' tmp <- findElbow(c(2, 4, 6, 8, 10, 12, 14, 16)^2,
#' 	plot = TRUE) # gradual, no obvious break
#' 
#' # Not the usual way to choose the number of PCs:
#' library("chemometrics")
#' data(glass)
#' pca <- prcomp(glass)
#' eigensum <- sum(pca$sdev * pca$sdev)
#' vv <- 100 * (pca$sdev * pca$sdev/eigensum)
#' cs <- cumsum(vv)
#' tmp <- findElbow(vv, plot = TRUE)
#' tmp <- findElbow(cs, plot = TRUE)
#'

findElbow <- function(y, plot = FALSE, returnIndex = TRUE) {
  
  # The following helper functions were found at
  # paulbourke.net/geometry/pointlineplane/pointline.r
  # via the SO reference below.  The short segment check
  # was modified to permit vectorization.
  
  ##========================================================
  ##
  ##  Credits:
  ##  Theory by Paul Bourke http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
  ##  Based in part on C code by Damian Coventry Tuesday, 16 July 2002
  ##  Based on VBA code by Brandon Crosby 9-6-05 (2 dimensions)
  ##  With grateful thanks for answering our needs!
  ##  This is an R (http://www.r-project.org) implementation by Gregoire Thomas 7/11/08
  ##
  ##========================================================
  
  distancePointLine <- function(x, y, slope, intercept) {
    ## x, y is the point to test.
    ## slope, intercept is the line to check distance.
    ##
    ## Returns distance from the line.
    ##
    ## Returns 9999 on 0 denominator conditions.
    x1 <- x-10
    x2 <- x+10
    y1 <- x1*slope+intercept
    y2 <- x2*slope+intercept
    distancePointSegment(x,y, x1,y1, x2,y2)
  }
  
  distancePointSegment <- function(px, py, x1, y1, x2, y2) {
    ## px,py is the point to test.
    ## x1,y1,x2,y2 is the line to check distance.
    ##
    ## Returns distance from the line, or if the intersecting point on the line nearest
    ## the point tested is outside the endpoints of the line, the distance to the
    ## nearest endpoint.
    ##
    ## Returns 9999 on 0 denominator conditions.
    lineMagnitude <- function(x1, y1, x2, y2) sqrt((x2-x1)^2+(y2-y1)^2)
    ans <- NULL
    ix <- iy <- 0   # intersecting point
    lineMag <- lineMagnitude(x1, y1, x2, y2)
    if(any(lineMag < 0.00000001)) { # modified for vectorization by BAH
      #warning("short segment")
      #return(9999)
      warning("At least one line segment given by x1, y1, x2, y2 is very short.")
    }
    u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
    u <- u / (lineMag * lineMag)
    if(any((u < 0.00001) || (u > 1))) { # BAH added any to vectorize
      ## closest point does not fall within the line segment, take the shorter distance
      ## to an endpoint
      ix <- lineMagnitude(px, py, x1, y1)
      iy <- lineMagnitude(px, py, x2, y2)
      if(ix > iy)  ans <- iy
      else ans <- ix
    } else {
      ## Intersecting point is on the line, use the formula
      ix <- x1 + u * (x2 - x1)
      iy <- y1 + u * (y2 - y1)
      ans <- lineMagnitude(px, py, ix, iy)
    }
    ans
  }
  
  # End of helper functions by PB
  
  ### Now for the actual findElbow function!
  
  # Find the elbow using the method described in
  # stackoverflow.com/a/2022348/633251
  # but translated to R (see above).
  
  # Add an index to argument values for easy plotting
  DF <- data.frame(x = 1:length(y), y = y)
  fit <- lm(y ~ x, DF[c(1,nrow(DF)),]) # 2 point 'fit'
  m <- coef(fit)[2]
  b <- coef(fit)[1]
  
  # Check to make sure the data is concave as described
  # in the documentation, as arbitrary trends could give
  # misleading answers.  The following approach simply
  # checks to make sure all values are either above or
  # below the reference line.  This allows the values
  # to vary quite a bit and still return an answer.
  
  concave <- FALSE
  use <- 2:(nrow(DF)-1)
  refpts <- m*DF$x[use] + b
  if (all(refpts > DF$y[use]) | all(refpts < DF$y[use])) concave <- TRUE
  if (!concave) stop("Your curve doesn't appear to be concave")
  
  # Calculate the orthogonal distances
  use <- 2:(nrow(DF)-1)
  elbowd <- distancePointLine(DF$x[use], DF$y[use], coef(fit)[2], coef(fit)[1])
  DF$dist <- c(NA, elbowd, NA) # first & last points don't have a distance
  
  if (plot) {
    edm <- which.max(DF$dist)
    plot(DF[,1:2], type = "b", xlab = "index", ylab = "y values",
         main = "Looking for the Elbow")
    segments(DF$x[1], DF$y[1],
             DF$x[nrow(DF)], DF$y[nrow(DF)], col = "red")
    points(DF$x[edm], DF$y[edm], cex = 1.5, col = "red")	
    points(DF$x[edm], DF$y[edm], pch = 20)
  }
  
  if (returnIndex) return(which.max(DF$dist))
  if (!returnIndex) return(DF)
  
} # end of findElbow

colPercents <- function(tab, digits=1){
  dim <- length(dim(tab))
  if (is.null(dimnames(tab))){
    dims <- dim(tab)
    dimnames(tab) <- lapply(1:dim, function(i) 1:dims[i])
  }
  sums <- apply(tab, 2:dim, sum)
  per <- apply(tab, 1, function(x) x/sums)
  dim(per) <- dim(tab)[c(2:dim,1)]
  per <- aperm(per, c(dim, 1:(dim-1)))
  dimnames(per) <- dimnames(tab)
  per <- round(100*per, digits)
  result <- abind(per, Total=apply(per, 2:dim, sum), Count=sums, along=1)
  names(dimnames(result)) <- names(dimnames(tab))
  result
} # function from Rcmdr

#End loading helper functions

############################### Begin by importing input files and normalizing data for comparison of samples ################################

#############################
## Create Phyloseq object  ##
#############################

setwd("/Volumes/Bay_2/avehi/microbiome_final")
otufile <- read.table("AS_sample_filtered_OTU_table.txt", sep = "\t", header = TRUE, row.names=1)
## just import and then tell phyloseq that these are the files that it is looking for

#read in metadata, and calculate body-condition and testosterone measures corrected by length.
mapfile <- read.csv("AS_metadata_filt_newT_20181205.csv", header = TRUE, row.names=1)
t_byLeng <- mapfile$raw_T_pg.hr / (mapfile$length / 10)
bc.lm <- lm (log(weight)~length, data = mapfile)
bc.res <- resid(bc.lm)
as.matrix(bc.res)

# phyloseq needs rownames
meta_table<-sample_data(cbind(mapfile, t_byLeng, bc.res))
head(mapfile)
head(meta_table) ## they look the same

#Load taxa table
ep_tax<-read.table("AS_tax.txt", header=TRUE, sep = "\t", row.names=1) 
# phyloseq needs rownames
ep_tax <- as.matrix(ep_tax)  # coerce to matrix
ep_tax <- tax_table(ep_tax)

#Load phylogenetic tree
ep_tree<-read.tree("AS_phylogeny.tre")  # must be a rooted tree for phyloseq
is.rooted(phy_tree(ep_tree)) 
ep_tree_rooted <-midpoint(ep_tree) # this is from phanghorn package that sets root to midpoint
is.rooted(phy_tree(ep_tree_rooted)) 

###########################
# For unrarified raw data #
###########################

otu <- otu_table(otufile, taxa_are_rows = TRUE)  # transpose back just to use commands
phyloseq_data_norarify <- phyloseq(otu, meta_table, ep_tax, ep_tree_rooted)

# Subset to filter out Cyanobacteria phylum.
phyloseq.noCy <- subset_taxa(phyloseq_data_norarify, Phylum!="Cyanobacteria")
# Remove the samples that have less than 3k total reads from filtered data
phyloseq_data_norarify <- prune_samples(sample_sums(phyloseq.noCy)>=3000, phyloseq.noCy)

#########################
# For rarified raw data #
#########################
otu_filt <- otu_table(phyloseq_data_norarify)
#rarefy data! This is optional and related to normalization when different samples have different read depths 
#these argumenents mattered more for 454
# subsample to minimum count for samples passing filters.
# make this number match your read count/filtered datasets
otur <- rrarefy(t(otu_filt),3879)  # from vegan package (set number for dataset) (transpose)
rowSums(otur)
colSums(otu_filt)
colSums(otufile)

otu <- otu_table(t(otur), taxa_are_rows = TRUE)  # transpose back just to use commands
#write.table(otu, file = "rarified_otu_table.txt", sep = "\t")

#lastly, turn all your inputs into phyloseq object
phyloseq_data <- phyloseq(otu, meta_table, ep_tax, ep_tree_rooted)

##################################################
#     Variance stabilizing normalization and     #
# differential abundance analysis through DESeq2 #
##################################################

# DESeq2
# *** MAKE SURE LOADED BEFORE PHYLOSEQ ***
phyloseq_data_ds <- subset_samples(phyloseq_data_norarify, prev_type!="NA")

#Conversion to DESeq data
phyloseq_ds = phyloseq_to_deseq2(phyloseq_data_ds, ~ prev_type)
#Run negative bionomial Wald test to test for differential abundance of specific taxa
diagdds = DESeq(phyloseq_ds, test="Wald", fitType="parametric", betaPrior=FALSE)
#Or run LRT to use more complex models
#diagdds = DESeq(phyloseq_ds, test="LRT", reduced = ~ fish)

# To extract results, choose the variable of interest and the pair of groups you want to test
# Cooks Cutoff dictates whether or not you want to remove outliers
res <- results(diagdds, cooksCutoff = TRUE )
res <- results(diagdds, contrast = c("prev_type", "D", "S"), cooksCutoff = TRUE )

# Filter for your desired alpha value, and extract relevant taxon names
alpha = 0.05
sigtab = res[which(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_data_norarify)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

# List results
sigtab

#### Plot differentially abundant taxa from pair of groups ####
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
sigtab$names <- rownames(sigtab)
ggplot(sigtab, aes(x=reorder(names, as.numeric(Genus)), y=log2FoldChange, color=Genus)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(x = "OTU", y = 'Log2 fold change')

#####################################################
# Convert DESeq2 table to usable table for phyloseq #
# then substitute into the previous phyloseq object #
#####################################################

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)
dim(diagvst)

# duplicate old phyloseq object for backup
phyloseq_data_ds2 <- phyloseq_ds

# Recommended to run the next command to floor data at 0 for input 
# into ordination or other downstream analyses
diagvst[diagvst < 0.0] <- 0.0

# Insert DESeq2 normalized data into phyloseq object
phyloseq_data_ds2_count <- otu_table(diagvst, taxa_are_rows = TRUE)
phyloseq_data_ds2 <- phyloseq(phyloseq_data_ds2_count, meta_table, ep_tax, ep_tree_rooted)

#optional output from DESeq2 normalized data
#write.table(diagvst, file = "DESeq2_table_correlational_analysis_floored.txt", sep = "\t")

####################################################################
# Subset the dataset for downstream analyses by metadata variables #
####################################################################
#subsetting used for exploratory analysis
#use this subset for correlational analysis
correlation_subset <- subset_samples(phyloseq_data_ds2, outlier=="ok")
correlation_subset2 <- subset_samples(correlation_subset, filter=="y")

#use this subset for correlational analysis
correlation_subset <- subset_samples(phyloseq_data, outlier=="ok")
correlation_subset2 <- subset_samples(correlation_subset, filter=="y")

#use this subset for correlational analysis
correlation_subset <- subset_samples(phyloseq_data_norarify, outlier=="ok")
correlation_subset2 <- subset_samples(correlation_subset, filter=="y")

#define subset to test
microbiome_subset <- phyloseq_data
microbiome_subset <- phyloseq_data_norarify
microbiome_subset <- phyloseq_data_ds2  #This is the subset we used for beta diversity analysis. Do not run with bray or jaccard distance methods unless counts floored at 0
microbiome_subset <- correlation_subset
microbiome_subset <- correlation_subset2

################################################## BETA DIVERSITY #######################################################

# Ordinations! try different methods and distance metrics
# Do not use OTU tables with negative values
# Interpret unifrac and jaccard ordinations with caution if using DESeq2 normalized count tables
# These two distance methods are presence/absence based and inappropriate for VST datasets.
ordu <- ordinate(microbiome_subset,"NMDS", "bray")
ordu2 <- ordinate(microbiome_subset,"NMDS", "unifrac")
ordu3 <- ordinate(microbiome_subset,"NMDS", "wunifrac")
ordu4 <- ordinate(microbiome_subset,"NMDS", "jaccard")
ordu5 <- ordinate(microbiome_subset,"PCoA", "bray")
ordu6 <- ordinate(microbiome_subset,"PCoA", "unifrac")
ordu7 <- ordinate(microbiome_subset,"PCoA", "wunifrac")
ordu8 <- ordinate(microbiome_subset,"PCoA", "jaccard")
ordu9 <- ordinate(microbiome_subset,"DPCoA")

#some extra plots
plot_ordination(microbiome_subset, ordu5, type = "scree")
stressplot(ordu)

####################################################
#         Exploratory beta diversity plots         #
####################################################
#Try looking at previous phenotype using a variety of distance methods
plot_ordination(microbiome_subset, ordu, color="prev_type", title=paste("Bray-curtis NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu2, color="prev_type", title=paste("UniFrac NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu3, color="prev_type", title=paste("WuniFrac NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu4, color="prev_type", title=paste("Jaccard NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu5, color="prev_type", title=paste("Bray-curtis PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu6, color="prev_type", title=paste("UniFrac PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu7, color="prev_type", title=paste("WuniFrac PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu8, color="prev_type", title=paste("Jaccard PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu9, color="prev_type", title=paste("DPCoA")) + geom_point(size=4) + coord_fixed()

# Check plotting of different variables
plot_ordination(microbiome_subset, ordu2, color="type", title=paste("UniFrac NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu3, color="type", title=paste("WuniFrac NMDS")) + geom_point(size=4) + coord_fixed()

#Check effect of sequencing depth
#Do this after adding a read count column labeled "count" to metadata table
plot_ordination(microbiome_subset, ordu, color="count", shape = "group",title=paste("Bray-curtis NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu2, color="count", shape = "group",title=paste("UniFrac NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu3, color="count", shape = "group",title=paste("WuniFrac NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu4, color="count", shape = "group",title=paste("Jaccard NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu5, color="count", shape = "group",title=paste("Bray-curtis PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu6, color="count", shape = "group",title=paste("UniFrac PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu7, color="count", shape = "group",title=paste("WuniFrac PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu8, color="count", shape = "group",title=paste("Jaccard PCoA")) + geom_point(size=4) + coord_fixed()

plot_ordination(microbiome_subset, ordu, color="group", shape = "prev_type", title=paste("Bray-curtis NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu2, color="group", shape = "prev_type", title=paste("UniFrac NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu3, color="group", shape = "prev_type",  title=paste("WuniFrac NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu4, color="fish", title=paste("Jaccard NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu5, color="fish", title=paste("Bray-curtis PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu6, color="fish", title=paste("UniFrac PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu7, color="fish", title=paste("WuniFrac PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu8, color="fish", title=paste("Jaccard PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microbiome_subset, ordu9, color="fish", title=paste("DPCoA")) + geom_point(size=4) + coord_fixed()


#######################################
#         Plots in Manuscript         #
#######################################
plot_ordination(microbiome_subset, ordu3, color="prev_type", shape = "type", title=paste("WuniFrac NMDS")) + geom_point(size=5) + coord_fixed() + 
    scale_color_manual(values = c('darkgoldenrod2', 'turquoise4')) + 
    scale_shape_manual(values=c(19,1))

############################################
#beta diversity statistsical tests
#CHOOSE ONE OF THESE SUBSETS AT A TIME TO RUN WITH DOWNSTREAM ANALYSES

#use this subset for VST normalized data with no "NA" values in metadata
#This is the one we recommend
correlation_subset <- subset_samples(phyloseq_data_ds2, outlier=="ok")
correlation_subset2 <- subset_samples(correlation_subset, t_byLeng > "0")
microbiome_subset <- correlation_subset2

#use this subset for rarefied data with no "NA" values in metadata
correlation_subset <- subset_samples(phyloseq_data, outlier=="ok")
correlation_subset2 <- subset_samples(correlation_subset, filter=="y")
microbiome_subset <- correlation_subset2

#use this subset for unrarefied data with no "NA" values in metadata
#(not an appropriate subset for between sample comparison but just an exploratory analysis)
correlation_subset <- subset_samples(phyloseq_data_norarify, outlier=="ok")
correlation_subset2 <- subset_samples(correlation_subset, filter=="y")
microbiome_subset <- correlation_subset2

#or define a full dataset to test
microbiome_subset <- phyloseq_data
microbiome_subset <- phyloseq_data_norarify
microbiome_subset <- phyloseq_data_ds2  #This is the subset we used for beta diversity analysis. Do not run with bray or jaccard distance methods unless counts floored at 0


#We tested weighted-UniFrac and Bray-Curtis, so change distance method for each
metaD = as(sample_data(microbiome_subset), "data.frame")
#use either one of the following
phydist = phyloseq::distance(microbiome_subset, "wunifrac")
phydist = phyloseq::distance(microbiome_subset, "bray")

#Run a Goodness-of-fit analysis to determine which variables best explain data in multi-dimensional space
#We used significant variables identified by this analysis for follow-up analysis with Adonis
metaD <- cbind.data.frame(row.names = row.names(metaD),prev_type = metaD$prev_type,type = metaD$type,timepoint = metaD$timepoint,tank = metaD$tank, fish = metaD$fish, paired_with = metaD$paired_with)
ef <- envfit(phydist, metaD, permu = 999, na.rm = TRUE)
ef

#Reset metaD to have all exploratory variables
metaD = as(sample_data(microbiome_subset), "data.frame")

#Prior to running Adonis, check to see whether the groups have roughly equal dispersion
#This is an assumption of PERMANOVA, and you need to report if data doesn't meet this assumption
beta <- betadisper(phydist, metaD$prev_type)
permutest(beta)

beta <- betadisper(phydist, metaD$type)
permutest(beta)

#PERMANOVA by different variable groups
#This tests for differences in centroids in multidimensional matrix
#Try different models using goodness-of-fit analysis as a guide for model selection
#Although "fish" (individual ID) wasn't found to be significant, include it anyway to account for repeated measures of the same fish
#Ideally, fish could be used as a random effect variable, but Adonis is not compatible with random effect models.
#Also, try adonis with current type (although not significant) to test our original hypothesis and confirm revised hypothesis better fits data
adonis_res <- adonis2(formula = phydist ~ metaD$prev_type, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$type, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$prev_type + metaD$fish, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$type + metaD$fish, data = metaD) 
adonis_res

#check other variables/models to confirm they are non-significant
adonis_res <- adonis2(formula = phydist ~ metaD$group, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$bc.res, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$timepoint, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$group:metaD$timepoint, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$group + metaD$timepoint + metaD$group:metaD$timepoint, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$behavior_scores, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$tank, data = metaD) 
adonis_res

################################################### alpha diversity ############################################################
# Initialize matrices to store richness and evenness estimates
# Make sure subset is your unrarified dataset
meta_table2 <- data.frame(sample_data(phyloseq_data_norarify))
min_lib <- min(sample_sums(phyloseq_data_norarify))
nsamp = nsamples(phyloseq_data_norarify)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(phyloseq_data_norarify)

# It is always important to set a seed when you subsample so your result is replicable 
set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(phyloseq_data_norarify, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
}

# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
SR_ave <- apply(richness, 1, mean)
SR_se <- apply(richness, 1, ste)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, SR_ave, SR_se)

#alpha diversity
#test for cvorrelations between differentially abundant OTUs, alpha diversity, and metadata variables
richness <- estimate_richness(phyloseq_data_norarify, split = TRUE, measures = NULL)
otufile_norar <- t(data.frame(otu_table(phyloseq_data_norarify)))
faith <- pd(otufile_norar, ep_tree, include.root=FALSE)
faith <- cbind(faith, prev_type = meta_table2$prev_type, type = meta_table2$type) 

# create table for correlational analysis
alphadiv_meta <- cbind(SR_ave, SR_se, richness, Faith_norar = faith$PD, behavior_scores = meta_table2$behavior_scores, weight = meta_table2$weight, length = meta_table2$length, bc.res = meta_table2$bc.res, t_byLeng = meta_table2$t_byLeng, otu4_deseq = meta_table2$deseq_norm_OTU_4, otu6_deseq = meta_table2$deseq_norm_OTU_6, otu21_deseq = meta_table2$deseq_norm_OTU_21, otu143_deseq = meta_table2$deseq_norm_OTU_143, otu474_deseq = meta_table2$deseq_norm_OTU_474, otu438_deseq = meta_table2$deseq_norm_OTU_438, otu376_deseq = meta_table2$deseq_norm_OTU_376, otu155_deseq = meta_table2$deseq_norm_OTU_155, otu125_deseq = meta_table2$deseq_norm_OTU_125, otu4_rar = meta_table2$OTU_4, otu6_rar = meta_table2$OTU_6, otu21_rar = meta_table2$OTU_21, otu143_rar = meta_table2$OTU_143, otu474_rar = meta_table2$OTU_474, otu438_rar = meta_table2$OTU_438, otu376_rar = meta_table2$OTU_376, otu155_rar = meta_table2$OTU_155, otu125_rar = meta_table2$OTU_125)
alphadiv_meta
nrow(alphadiv_meta)
# write table with diversity stats and metadata
richness_meta <- cbind(SR_ave, SR_se, richness, Faith_norar = faith$PD, meta_table2)
# write.table(richness_meta, file = "richness+meta.txt", sep = "\t")


#to actually test correlations between all fields in the diversity table, make sure your dataset has no "NA" values
#filter samples by subset function prior to correlation check
alphadiv_meta_mat <- as.matrix(alphadiv_meta)
alphacorr <- rcorr(alphadiv_meta_mat)
alphacorr

# write tables with correlation results
#write.table(alphacorr$P, file = "alpha_div_correlations_wOTUs_Pval.txt", sep = "\t")
#write.table(alphacorr$r, file = "alpha_div_correlations_wOTUs_rval.txt", sep = "\t")

# Make plots of correlations
# First, define dataset to plot by calling it alpha_subset. Unrarified better for most measures, but you may
# want to substiture rarified to better compare samples (since they are sampled at same depth).
alpha_subset <- phyloseq_data_norarify

#plot alpha diversity for Previous Type
p = plot_richness(alpha_subset, x = "prev_type",
                  title = NULL, scales = "free_y", nrow = 1, shsi = NULL,
                  measures = c("Observed", "Shannon"), sortby = NULL)
p + geom_boxplot(data = p$alpha_subset, aes(x = prev_type, y = value, color = NULL), alpha = 0.1) + labs(x = "Previous Phenotype")

#plot alpha diversity for group:timepoint
#optional filter for timepoint
#alpha_subset <- subset_samples(alpha_subset, timepoint=="second")
p = plot_richness(alpha_subset, x = "group:timepoint",
                  title = NULL, scales = "free_y", nrow = 1, shsi = NULL,
                  measures = c("Observed", "Shannon"), sortby = NULL)
p + geom_boxplot(data = p$alpha_subset, aes(x = group:timepoint, y = value, color = NULL), alpha = 0.1) + labs(x = "Group:Timepoint")


#######################################
#         Plots in Manuscript         #
#######################################

#make boxplot with colours!
#previous type plot
p = plot_richness(alpha_subset, x = "prev_type",
                  title = NULL, scales = "free_y", nrow = 1, shsi = NULL,
                  measures = c("Observed", "Shannon"), sortby = NULL)
p + geom_boxplot(data = p$alpha_subset, aes(x = prev_type, y = value, fill = prev_type)) + 
  labs(x = "Previous Phenotype")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))

#current type plot
p = plot_richness(alpha_subset, x = "type",
                  title = NULL, scales = "free_y", nrow = 1, shsi = NULL,
                  measures = c("Observed", "Shannon"), sortby = NULL)
p + geom_boxplot(data = p$alpha_subset, aes(x = type, y = value, fill = type)) + 
  labs(x = "Current Phenotype")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))

#Plot faith's PD (and SR calculated through Picante)
faith_plot <- faith[!is.na(faith$prev_type),]

#previous type
p = ggplot(faith_plot, aes(x = prev_type, y = PD, group = prev_type)) 
p + geom_boxplot(data = faith_plot, aes(x = prev_type, fill = prev_type)) + 
  labs(x = "Previous Phenotype")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))+
  geom_point()

#current type
p = ggplot(faith_plot, aes(x = type, y = PD, group = type)) 
p + geom_boxplot(data = faith_plot, aes(x = type, fill = type)) + 
  labs(x = "Current Phenotype")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))+
  geom_point()

#Plot SR and Shannon for all non-NA prev_type fish
richness_plot <- richness_meta[!is.na(richness_meta$prev_type),]

#Plot both current type and prev type on the same graph for alpha diversity stats
#Species Richness
p = ggplot(richness_plot, aes(x = prev_type:type, y = SR_ave, group = prev_type:type)) 
p1 = p + geom_boxplot(data = richness_plot, aes(x = prev_type:type, fill = prev_type:type)) + 
  labs(x = "Previous:Current Phenotype")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'darkgoldenrod2','turquoise4', 'turquoise4'))+
  geom_point() + 
  theme(legend.position="none")
#Shannon Index
p = ggplot(richness_plot, aes(x = prev_type:type, y = Shannon, group = prev_type:type)) 
p2 = p + geom_boxplot(data = richness_plot, aes(x = prev_type:type, fill = prev_type:type)) + 
  labs(x = "Previous:Current Phenotype", y = "Shannon Index")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'darkgoldenrod2','turquoise4', 'turquoise4'))+
  geom_point() + 
  theme(legend.position="none")
#Faith's Phylogenetic Diversity
p = ggplot(richness_plot, aes(x = prev_type:type, y = Faith_norar, group = prev_type:type)) 
p3 = p + geom_boxplot(data = richness_plot, aes(x = prev_type:type, fill = prev_type:type)) + 
  labs(x = "Previous:Current Phenotype", y = "Faith's PD")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'darkgoldenrod2','turquoise4', 'turquoise4'))+
  geom_point() + 
  theme(legend.position="none")

grid.arrange(p1, p2, p3, nrow = 1)

# Comparison of Species Richness based on 100 iterations of rarefication (SR_ave) and raw counts (Observed)
p = ggplot(richness_plot, aes(x = prev_type:type, y = Observed, group = prev_type:type)) 
p4 = p + geom_boxplot(data = richness_plot, aes(x = prev_type:type, fill = prev_type:type)) + 
  labs(x = "Previous:Current Phenotype")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'darkgoldenrod2','turquoise4', 'turquoise4'))+
  geom_point() + 
  theme(legend.position="none")

grid.arrange(p1, p4, nrow = 1)

##############################################################
#         Calculate Alpha-Div for S and D separately         #
##############################################################

#subset table to separate subordinate and dominant and test for OTU abundances and test for correlations
alphadiv_for_subset <- cbind(SR_ave, SR_se, richness, Faith_norar = faith$PD, prev_type = meta_table2$prev_type, behavior_scores = meta_table2$behavior_scores, weight = meta_table2$weight, length = meta_table2$length, bc.res = meta_table2$bc.res, t_byLeng = meta_table2$t_byLeng, otu4_deseq = meta_table2$deseq_norm_OTU_4, otu6_deseq = meta_table2$deseq_norm_OTU_6, otu21_deseq = meta_table2$deseq_norm_OTU_21, otu143_deseq = meta_table2$deseq_norm_OTU_143, otu474_deseq = meta_table2$deseq_norm_OTU_474, otu438_deseq = meta_table2$deseq_norm_OTU_438, otu376_deseq = meta_table2$deseq_norm_OTU_376, otu155_deseq = meta_table2$deseq_norm_OTU_155, otu125_deseq = meta_table2$deseq_norm_OTU_125, otu4_rar = meta_table2$OTU_4, otu6_rar = meta_table2$OTU_6, otu21_rar = meta_table2$OTU_21, otu143_rar = meta_table2$OTU_143, otu474_rar = meta_table2$OTU_474, otu438_rar = meta_table2$OTU_438, otu376_rar = meta_table2$OTU_376, otu155_rar = meta_table2$OTU_155, otu125_rar = meta_table2$OTU_125)
subordinate <- subset(alphadiv_for_subset, prev_type == 'S')
dominant <- subset(alphadiv_for_subset, prev_type == 'D')

#make correlations for subordinate and dominant separately

subordinate <- subset(subordinate, select = -c(prev_type) )
sub_alphacorr <- rcorr(as.matrix(subordinate))
sub_alphacorr
#write.table(sub_alphacorr$P, file = "alpha_div_correlations_subordinates_Pval.txt", sep = "\t")
#write.table(sub_alphacorr$r, file = "alpha_div_correlations_subordinates_rval.txt", sep = "\t")

dominant <- subset(dominant, select = -c(prev_type) )
dom_alphacorr <- rcorr(as.matrix(dominant))
dom_alphacorr
#write.table(dom_alphacorr$P, file = "alpha_div_correlations_dominants_Pval.txt", sep = "\t")
#write.table(dom_alphacorr$r, file = "alpha_div_correlations_dominants_rval.txt", sep = "\t")

#par(mfrow=c(1,2))
richness_manuscript <- cbind(richness_meta, prev_type_col = meta_table2$prev_type)
richness_ms <- as.matrix(richness_manuscript)
richness_ms <- richness_manuscript

# colour by prev_type 
richness_ms$prev_type_col <- as.character(richness_ms$prev_type_col)
richness_ms$prev_type_col <- replace(richness_ms$prev_type_col, richness_ms$prev_type_col =="D", "darkgoldenrod2")
richness_ms$prev_type_col <- replace(richness_ms$prev_type_col, richness_ms$prev_type_col =="S", "turquoise4")

par(mfrow=c(2,3))

# Plot correlations of OTUs 4 and 6 with Species Richness, Shannon Index, and Faith's PD using newly generated colours
#scatter plots of alpha diversity vs. differentially abundant OTU counts
plot(richness_ms$SR_ave,richness_ms$deseq_norm_OTU_4, xlab = "Species Richness", ylab = "OTU_4 normalised abundance", col=richness_ms$prev_type_col, pch = 19, cex = 1.5)
abline(lm(richness_ms$deseq_norm_OTU_4~richness_ms$SR_ave), col="black") # regression line (y~x) 
abline(lm(subordinate$otu4_deseq~subordinate$SR_ave), col="turquoise4") # regression line (y~x) 
abline(lm(dominant$otu4_deseq~dominant$SR_ave), col="darkgoldenrod2") # regression line (y~x) 
#lines(lowess(richness_ms$deseq_norm_OTU_4~richness_ms$SR_ave), col="blue") # lowess line (x,y)

plot(richness_ms$Shannon,richness_ms$deseq_norm_OTU_4, xlab = "Shannon Index", ylab = "OTU_4 normalised abundance", col=richness_ms$prev_type_col, pch = 19, cex = 1.5)
abline(lm(richness_ms$deseq_norm_OTU_4~richness_ms$Shannon), col="black") # regression line (y~x) 
abline(lm(subordinate$otu4_deseq~subordinate$Shannon), col="turquoise4") # regression line (y~x) 
abline(lm(dominant$otu4_deseq~dominant$Shannon), col="darkgoldenrod2") # regression line (y~x) 
#lines(lowess(richness_ms$deseq_norm_OTU_4~richness_ms$Shannon), col="blue") # lowess line (x,y)

plot(richness_ms$Faith_norar,richness_ms$deseq_norm_OTU_4, xlab = "Faith's PD", ylab = "OTU_4 normalised abundance", col=richness_ms$prev_type_col, pch = 19, cex = 1.5)
abline(lm(richness_ms$deseq_norm_OTU_4~richness_ms$Faith_norar), col="black") # regression line (y~x) 
abline(lm(subordinate$otu4_deseq~subordinate$Faith_norar), col="turquoise4") # regression line (y~x) 
abline(lm(dominant$otu4_deseq~dominant$Faith_norar), col="darkgoldenrod2") # regression line (y~x) 
#lines(lowess(richness_ms$deseq_norm_OTU_4~richness_ms$Faith_norar), col="blue") # lowess line (x,y)

plot(richness_ms$SR_ave,richness_ms$deseq_norm_OTU_6, xlab = "Species Richness", ylab = "OTU_6 normalised abundance", col=richness_ms$prev_type_col, pch = 19, cex = 1.5)
abline(lm(richness_ms$deseq_norm_OTU_6~richness_ms$SR_ave), col="black") # regression line (y~x) 
abline(lm(subordinate$otu6_deseq~subordinate$SR_ave), col="turquoise4") # regression line (y~x) 
abline(lm(dominant$otu6_deseq~dominant$SR_ave), col="darkgoldenrod2") # regression line (y~x) 
#lines(lowess(richness_ms$deseq_norm_OTU_6~richness_ms$SR_ave), col="blue") # lowess line (x,y)

plot(richness_ms$Shannon,richness_ms$deseq_norm_OTU_6, xlab = "Shannon Index", ylab = "OTU_6 normalised abundance", col=richness_ms$prev_type_col, pch = 19, cex = 1.5)
abline(lm(richness_ms$deseq_norm_OTU_6~richness_ms$Shannon), col="black") # regression line (y~x) 
abline(lm(subordinate$otu6_deseq~subordinate$Shannon), col="turquoise4") # regression line (y~x) 
abline(lm(dominant$otu6_deseq~dominant$Shannon), col="darkgoldenrod2") # regression line (y~x) 
#lines(lowess(richness_ms$deseq_norm_OTU_6~richness_ms$Shannon), col="blue") # lowess line (x,y)

plot(richness_ms$Faith_norar,richness_ms$deseq_norm_OTU_6, xlab = "Faith's PD", ylab = "OTU_6 normalised abundance", col=richness_ms$prev_type_col, pch = 19, cex = 1.5)
abline(lm(richness_ms$deseq_norm_OTU_6~richness_ms$Faith_norar), col="black") # regression line (y~x) 
abline(lm(subordinate$otu6_deseq~subordinate$Faith_norar), col="turquoise4") # regression line (y~x) 
abline(lm(dominant$otu6_deseq~dominant$Faith_norar), col="darkgoldenrod2") # regression line (y~x) 
#lines(lowess(richness_ms$deseq_norm_OTU_6~richness_ms$Faith_norar), col="blue") # lowess line (x,y)

par(mfrow=c(1,1))


###########################################
# Test whether alpha diversity and BC/GSI #
#  stats differ among treatment groups    #
###########################################
alphadiv <- cbind(richness_plot,length_mm = richness_plot$length/10)
dim(alphadiv)
as.matrix(colnames(alphadiv))  # just to see colnames

plot(log(weight)~length_mm, col = ifelse(prev_type=="S","turquoise4", "darkgoldenrod2"), pch = 19, data = alphadiv) ## looks linear
bc.lm <- lm(log(weight)~length_mm, data = alphadiv)
plot(log(alphadiv$weight)~alphadiv$length_mm, col = alphadiv$type )
abline(lm(log(weight)~length_mm, data = alphadiv ))
bc.res <- resid(bc.lm)  # these are actuall in alphadiv already

bc.lmer<-lmer(alphadiv$bc.res~alphadiv$type+(1|alphadiv$fish));summary(bc.lmer) 
# obs = 34, t = -1.493, P =  0.145 # stuff I did for paper had 35 fish must be one that died ; will change text to these stats

alphadiv.t <- alphadiv[!is.na(alphadiv$t_byLeng),]  # need to get rid of the two NAs, this column is raw T/ length mm
plot(alphadiv.t$t_byLeng ~ alphadiv.t$length_mm, col = ifelse(alphadiv.t$type=="S","turquoise4", "darkgoldenrod2"), pch = 19)
plot(alphadiv.t$t_byLeng ~ alphadiv.t$type, col = c("darkgoldenrod2", "turquoise4"))

# Test experiemtal design: Current dominant and subordinate fish have significantly different testosterone levels
t.lmer<-lmer(alphadiv.t$t_byLeng~alphadiv.t$type+(1|alphadiv.t$fish));summary(t.lmer) 


#MODELS TO INCLUDE IN RESUBMISSION START HERE
#first need to remove that one fish that doesn't have a prev_types observation:
alphadiv<-alphadiv%>%filter(prev_type!="NA")

#now SR_ave:
#first normalize the residuals if necessary:
SR_ave_lm<-lm(data=alphadiv, SR_ave~prev_type*type);plot(resid(SR_ave_lm));shapiro.test(resid(SR_ave_lm));hist(resid(SR_ave_lm))
boxcox(SR_ave_lm,lambda=seq(-1,1,0.1))
SR_ave_lm<-lm(data=alphadiv, SR_ave^(-0.25)~prev_type*type);plot(resid(SR_ave_lm));shapiro.test(resid(SR_ave_lm));hist(resid(SR_ave_lm))
SR_ave_lmer<-lmer(data=alphadiv, SR_ave^(-0.25)~prev_type*type+(1|fish)+(1|paired_with));summary(SR_ave_lmer)
#employ model selection with MuMIn:
SR_ave_lm_global<-lm(data=alphadiv, SR_ave^(-0.25)~prev_type*type)
options(na.action="na.fail");best.subsets<-dredge(SR_ave_lm_global,rank="AICc",trace=TRUE);print(best.subsets, abbrev.names=FALSE);options(na.action="na.omit")
#one top-fitting model only:
get.models(best.subsets,subset=1:1)
SR_ave_best_model<-lmer(data=alphadiv,SR_ave^(-0.25)~prev_type+(1|fish)+(1|paired_with));summary(SR_ave_best_model)
#prev_type S has a negative effect on SR, p = 0.0177  

#now Shannon:
#first normalize the residuals if necessary:
Shannon_lm<-lm(data=alphadiv, Shannon~prev_type*type);plot(resid(Shannon_lm));shapiro.test(resid(Shannon_lm));hist(resid(Shannon_lm))
#not necessary
#employ model selection with MuMIn:
Shannon_lm_global<-lm(data=alphadiv, Shannon~prev_type*type)
options(na.action="na.fail");best.subsets<-dredge(Shannon_lm_global,rank="AICc",trace=TRUE);print(best.subsets, abbrev.names=FALSE);options(na.action="na.omit")
#two top-fitting models:
get.models(best.subsets,subset=1:2)
#prev_type and the null model are the two top-fitters - prev_type a little better AICc than the null!
fmList<-get.models(best.subsets,subset=1:2); summary(model.avg(fmList))
#conditional model avg gives a p of 0.044 for the negative effect of prev_typeS on Shannon's
#conditional model avg uses only those two top-fitted models, so this seems appropriate
#full model avg gives only p = 0.238, but incorporates all models (weighted by position in heirarchy)
Shannon_best_model<-lmer(data=alphadiv,Shannon~prev_type+(1|fish)+(1|paired_with));summary(Shannon_best_model)
#prev_type S has a negative effect on Shannon, p = 0.0416

#now Faith_norar:
#first normalize the residuals if necessary:
Faith_norar_lm<-lm(data=alphadiv,Faith_norar~prev_type*type);plot(resid(Faith_norar_lm));shapiro.test(resid(Faith_norar_lm));hist(resid(Faith_norar_lm))
#not necessary
#employ model selection with MuMIn:
Faith_norar_lm_global<-lm(data=alphadiv, Faith_norar~prev_type*type)
options(na.action="na.fail");best.subsets<-dredge(Faith_norar_lm_global,rank="AICc",trace=TRUE);print(best.subsets, abbrev.names=FALSE);options(na.action="na.omit")
#two top-fitting models:
get.models(best.subsets,subset=1:2)
#prev_type and the null model are the two top-fitters - prev_type a little better AICc than the null!
fmList<-get.models(best.subsets,subset=1:2); summary(model.avg(fmList))
#conditional model avg gives a p of 0.15 for the negative effect of prev_typeS on Faith's
#conditional model avg uses only those two top-fitted models, so this seems appropriate
#full model avg gives only p = 0.485, but incorporates all models (weighted by position in heirarchy)
Faith_norar_best_model<-lmer(data=alphadiv,Faith_norar~prev_type+(1|fish)+(1|paired_with));summary(Faith_norar_best_model)
#prev_type S has a non-sig negative effect on Faith's, p = 0.185

#now Chao1:
# Although this statistic was not used in Manuscript, it is a good measure of species richness that can be compared among unrarified samples with different sampling depths.
#first normalize the residuals if necessary:
Chao1_lm<-lm(data=alphadiv, Chao1~prev_type*type);plot(resid(Chao1_lm));shapiro.test(resid(Chao1_lm));hist(resid(Chao1_lm))
#not necessary
#employ model selection with MuMIn:
Chao1_lm_global<-lm(data=alphadiv, Chao1~prev_type*type)
options(na.action="na.fail");best.subsets<-dredge(Chao1_lm_global,rank="AICc",trace=TRUE);print(best.subsets, abbrev.names=FALSE);options(na.action="na.omit")
#three top-fitting models:
get.models(best.subsets,subset=1:3)
#huh, type makes an appearance this time
fmList<-get.models(best.subsets,subset=1:3); summary(model.avg(fmList))
#conditional model avg gives a p of 0.0777 for the negative effect of prev_typeS on Chao1
#type S has a neg effect on Chao1 with p = 0.3985, then down to p = 0.762 in the full average - seems discardable
#conditional model avg uses only those two top-fitted models, so this seems appropriate
#full model avg gives only p = 0.295, but incorporates all models (weighted by position in heirarchy)
Chao1_best_model<-lmer(data=alphadiv,Chao1~prev_type+type+(1|fish)+(1|paired_with));summary(Chao1_best_model)
#prev_type S has a non-sig negative effect on Chao1, p = 0.0777



#############################################
#         Physiological comparisons         #
#############################################
data <- read.table(file = "AS_metadata_filt_newT_20181205.csv", header = TRUE, sep = "," )
plot(behavior_scores~type, data = data)  #### Make this one in ggplot
pbs = ggplot(data, aes(x = type, y = behavior_scores, group = type)) 
pbs + geom_boxplot(data = data, aes(x = type, fill = type)) + 
  labs(x = "Current Phenotype", y = "Behavior Score")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))+
  geom_point()

summary(aov(behavior_scores~type, data = data))  # df = 1,33, F = 330.1 p<2e-16
wilcox.test(behavior_scores~type, data = data, paired = FALSE )

#plot(weight~length, data = data)  ## roughly linear, log linear?
#plot(log(weight)~length, col = ifelse(prev_type=="S","cyan", "gold"), pch = 19, data = data) ## looks linear

bc.lm <- lm (log(weight)~length, data = data)
bc.res <- resid(bc.lm)
as.matrix(bc.res)  # These are the BC values that Josh should use for correlations

#### they are added to the file "AS_metadata_filt_newT_20181205.csv"
# cor.test (bc.res, data$behavior_scores, method = "pearson") # R2 = 0.311125, t = 1.8806, df = 33, p-value = 0.06887
cor.test (bc.res, data$behavior_scores, method = "spearman") # rho = 0.4629 S = 3834.3, p-value = 0.005099
# plot(bc.res~data$type, col = c("gold", "cyan"))

data.t <- data[!is.na(data$raw_T_pg.hr),]  # need to get rid of the two NAs

# plot(data.t$raw_T_pg.hr ~ data.t$length, col = ifelse(data.t$type=="S","cyan", "gold"), pch = 19)
# plot(data.t$raw_T_pg.hr ~ data.t$type, col = c("gold", "cyan"))

# t.test(data.t$raw_T_pg.hr~data.t$type ) # p-value = 0.0002582  but need to correct for size

t_byLeng <- data.t$raw_T_pg.hr/(data.t$length/10)  ## These are the T values to be used in correlations
data.t <- cbind(data.t, t_byLeng)

# plot(t_byLeng ~ data.t$length, col = ifelse(data.t$type=="S","cyan", "gold"), pch = 19)
# plot(t_byLeng ~ data.t$type, col = c("gold", "cyan"))  ## make this one in ggplot

#Plot corrected testosterone by type to confirm experimental paradigm
pte = ggplot(data.t, aes(x = type, y = t_byLeng, group = type)) 
pte + geom_boxplot(data = data.t, aes(x = type, fill = type)) + 
  labs(x = "Current Phenotype", y = "Testosterone pg/hr/mm")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))+
  geom_point()


#Compare behavior scores by type to confirm experimental paradigm
summary(aov(data.t$t_byLeng  ~data.t$type )) # p-value = 0.0424 F = 4.482 df = 1,31 
data <- read.table(file = "AS_metadata_filt_newT_20181205.csv", header = TRUE, sep = "," )
plot(behavior_scores~type, data = data)  #### Make this one in ggplot
p = ggplot(data, aes(x = type, y = behavior_scores, group = type)) 
p + geom_boxplot(data = data, aes(x = type, fill = type)) + 
  labs(x = "Current Phenotype", y = "Behavior Score")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))+
  geom_point()

#Paramertic and nonparametric tests
summary(aov(behavior_scores~type, data = data))
wilcox.test(behavior_scores~type, data = data, paired = FALSE )

#plot(weight~length, data = data)  ## roughly linear, log linear?

#plot(log(weight)~length, col = ifelse(prev_type=="S","cyan", "gold"), pch = 19, data = data) ## looks linear

bc.lm <- lm (log(weight)~length, data = data)
bc.res <- resid(bc.lm)
as.matrix(bc.res)  # These are the BC values that Josh should use for correlations
#### they are added to the file "AS_metadata_filt_newT_20181205.csv"

# cor.test (bc.res, data$behavior_scores, method = "pearson") # R2 = 0.311125, t = 1.8806, df = 33, p-value = 0.06887
cor.test (bc.res, data$behavior_scores, method = "spearman") # rho = 0.4629 S = 3834.3, p-value = 0.005099

# plot(bc.res~data$type, col = c("gold", "cyan"))

data.t <- data[!is.na(data$raw_T_pg.hr),]  # need to get rid of the two NAs

# plot(data.t$raw_T_pg.hr ~ data.t$length, col = ifelse(data.t$type=="S","cyan", "gold"), pch = 19)
# plot(data.t$raw_T_pg.hr ~ data.t$type, col = c("gold", "cyan"))
# t.test(data.t$raw_T_pg.hr~data.t$type ) # p-value = 0.0002582  but need to correct for size

t_byLeng <- data.t$raw_T_pg.hr/data.t$length  ## These are the T values to be used in correlations
data.t <- cbind(data.t, t_byLeng)
# plot(t_byLeng ~ data.t$length, col = ifelse(data.t$type=="S","cyan", "gold"), pch = 19)
# plot(t_byLeng ~ data.t$type, col = c("gold", "cyan"))  ## make this one in ggplot
p = ggplot(data.t, aes(x = type, y = t_byLeng, group = type)) 
p + geom_boxplot(data = data.t, aes(x = type, fill = type)) + 
  labs(x = "Current Phenotype", y = "Testosterone pg/hr/mm")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))+
  geom_point()
summary(aov(data.t$t_byLeng  ~data.t$type )) # p-value = 0.0424 F = 4.482 df = 1,31 


# Plot all social hierarchy proxy measures for Current Type
p = ggplot(sample_data(phyloseq_data), aes(x = type, y = behavior_scores, group = type)) 
pbs = p + geom_boxplot(data = sample_data(phyloseq_data), aes(x = type, fill = type)) + 
  labs(x = "Current Phenotype", y = "Behavior Score")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))+
  geom_point() + 
  theme(legend.position="none")

p = ggplot(sample_data(phyloseq_data), aes(x = type, y = t_byLeng, group = type)) 
pte = p + geom_boxplot(data = sample_data(phyloseq_data), aes(x = type, fill = type)) + 
  labs(x = "Current Phenotype", y = "Testosterone pg/hr/mm")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))+
  geom_point() + 
  theme(legend.position="none")

p = ggplot(sample_data(phyloseq_data), aes(x = type, y = bc.res, group = type)) 
pbc = p + geom_boxplot(data = sample_data(phyloseq_data), aes(x = type, fill = type)) + 
  labs(x = "Current Phenotype", y = "Residual Body Condition")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))+
  geom_point() + 
  theme(legend.position="none")

grid.arrange(pbs, pte, pbc, nrow = 1)

################################################### PhyloFactor ##############################################################

##########################################################################
# In addition to the DESeq2 analysis to test for differential abundance, #
# this analysis can detect either individual OTUs or entire clades that  #
# explain differences among groups                                       #
##########################################################################

#Import abundance table
#OTUs in table must be sorted to reflect the order of OTUs in the FastTree phylogeny
#This table is also sample-sorted by treatment for ease of making figures
otufile <- read.table("AS_sample_PF_ordered_OTU_table.txt", header = TRUE, row.names=1)
#OPTIONAL: remove columns with NA values for prev_type if that's the metadata variable being used for analysis
otufile <- select(otufile, -NA.AS12.1)
otufile <- select(otufile, -D.AS10.2)
otu <- data.matrix(otufile)

# rarefy data!  ### recommended normalization for when different samples
# have different read depths. Keep checking Phylofactor documentation to see
# if this changes, since rarification is not always ideal, but often good enough
otur <- rrarefy(t(otu), 3879)  # from vegan package (set number for dataset) (transpose)
otu <- t(otur)  # transpose back just to use commands

#Import unrooted 16S tree for OTUs created by FastTree
pf_tree<-read.tree("AS_phylogeny.tre")  # must be an unrooted tree for phylofactor

#Import metadata
mapfile <- read.table("AS_metadata_filt_PF_ordered.txt", sep = "\t", header = TRUE, row.names=1)
#OPTIONAL: remove row with NA value for prev_type if that's the metadata variable being used for analysis
mapfile <- subset(mapfile, prev_type!="NA")
mapfile <- subset(mapfile, outlier!="bad")

# phyloseq needs rownames
meta_table_pf<-sample_data(mapfile)
head(mapfile)
head(meta_table_pf) ## they look the same
#pull out metadata variables to perform phylofactorization
#we used meta_list1 for final analysis
meta_list1 <- meta_table_pf$prev_type
meta_list2 <- meta_table_pf$group
meta_list3 <- meta_table_pf$timepoint
meta_list4 <- meta_table_pf$type

#Import taxa designations for OTUs
pf_tax<-data.frame(read.table("AS_pf_tax2.txt",header=TRUE, sep = "\t"))

###OPTIONAL, filter out OTUs that are in more than 5 samples###
# Choose your own number based on portion of samples you're analyzizng
# This set of commands also prunes the phylogeny, you can prune after
# running phylofactor as well to make heatmaps more readable with common
# taxa only
common.otus <- which(rowSums(otu>0)>5)
otu <- otu[common.otus,]
pf_tree <- ape::drop.tip(pf_tree,setdiff(pf_tree$tip.label,rownames(otu)))
pf_tax <- pf_tax[match(pf_tree$tip.label,pf_tax[,1]),]
dim(otu)

#Perform Phylofactorization
#you can change number of factors to look at additional phylogenetic dimensions within dataset
#you can also change meta_list based on which variables you want to consider in PF analysis

#first, set order of otu to order of tip labels in tree for easy downstream processing
pf_tree <- ape::drop.tip(pf_tree,setdiff(pf_tree$tip.label,rownames(otu)))
otu <- otu[pf_tree$tip.label,]
#for first pass, set nfactors to a large number, ex. >30
#this will allow you to plot a scree-plot of explained variation
#from each factor, then look for elbow
#pick meta_list that best defines the comparison you want to make
PF <- PhyloFactor(otu,pf_tree,meta_list1,nfactors=30)

#### Look at summary and determine the number of factors you want to look at for your second pass based on ExpVar and p-values #####
#### For this study, we looked for the point at which we saw 5 consecutive non-significant factors, and stopped factorization there #####

#for group
#PF <- PhyloFactor(otu,pf_tree,meta_list2,nfactors=30)
#Using Variance Stabilized transformations
#PF <- PhyloFactor(diagvst,pf_tree,meta_list,nfactors=30)

#Look at summary of PF results, specifically which phylogenetic clades and OTUs explain variable of interest 
#This table has many of the main results from PF like ExpVar, F, and Pr(>F)
PF$factors

#####################################################################################
# creates full data summary for factor of interest, and then allows you to look at  #
# relative abundance ratios of the OTU groups defined by each factor on a sample to #
# sample basis                                                                      #
#####################################################################################

# pick factor to examine more in-depth, it will be labeled as group 1 in results
smry <- pf.summary(PF,pf_tax,factor=1)
pf.tidy(smry)

#Allows you to visualize the abundance ratios between OTU groups for the summarized 
#factor chosen above by "pf.summary" function in each of your sample variable categories

# names(smry)
par(mfrow=c(1,2))
plot(meta_list1,smry$ilr,main='Isometric log-ratio',ylab='ILR balance')
plot(meta_list1,smry$MeanRatio,main='Ratio of Geometric Means',ylab='Previous Phenotype')

#Look at the OTU group from the summarized factors 
#find taxonomy
smry$TaxaSplit %>%
  lapply(.,FUN=function(x) unique(x$TaxaIDs)) 
#find OTU names
Group1.otus <- smry$group1$IDs$otuIDs %>% as.list %>% sapply(.,toString)
Group1.otus

##########################################################
# Create heatmap of OTU and clade abundance vs predicted #
##########################################################

# To look at phylogeny only:
#par(mfrow=c(1,1))
#plot.phylo(pf_tree,use.edge.length=FALSE,main='Community Phylogeny')

#create single heatmap of log relative abundance of OTUs based on phylogenetic relationships.
#it takes a while to run with all OTUs, so pruning tree is recommended for time and readability 
colors<-colorRampPalette(colors=c("white", "blue"))(100)
clr <- function(Matrix) apply(Matrix,MARGIN=2,FUN=function(x) log(x)-mean(log(x)))
par(mfrow=c(1,1))
# heatmap for IRL values (log abundance relative to other OTUs within sample)
phylo.heatmap(pf_tree,colors=colors,clr(PF$Data))
# heatmap for raw abundance, but tends to be skewed from high adundance taxa, not very readable
#phylo.heatmap(pf_tree,colors=colors,clr(PF$Data))


#######################################
#         Plots in Manuscript         #
#######################################

# Create multiple heatmaps, showing raw log relative abundance and averages for factors and groups
# Predictions are the average of factors within groups
####### to make smaller heatmap with unfiltered OTU dataset, run PF for full dataset then prune phylogeny for input into phylo.heatmap #########
clr <- function(Matrix) apply(Matrix,MARGIN=2,FUN=function(x) log(x)-mean(log(x)))
colors<-colorRampPalette(colors=c("white", "blue","blue4", "darkblue", "navy"))(100)
prediction <- pf.predict(PF)
dim(prediction)
prediction <- pf.predict(PF,factors=21)
par(mfrow=c(2,1))
phylo.heatmap(pf_tree,clr(PF$Data),colors=colors, fsize=0.5)
phylo.heatmap(pf_tree,clr(prediction),colors=colors, fsize=0.5)

################################################ Other useful plots ############################################################

#######################################
#         Plots in Manuscript         #
#######################################

#for manuscript: make phylum abundance plot by previous type
#be sure to load the phylofactor data before creating this phyloseq object; 
#Change taxonomic level or variable as desired for exploration of data

plot_subset <- subset_samples(phyloseq_data_norarify, prev_type!="NA")
plot_bar(plot_subset, fill="Phylum")

#group taxa by phylum for easier visualization, less complexity
phylumGlommed = tax_glom(phyloseq_data_norarify, "Phylum")
plot_bar(phylumGlommed, fill ="Phylum")

otu <- otu_table(t(otur), taxa_are_rows = TRUE) 
pfdata <- phyloseq(otu, meta_table_pf, ep_tax, ep_tree_rooted)

phylumGlommed = tax_glom(pfdata, "Phylum")
plot_bar(phylumGlommed, fill ="Phylum")

#convert unrarified data into percents for easier plotting
otu_percent <- otu_table(colPercents(otu_table(phyloseq_data_norarify),digits = 3), taxa_are_rows = TRUE)
pfdata <- phyloseq(otu_percent, meta_table, ep_tax, ep_tree_rooted)
plot_subset <- subset_samples(pfdata, prev_type!="NA")
plot_bar(plot_subset, fill="Phylum")

#group taxa by phylum for easier visualization, less complexity
phylumGlommed = tax_glom(pfdata, "Phylum")
plot_bar(phylumGlommed, fill ="Phylum")

#plot DESeq2 transformed data to percent for plotting
vst_percent <- otu_table(colPercents(diagvst,digits = 3), taxa_are_rows = TRUE)
meta_table_noNA <- subset(meta_table, prev_type!="NA")
pfdata <- phyloseq(vst_percent, meta_table_noNA, ep_tax, ep_tree_rooted)
plot_subset <- subset_samples(pfdata, prev_type!="NA")
plot_bar(plot_subset, fill="Phylum")

#group taxa by phylum for easier visualization, less complexity
phylumGlommed = tax_glom(pfdata, "Phylum")
plot_bar(phylumGlommed, fill ="Phylum")

