#############################################################################
#       Microbiome analysis for dominant and subordinate A. burtoni         #
#     Last updated 10/12/2018 by Josh Faber-Hammond and Avehi Singh         #
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

theme_set(theme_bw())
set.seed(42)

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



############################### Begin by importing input files and normalizing data for comparison of samples ################################

#############################
## Create Phyloseq object  ##
#############################

setwd("/Volumes/Bay_2/avehi/microbiome_final")
otufile <- read.table("AS_sample_filtered_OTU_table.txt", sep = "\t", header = TRUE, row.names=1)
## just import and then tell phyloseq that these are the files that it is looking for

mapfile <- read.table("AS_metadata_filt.txt", sep = "\t", header = TRUE, row.names=1)
# phyloseq needs rownames
meta_table<-sample_data(mapfile)
head(mapfile)
head(meta_table) ## they look the same

ep_tax<-read.table("AS_tax.txt", header=TRUE, sep = "\t", row.names=1) 
# phyloseq needs rownames
ep_tax <- as.matrix(ep_tax)  # coerce to matrix
ep_tax <- tax_table(ep_tax)

ep_tree<-read.tree("AS_phylogeny.tre")  # must be a rooted tree for phyloseq
is.rooted(phy_tree(ep_tree)) 
ep_tree_rooted <-midpoint(ep_tree) # this is from phanghorn package that sets root to midpoint
is.rooted(phy_tree(ep_tree_rooted)) 

#rarefy data! This is optional and related to normalization when different samples have different read depths 
#these argumenents mattered more for 454
# subsample to minimum count for samples passing filters.
# make this number match your read count/filtered datasets
otur <- rrarefy(t(otufile),3879)  # from vegan package (set number for dataset) (transpose)

#write.table(otur, file = "rarified_otu_table.txt", sep = "\t")

otu <- otu_table(t(otur), taxa_are_rows = TRUE)  # transpose back just to use commands

#lastly, turn all your inputs into phyloseq object
phyloseq_data <- phyloseq(otu, meta_table, ep_tax, ep_tree_rooted)

###########################
# For unrarified raw data #
###########################

otu <- otu_table(otufile, taxa_are_rows = TRUE)  # transpose back just to use commands
phyloseq_data_norarify <- phyloseq(otu, meta_table, ep_tax, ep_tree_rooted)

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
res <- results( diagdds, cooksCutoff = TRUE )
res <- results(diagdds, contrast = c("prev_type", "D", "S"), cooksCutoff = TRUE )

# Filter for your desired alpha value, and extract relevant taxon names
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
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
#phyloseq_data_ds2 <- phyloseq(phyloseq_data_ds2_count, meta_table, ep_tax, ep_tree_rooted)

#optional output from DESeq2 normalized data
#write.table(diagvst, file = "DESeq2_table_correlational_analysis_floored.txt", sep = "\t")

####################################################################
# Subset the dataset for downstream analyses by metadata variables #
####################################################################
#subsetting used for exploratory analysis
subset_dd <- subset_samples(phyloseq_data,group=="DD")
subset_ss <- subset_samples(phyloseq_data,group=="SS")
subset_ds <- subset_samples(phyloseq_data,group=="DS")
subset_sd <- subset_samples(phyloseq_data,group=="SD")
dd_no <- subset_samples(subset_dd,outlier=="ok")
ds_no <- subset_samples(subset_ds,outlier=="ok")
sd_no <- subset_samples(subset_sd,outlier=="ok")
ss_no <- subset_samples(subset_ss,outlier=="ok")
dd_paired <- subset_samples(subset_dd,pair=="y")

#this is unrarified data without the cyanobacteria sample
nooutlier <- subset_samples(phyloseq_data_norarify,outlier=="ok")
noout_nocyan <- subset_samples(nooutlier,pair=="y")

#here we are removing unknown previous type and removing outliers with cyanobacteria
noprevtype <- subset_samples(phyloseq_data,prev_type!="NA")
noprev_nocyan <- subset_samples(noprevtype,outlier=="ok")

#use this subset for correlational analysis
correlation_subset <- subset_samples(phyloseq_data_ds2, filter=="y")

#define subset to test
microbiome_subset <- noprev_nocyan #this is the subset we used for beta diversity analysis
microbiome_subset <- phyloseq_data
microbiome_subset <- phyloseq_data_norarify
microbiome_subset <- phyloseq_data_ds2 #do not run with bray or jaccard distance methods unless counts floored at 0
microbiome_subset <- noout_nocyan
microbiome_subset <- correlation_subset
microbiome_subset <- subset_dd
microbiome_subset <- subset_ss
microbiome_subset <- subset_ds
microbiome_subset <- subset_sd
microbiome_subset <- dd_no
microbiome_subset <- ds_no
microbiome_subset <- sd_no
microbiome_subset <- ss_no

################################################## BETA DIVERSITY #######################################################

#ordinations! try different methods and distance metrics
#do not use OTU tables with negative values
ordu <- ordinate(microbiome_subset,"NMDS", "bray")
ordu2 <- ordinate(microbiome_subset,"NMDS", "unifrac")
ordu3 <- ordinate(microbiome_subset,"NMDS", "wunifrac")
ordu4 <- ordinate(microbiome_subset,"NMDS", "jaccard")
ordu5 <- ordinate(microbiome_subset,"PCoA", "bray")
ordu6 <- ordinate(microbiome_subset,"PCoA", "unifrac")
ordu7 <- ordinate(microbiome_subset,"PCoA", "wunifrac")
ordu8 <- ordinate(microbiome_subset,"PCoA", "jaccard")
ordu9 <- ordinate(microbiome_subset,"DPCoA")

plot_ordination(microbiome_subset, ordu2, type = "scree")
stressplot(ordu)


#######################################
#         Plots in Manuscript         #
#######################################

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

############################################
#beta diversity statistsical tests
#I tested unifrac and wunifrac, so change distance method for each
metaD = as(sample_data(microbiome_subset), "data.frame")
phydist = distance(microbiome_subset, "wunifrac")

#ANOSIM, fit data to single variable groupings/distributions
#This tests whether two or more groups of samples are significantly different
anosim_res = anosim(phydist, metaD$bc_linear, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$group, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$type, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$timepoint, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$group_timepoint, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$behavior_scores, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$prev_type, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$testosterone_nc, distance = "wunifrac")
anosim_res

#PERMANOVA by different variable groups
#This tests for differences in centroids in multidimensional matrix
adonis_res <- adonis(formula = phydist ~ metaD$type, data = metaD) 
adonis_res
adonis_res <- adonis(formula = phydist ~ metaD$group, data = metaD) 
adonis_res
adonis_res <- adonis(formula = phydist ~ metaD$bc_linear, data = metaD) 
adonis_res
adonis_res <- adonis(formula = phydist ~ metaD$timepoint, data = metaD) 
adonis_res
adonis_res <- adonis(formula = phydist ~ metaD$group:metaD$timepoint, data = metaD) 
adonis_res
adonis_res <- adonis(formula = phydist ~ metaD$group + metaD$timepoint + metaD$group:metaD$timepoint, data = metaD) 
adonis_res
adonis_res <- adonis(formula = phydist ~ metaD$behavior_scores, data = metaD) 
adonis_res
adonis_res <- adonis(formula = phydist ~ metaD$prev_type, data = metaD) 
adonis_res
adonis_res <- adonis(formula = phydist ~ metaD$testosterone_nc, data = metaD) 
adonis_res
adonis_res <- adonis(formula = phydist ~ metaD$testosterone_nc + metaD$weight, data = metaD) 
adonis_res
adonis_res <- adonis(formula = phydist ~ metaD$tank:metaD$timepoint, data = metaD) 
adonis_res

#This tests for differences in variances of clusters based on variables
beta <- betadisper(phydist, metaD$type)
permutest(beta)

beta <- betadisper(phydist, metaD$group)
permutest(beta)

beta <- betadisper(phydist, metaD$timepoint)
permutest(beta)

beta <- betadisper(phydist, metaD$prev_type)
permutest(beta)

################################################### alpha diversity ############################################################

#estimate species richness/diversity
#subset data so that there are no NA values in the metadata column for correlations
alpha_subset <- subset_samples(phyloseq_data_norarify,filter == "y")
mapfile2 <- read.table("AS_metadata_noNA_20180814.txt", sep = "\t", header = TRUE, row.names=1)
meta_table2<-sample_data(mapfile2)
meta_table3 <- data.frame(meta_table2)

#alpha diversity
#test for cvorrelations between differentially abundant OTUs, alpha diversity, and metadata variables
richness <- estimate_richness(alpha_subset, split = TRUE, measures = NULL)
richness <- cbind(richness, prev_type = meta_table2$prev_type, behavior_scores = meta_table2$behavior_scores, weight = meta_table2$weight, bc_linear = meta_table2$bc_linear,  testosterone_nc = meta_table2$testosterone_nc, otu4_deseq = meta_table2$deseq_norm_otu4, otu6_deseq = meta_table2$deseq_norm_otu6, otu125_deseq = meta_table2$OTU_125, otu438_deseq = meta_table2$OTU_438, otu376_deseq = meta_table2$OTU_376, otu474_deseq = meta_table2$OTU_474, otu155_deseq = meta_table2$OTU_155, otu190_deseq = meta_table2$OTU_190, otu322_deseq = meta_table2$OTU_322, otu80_deseq = meta_table2$OTU_80, otu125_nf_deseq = meta_table2$OTU_125_nf, otu438_nf_deseq = meta_table2$OTU_438_nf, otu376_nf_deseq = meta_table2$OTU_376_nf, otu474_nf_deseq = meta_table2$OTU_474_nf, otu155_nf_deseq = meta_table2$OTU_155_nf, otu190_nf_deseq = meta_table2$OTU_190_nf, otu322_nf_deseq = meta_table2$OTU_322_nf, otu80_nf_deseq = meta_table2$OTU_80_nf)
richness
# write table with diversity stats and metadata
write.table(richness, file = "richness+meta.txt", sep = "\t")

#to actually test correlations between all fields in the diversity table, make sure your dataset has no "NA" values
#filter samples by subset function prior to correlation check
alphacorr <- rcorr(as.matrix(richness))
alphacorr

# write tables with correlation results
write.table(alphacorr$P, file = "alpha_div_correlations_wOTUs_Pval.txt", sep = "\t")
write.table(alphacorr$r, file = "alpha_div_correlations_wOTUs_rval.txt", sep = "\t")

#make plots of correlations
richness_df <- data.frame(richness) #just need to run once before making plots

#plot alpha diversity for metadata variables
p = plot_richness(alpha_subset, x = "prev_type",
                  title = NULL, scales = "free_y", nrow = 1, shsi = NULL,
                  measures = c("Observed", "Shannon"), sortby = NULL)
p + geom_boxplot(data = p$alpha_subset, aes(x = prev_type, y = value, color = NULL), alpha = 0.1) + labs(x = "Previous Phenotype")



#######################################
#         Plots in Manuscript         #
#######################################

#make boxplot with colours!
p = plot_richness(alpha_subset, x = "prev_type",
                  title = NULL, scales = "free_y", nrow = 1, shsi = NULL,
                  measures = c("Observed", "Shannon"), sortby = NULL)
p + geom_boxplot(data = p$alpha_subset, aes(x = prev_type, y = value, fill = prev_type)) + 
  labs(x = "Previous Phenotype")+ 
  scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))

#######################################
#         Plots in Manuscript         #
#######################################

#subset table to separate subordinate and dominant and test for OTU abundances and test for correlations
subordinate <- subset(richness, prev_type == 'S')
dominant <- subset(richness, prev_type == 'D')

#make correlations for subordinate and dominant separately
subordinate <- subset(subordinate, select = -c(prev_type) )
sub_alphacorr <- rcorr(as.matrix(subordinate))
sub_alphacorr
write.table(sub_alphacorr$P, file = "alpha_div_correlations_subordinates_Pval.txt", sep = "\t")
write.table(sub_alphacorr$r, file = "alpha_div_correlations_subordinates_rval.txt", sep = "\t")

dominant <- subset(dominant, select = -c(prev_type) )
dom_alphacorr <- rcorr(as.matrix(dominant))
dom_alphacorr
write.table(dom_alphacorr$P, file = "alpha_div_correlations_dominants_Pval.txt", sep = "\t")
write.table(dom_alphacorr$r, file = "alpha_div_correlations_dominants_rval.txt", sep = "\t")

#par(mfrow=c(1,2))
richness_manuscript <- cbind(richness, prev_type_col = meta_table2$prev_type)
richness_ms <- as.matrix(richness_manuscript)
richness_ms <- richness_manuscript

# colour by prev_type 
richness_ms$prev_type_col <- as.character(richness_ms$prev_type_col)
richness_ms$prev_type_col <- replace(richness_ms$prev_type_col, richness_ms$prev_type_col =="D", "darkgoldenrod2")
richness_ms$prev_type_col <- replace(richness_ms$prev_type_col, richness_ms$prev_type_col =="S", "turquoise4")

# Plot all points at once, using newly generated colours
#scatter plots of alpha diversity vs. differentially abundant OTU counts
plot(richness_ms$Observed,richness_ms$otu4_deseq, xlab = "Richness", ylab = "OTU4 normalized abundance", col=richness_ms$prev_type_col, pch = 19, cex = 1.5)
abline(lm(richness_ms$otu4_deseq~richness_ms$Observed), col="black") # regression line (y~x) 
abline(lm(subordinate$otu4_deseq~subordinate$Observed), col="turquoise4") # regression line (y~x) 
abline(lm(dominant$otu4_deseq~dominant$Observed), col="darkgoldenrod2") # regression line (y~x) 
#lines(lowess(richness_ms$otu4_deseq~richness_ms$Observed), col="blue") # lowess line (x,y)

plot(richness_ms$Observed,richness_ms$otu6_deseq, xlab = "Richness", ylab = "OTU6 normalized abundance", col=richness_ms$prev_type_col, pch = 19, cex = 1.5)
abline(lm(richness_ms$otu6_deseq~richness_ms$Observed), col="black") # regression line (y~x) 
abline(lm(subordinate$otu6_deseq~subordinate$Observed), col="turquoise4") # regression line (y~x) 
abline(lm(dominant$otu6_deseq~dominant$Observed), col="darkgoldenrod2") # regression line (y~x) 
#lines(lowess(richness_ms$otu6_deseq~richness_ms$Observed), col="blue") # lowess line (x,y)

###############################################
# Better way to estimate alpha diverstiy:     #
# Perform multiple iterations of rarification #
# before calulating diversity stats           #
###############################################

# Initialize matrices to store richness and evenness estimates
# Make sure subset is your unrarified dataset
microbiome_nonrarify_subset <- phyloseq_data_norarify
min_lib <- min(sample_sums(microbiome_nonrarify_subset))
nsamp = nsamples(microbiome_nonrarify_subset)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(microbiome_nonrarify_subset)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(microbiome_nonrarify_subset)

shannon <- matrix(nrow = nsamp, ncol = trials)
row.names(shannon) <- sample_names(microbiome_nonrarify_subset)

simpson <- matrix(nrow = nsamp, ncol = trials)
row.names(simpson) <- sample_names(microbiome_nonrarify_subset)

chao1 <- matrix(nrow = nsamp, ncol = trials)
row.names(chao1) <- sample_names(microbiome_nonrarify_subset)
chao1 <- rbind(chao1, chao1)

# It is always important to set a seed when you subsample so your result is replicable 
set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(microbiome_nonrarify_subset, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
  
  # Calculate shannon
  shan <- as.numeric(as.matrix(estimate_richness(r, measures = "Shannon")))
  shannon[ ,i] <- shan
  
  # Calculate simpson
  simp <- as.numeric(as.matrix(estimate_richness(r, measures = "Simpson")))
  simpson[ ,i] <- simp
  
  # Calculate Chao1
  chao <- as.numeric(as.matrix(estimate_richness(r, measures = "Chao1")))
  chao1[ ,i] <- chao
  
}

# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean_alpha <- apply(richness, 1, mean)
se <- apply(richness, 1, ste)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean_alpha, se, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness)
mean_alpha <- apply(evenness, 1, mean)
se <- apply(evenness, 1, ste)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean_alpha, se, measure)

# Create a new dataframe to hold the means and standard deviations of diversity estimates
SampleID <- row.names(simpson)
mean_alpha <- apply(simpson, 1, mean)
se <- apply(simpson, 1, ste)
measure <- rep("Simpson", nsamp)
simp_stats <- data.frame(SampleID, mean_alpha, se, measure)

# Create a new dataframe to hold the means and standard deviations of diveristy estimates
SampleID <- row.names(shannon)
mean_alpha <- apply(shannon, 1, mean)
se <- apply(shannon, 1, ste)
measure <- rep("Shannon", nsamp)
shan_stats <- data.frame(SampleID, mean_alpha, se, measure)

# Create a new dataframe to hold the means and standard deviations of diversity estimates
SampleID <- row.names(chao1)
mean_alpha <- apply(chao1, 1, mean)
se <- apply(chao1, 1, ste)
measure <- rep("Chao1", nsamp)
chao_stats <- data.frame(SampleID, mean_alpha, se, measure)

#combine your stats and export tables
alpha <- rbind(rich_stats, even_stats, simp_stats, shan_stats)

s <- data.frame(sample_data(microbiome_subset))
alphadiv <- merge(alpha, s, by = "row.names") 

write.csv(alphadiv, "alphadiv.csv")

#chao_stats contains both means and variances based on the initial calculation of chao1, 
#therefore we need to export separate since it is not a 1:1 join with n samples.
#be cautious with how you treat means and SD of the original chao1 SD calculations.  They may be meaningless depending on your goals.
#chao1 mean stats are written first, followed by SD calculations
write.csv(chao_stats, "chao1.csv")

#################################################
# Run ANOVA to test whether alpha diversity and #
# BC/GSI stats differ among treatment groups    #
#################################################

rich2 <- data.frame(merge(rich_stats, s, by = "row.names"))
even2 <- data.frame(merge(even_stats, s, by = "row.names"))
simp2 <- data.frame(merge(simp_stats, s, by = "row.names"))
shan2 <- data.frame(merge(shan_stats, s, by = "row.names"))

fit <- aov(mean_alpha ~ prev_type, data=shan2)
summary(fit)

fit <- aov(mean_alpha ~ prev_type, data=rich2)
summary(fit)

fit <- aov(mean_alpha ~ prev_type, data=even2)
summary(fit)

fit <- aov(mean_alpha ~ prev_type, data=simp2)
summary(fit)

# Run pairwise tests to see specifically which groups differ
fit <- pairwise.t.test(rich2$mean_alpha, rich2$group:rich2$timepoint, p.adj = "bonf")
fit

fit <- pairwise.t.test(even2$mean_alpha, even2$group:even2$timepoint, p.adj = "bonf")
fit

fit <- pairwise.t.test(simp2$mean_alpha, simp2$group:simp2$timepoint, p.adj = "bonf")
fit

fit <- pairwise.t.test(shan2$mean_alpha, shan2$group:shan2$timepoint, p.adj = "bonf")
fit

#test assumption of normality for ANOVAs
shapiro.test(shan2$mean_alpha)
shapiro.test(simp2$mean_alpha)
shapiro.test(even2$mean_alpha)
shapiro.test(rich2$mean_alpha)



#######################################
#         Plots in Manuscript         #
#######################################

#make bar plots for testosterone vs. current phenotype with colours
AVG<-aggregate(testosterone_nc ~ type  , meta_table3, mean, na.rm = T) 
SE<-aggregate(testosterone_nc ~ type  , meta_table3, ste) 
ptypetest <- data.frame(merge(AVG, SE, by = "type"))

p <- ggplot(ptypetest, aes(x=factor(type), y=testosterone_nc.x, fill = type)) + 
  stat_summary(fun.y="mean", geom="bar") + 
  geom_errorbar(aes(ymin=ptypetest$testosterone_nc.x - ptypetest$testosterone_nc.y, ymax=ptypetest$testosterone_nc.x + ptypetest$testosterone_nc.y),
                width=0.2,                    # Width of the error bars
                position=position_dodge(.9)) + labs(x = "Current Phenotype", y = "Testosterone (pg/hr)") 
p + scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))

#make bar plots for behaviour vs. current phenotype
AVG<-aggregate(behavior_scores ~ type  , meta_table3, mean, na.rm = T) 
SE<-aggregate(behavior_scores ~ type  , meta_table3, ste) 
bsvstype <- data.frame(merge(AVG, SE, by = "type"))

x <- ggplot(bsvstype, aes(x=factor(type), y=behavior_scores.x, fill = type)) + 
  stat_summary(fun.y="mean", geom="bar") + 
  geom_errorbar(aes(ymin=bsvstype$behavior_scores.x - bsvstype$behavior_scores.y, ymax=bsvstype$behavior_scores.x + bsvstype$behavior_scores.y),
                width=0.2,                    # Width of the error bars
                position=position_dodge(.9)) + labs(x = "Current Phenotype", y = "Aggression Score")
x + scale_fill_manual(values = c('darkgoldenrod2', 'turquoise4'))

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
#OPTIONAL: remove columns with NA values for prev_type if that's the matadata variable being used for analysis
#otufile <- select(otufile, -NA.AS12.1)
otu <- data.matrix(otufile)

# rarefy data!  ### recommended normalization for when different samples
# have different read depths. Keep checking Phylofactor documentation to see
# if this changes, since rarification is not always ideal, but often good enough
otur <- rrarefy(t(otufile), 3879)  # from vegan package (set number for dataset) (transpose)
otu <- otu_table(t(otur), taxa_are_rows = TRUE)  # transpose back just to use commands

#Import unrooted 16S tree for OTUs created by FastTree
pf_tree<-read.tree("AS_phylogeny.tre")  # must be an unrooted tree for phylofactor

#Import metadata
mapfile <- read.table("AS_metadata_filt_PF_ordered.txt", sep = "\t", header = TRUE, row.names=1)
#OPTIONAL: remove row with NA value for prev_type if that's the matadata variable being used for analysis
#mapfile <- subset(mapfile, prev_type!="NA")

# phyloseq needs rownames
meta_table<-sample_data(mapfile)
head(mapfile)
head(meta_table) ## they look the same
#pull out metadata variables to perform phylofactorization
#we used meta_list1 for final analysis
meta_list1 <- meta_table$prev_type
meta_list2 <- meta_table$group
meta_list3 <- meta_table$timepoint
meta_list4 <- meta_table$type

#Import taxa designations for OTUs
pf_tax<-data.frame(read.table("AS_pf_tax.txt",header=TRUE, sep = "\t"))

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
otu <- otu[pf_tree$tip.label,]
#for first pass, set nfactors to a large number, ex. >20
#this will allow you to plot a scree-plot of explained variation
#from each factor, then look for elbow
#pick meta_list that best defines the comparison you want to make
PF <- PhyloFactor(otu,pf_tree,meta_list1,nfactors=6)
#for group
#PF <- PhyloFactor(otu,pf_tree,meta_list2,nfactors=20)
#Using Variance Stabilized transformations
#PF <- PhyloFactor(diagvst,pf_tree,meta_list,nfactors=20)

#Look at summary of PF results, specifically which phylogenetic clades and OTUs explain variable of interest 
#This table has many of the main results from PF like ExpVar, F, and Pr(>F)
PF$factors

#######################################
# This next set of commands will help #
# you know when to stop factorization #
#######################################

# Find the elbow of the scree plot for ExpVar
# To do so, you will preform loess smoothing of the ExpVar plot
# then look for  the peak in the change of slope

# For this exploratory analysis you need to use a larger number 
# of factors than you expect to be biologically significant to 
# identify a somewhat empirical cutoff (ex: n=20)

elb <- PF$factors
elb$Factors <- 1:nrow(elb)
y <- elb$ExpVar
x <- elb$Factors
par(mfrow=c(1,2))
plot(x,y,type="l")
lo <- loess(y~x)
#xl <- seq(min(x),max(x), (max(x) - min(x))/19) # set the denominator to (nfactors - 1)
out = predict(lo,x)
lines(x, out, col='red', lwd=2)

findElbow(out, plot = TRUE, returnIndex = TRUE)
par(mfrow=c(1,1))

#This function finds the elbow of the curve. Take this value and rerun phylofactor 
#with this as the number of factors to report. This method is not perfect but the 
#best way at getting to the number of factors that may be biologically interesting. 
#The author of phylofactor recognizes that one of the main limitations of his program 
#is that there's no systematic way to know when to stop factorization, and this is 
#our best attempt at solving the problem. One remaining issue with this method is that
#it requires a non-noisy concave curve and our loess smoothing may over-smooth the 
#data to the point where it doesn't quite match the original biological data, so 
#proceed with caution. Still better than nothing though.

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
# to make smaller heatmap with unfiltered OTU dataset, run PF for full dataset then prune phylogeny for input into phylo.heatmap
clr <- function(Matrix) apply(Matrix,MARGIN=2,FUN=function(x) log(x)-mean(log(x)))
colors<-colorRampPalette(colors=c("white", "blue","blue4", "darkblue", "navy"))(100)
prediction <- pf.predict(PF)
dim(prediction)
prediction <- pf.predict(PF,factors=12)
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
pfdata <- phyloseq(otu, meta_table, ep_tax, ep_tree_rooted)
plot_subset <- subset_samples(pfdata, prev_type!="NA")
plot_bar(plot_subset, fill="Phylum")
