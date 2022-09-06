# ------------------------------------------------------------------------------------------------------------#
# Manuscript title: Mapping tropical forest functional variation at satellite remote sensing resolutions 
# depends on key traits
# ------------------------------------------------------------------------------------------------------------#
# Manuscript authors: Elsa M. Ordway, Gregory P. Asner, David F.R.P. Burslem, Simon L. Lewis, Reuben Nilus,
# Roberta Martin, Michael J. O’Brien, Oliver L. Phillips, Lan Qie, Nicholas R. Vaughn, Paul R. Moorcroft
# ------------------------------------------------------------------------------------------------------------#
# Code author: Elsa M. Ordway
# ------------------------------------------------------------------------------------------------------------#
# Date: April 28, 2022
# ------------------------------------------------------------------------------------------------------------#

# ------------------------------------------------------------------------------------------------------------#
# load libraries
# ------------------------------------------------------------------------------------------------------------#
library(dplyr); library(tidyr); library(ggplot2); library(viridis)
library(stringr); library(ggfortify); require(data.table)
library(raster); library(rgdal); library(sp); library(GISTools); library(sf)
library(fasterize); library(purrr); library(forcats); library(reshape2)
library(rgeos); library(grid); library(rasterVis);library(ggplotify)
library(cluster); library(ggspatial); library(mapproj); library(outliers)


# ------------------------------------------------------------------------------------------------------------#
# Structural attributes  
# ------------------------------------------------------------------------------------------------------------#
DNM_tch = raster("G:/My Drive/Manuscripts/Communications_Earth_&_Environment/Communications_Revision_Final/Code/Danum_TCH_2m.tif")
DNM_ph  = raster("G:/My Drive/Manuscripts/Communications_Earth_&_Environment/Communications_Revision_Final/Code/Danum_PH_5m.tif")
DNM_lad = brick("G:/My Drive/Manuscripts/Communications_Earth_&_Environment/Communications_Revision_Final/Code/Danum_LAD_50m.tif")
# ------------------------------------------------------------------------------------------------------------#

# ------------------------------------------------------------------------------------------------------------#
# reformat leaf area index / leaf area density (LAD) data
dat_lad <- as.data.frame(DNM_lad, xy = TRUE)
dat_lad <- unite(dat_lad,"xy", c("x","y")); head(dat_lad)
dat_df <- dat_lad[-1] #-c(1:2) excludes x,y columns

# transpose data to long-format
dnm_lad_df <- dat_lad %>% tidyr::gather("band","lad", -c(xy)); head(dnm_lad_df) 
# add height column based on band number
# multiply height * 2 (b/c LAD calculated in 2m bins)

dnm_lad_df$height <- rep(45:0, each=dim(dat_df)[1])*2; head(dnm_lad_df); tail(dnm_lad_df)

# shift/normalize all pixels to ground

# create list of data frames for each pixel
dnm_lad_list <- split(dnm_lad_df, f = dnm_lad_df$xy)

# remove all NAs
# [1] first omit NA rows
dnm_lad_list2 <- lapply(dnm_lad_list, na.omit) 
# [2] exclude empty data frames from list
dnm_lad_list2 <- dnm_lad_list2[sapply(dnm_lad_list2, nrow)>0] 
length(dnm_lad_list); length(dnm_lad_list2)

# use lapply() to remove bands/rows where lad == 0 in each data frame in list
dnm_lad_list_reduced <- lapply(names(dnm_lad_list2), function(x) dnm_lad_list2[[x]][dnm_lad_list2[[x]]$lad != 0, ])
length(dnm_lad_list2); length(dnm_lad_list_reduced)

# remove additional NAs
dnm_lad_list_reduced2 <- lapply(dnm_lad_list_reduced, na.omit) 
dnm_lad_list_reduced2 <- dnm_lad_list_reduced2[sapply(dnm_lad_list_reduced2, nrow)>0] 
length(dnm_lad_list_reduced); length(dnm_lad_list_reduced2)

# count n rows using lapply()
dnm_lad_list_reduced3 <- lapply(dnm_lad_list_reduced2, function(x) transform(x, heightv2 = rev((seq(1,nrow(x)*2, by=2)))))
length(dnm_lad_list_reduced2); length(dnm_lad_list_reduced3) 

# recombine into single data frame
new_dat_dnm <- do.call(rbind, dnm_lad_list_reduced3) # less efficient way of concatenating list elements into df, but works with spatial df's
dim(new_dat_dnm); head(new_dat_dnm)

# calculate peak height of LAI (Hpeak LAI) in meters
# first, get the height at which max LAD occurs for each xy location
result2 <- new_dat_dnm %>% group_by(xy) %>% slice(which.max(lad)) %>% as.data.frame
# separate xy and create raster
result_xy <- result2 %>% separate(xy, c("x","y"))
result_xy$x <- as.numeric(result_xy$x)
result_xy$y <- as.numeric(result_xy$y)
head(result_xy); str(result_xy)
coordinates(result_xy) <- result_xy[,1:2]
projection(result_xy) <- projection(danum_lad)
peakH <- result_xy[,6]
# coerce to gridded dataset
gridded(peakH) <- TRUE

# height at  which peak LAD occurs
peakLAI <- raster(peakH)
plot(peakLAI)
# ------------------------------------------------------------------------------------------------------------#
danum_LAD <- mean(danum_lad)
danum_LAI <- sum((danum_lad*2))

# exclude outliers: LAI values > 15
danum_LAI[danum_LAI > 15] <- NA
danum_LAD_v2 <- mask(danum_LAD, danum_LAI)

DNM_peakLAI <- peakLAI
DNM_LAD_LAI <- stack(danum_LAD_v2, danum_LAI)
DNM_LAD_LAI <- crop(DNM_LAD_LAI, peakLAI, snap="out") 

DNM_LAD_LAI_peakLAI <- stack(DNM_LAD_LAI, peakLAI)

plot(DNM_LAD_LAI_peakLAI)
# ------------------------------------------------------------------------------------------------------------#


# ------------------------------------------------------------------------------------------------------------#
# Canopy foliar traits
# ------------------------------------------------------------------------------------------------------------#
# Leaf Mass per Area (LMA), Leaf nitrogen (N), Leaf phosphorus (P)
# ------------------------------------------------------------------------------------------------------------#
DNM_nplma = brick("G:/My Drive/Manuscripts/Communications_Earth_&_Environment/Communications_Revision_Final/Code/Danum_chems_4m.tif")
DNM_lma    <- DNM_nplma[[1]]
DNM_N_mass <- DNM_nplma[[2]]
DNM_P_mass <- DNM_nplma[[3]]
rm(DNM_nplma)

####  MASK 0's FOR ALL VARIABLES
DNM_tch[DNM_tch <= 0]       <- NA
DNM_N_mass[DNM_N_mass <= 0] <- NA
DNM_P_mass[DNM_P_mass <= 0] <- NA
DNM_lma[DNM_lma <= 0]       <- NA

mask_stack <- stack(DNM_lma, DNM_N_mass, DNM_P_mass)
# ------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
# calculate N:P ratio
DNM_NP <- DNM_N_mass/DNM_P_mass
#----------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
# Combine all into data frame: TCH, N, P, LMA, NP into data frame
#----------------------------------------------------------------------------------------#
# aggregate / resample to different resolutions
#----------------------------------------------------------------------------------------#

# set aggregation factor size (resolution m)
# 2 (4m), 4 (8m), 5 (10m), 10 (20m), 15 (30m), 20 (40m), 25 (50m), 
# 30 (60m), 35 (70m), 40 (80m), 45 (90m), 50 (100m), 60 (120m), 
# 75 (150m), 85 (170m), 100 (200m)
factor_size = 15

tch_resamp     <- aggregate(DNM_tch, fact=factor_size)
res(tch_resamp)

LAD_LAI_resamp <- resample(DNM_LAD_LAI_peakLAI, tch_resamp)
PH_resamp      <- resample(DNM_ph, tch_resamp)
LMA_resamp     <- resample(DNM_lma, tch_resamp)
Nm_resamp      <- resample(DNM_N_mass, tch_resamp)
Pm_resamp      <- resample(DNM_P_mass, tch_resamp)
NP_resamp      <- resample(DNM_NP, tch_resamp)

# calculate max height (99th percentile)
tch_p99 <- aggregate(DNM_tch, fact=factor_size, fun = function(i,...) quantile(i, probs=0.99, na.rm=T))
res(tch_p99) # resample TCH to chem res (4x4) 

# calculate 'cover20' - see Asner et al 2018 Biological Conservation - Sabah carbon paper
# aggregate 2m TCH data to Xm (40) and in process calculate fraction of Xm pixel that's > 20 m
# sum(ifelse(DNM_tch > 20, 1, 0))/n_pixels
cover20 <- aggregate(DNM_tch, fact=factor_size, fun = function(i,...) sum(ifelse(i > 20, 1, 0))/factor_size^2); res(cover20) # resample TCH to chem res (4x4) 
summary(cover20) # should be b/w 0 and 1
#----------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
# stack all
#----------------------------------------------------------------------------------------#
trait_stack_full <- stack(tch_p99, cover20, LAD_LAI_resamp, PH_resamp, LMA_resamp, Nm_resamp, Pm_resamp, NP_resamp)
trait_stack <- mask(trait_stack_full, DNM_50)
trait_stack <- crop(trait_stack, DNM_50, snap="out")

dat_stk <- as.data.frame(trait_stack, xy=T)
colnames(dat_stk) <- c("x","y","maxH","cover20","lad","lai","peakLAI_H","PH","lma","N_mass","P_mass","NP")
summary(dat_stk)

# change digit options for summary() 
options(digits=10)
summary(dat_stk)

#----------------------------------------------------------------------------------------#
# Estimate Vcmax using using the equation in Table 3, model 1 from Walker et al. 2014
#----------------------------------------------------------------------------------------#
# Walker et al. 2014. The relationship of leaf photosynthetic traits–Vcmax and Jmax–to leaf 
# nitrogen, leaf phosphorus, and specific leaf area: a meta‐analysis and modeling study.
# Ecology and evolution, 4(16), pp.3218-3235.
#----------------------------------------------------------------------------------------#
N   <- log(rnorm(3000, 1.75, 0.35)); hist(N)
P   <- log(rnorm(3000, 0.3, 0.05)); hist(P)
N_P <- N*P

vcmax_est <- 3.946 + 0.921*N + 0.121*P + 0.282*N_P
summary(vcmax_est); summary(exp(vcmax_est))

samp_dat <- data.frame(cbind(N,P,N_P, vcmax_est))

# create model object
fit <- glm(vcmax_est ~ N + P + N:P, data = samp_dat); summary(fit)
fit$coef

# site-level Vcmax prediction
pred <- predict(fit, data.frame(N = log(dat_stk$N_area), P = log(dat_stk$P_area)))
dat_stk$vcmax <- exp(pred); summary(dat_stk$vcmax)
dat_stk$vm0 <- dat_stk$vcmax/2.4
summary(dat_stk)

#----------------------------------------------------------------------------------------#
# remove large datasets
#----------------------------------------------------------------------------------------#
rm(danum_lad); rm(danum_LAD); rm(danum_LAI); rm(danum_tch); rm(DNM_LAD_LAI)
rm(danum_p_h); rm(DNM_N_mass); rm(DNM_P_mass); rm(DNM_NP)
rm(DNM_ph);  rm(PH_resamp); rm(LAD_LAI_resamp); rm(DNM_tch)
#----------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------#
# exclude all observations where tch <= 5 m (~16.4 feet)
#----------------------------------------------------------------------------------------#
dat_stk_select_dnm <- subset(dat_stk, maxH > 5)
dat_stk_select_dnm <- subset(dat_stk_select_dnm, NP < 100)

#----------------------------------------------------------------------------------------#
analysis_dat_dnm <- dplyr::select(dat_stk_select_dnm, maxH, cover20, lai, peakLAI_H, PH, lma, N_mass, P_mass, NP, vcmax)
#----------------------------------------------------------------------------------------#
analysis_dat_dnm <- na.omit(analysis_dat_dnm)
#----------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------#
# transform skewed variables 
#----------------------------------------------------------------------------------------#
# right-skewed:sla, NP (log)
hist(dat_stk_select_dnm$NP, breaks=50)
hist(log(dat_stk_select_dnm$NP), breaks=50)
# left-skewed: PH
outlr <- outlier(dat_stk_select_dnm$PH)

analysis_dat_dnm <- subset(analysis_dat_dnm, PH > outlr)

hist(analysis_dat_dnm$PH, breaks=50)
hist(analysis_dat_dnm$PH^(1/3), breaks=50)
hist(log(analysis_dat_dnm$PH), breaks=50)

analysis_dat_dnm$NP <- log(analysis_dat_dnm$NP) 
analysis_dat_dnm$lma <- log(analysis_dat_dnm$lma) 
analysis_dat_dnm$PH <- analysis_dat_dnm$PH^(1/3)

#----------------------------------------------------------------------------------------#
# center & scale data for clustering
#----------------------------------------------------------------------------------------#
kdatS_dnm <- data.frame(apply(analysis_dat_dnm, 2, scale)) 
#----------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
# PCA 
#----------------------------------------------------------------------------------------#
pca_dat_dnm <- prcomp(kdatS_dnm); pca_dat_dnm
summary(pca_dat_dnm)
pca_dat_dnm$rotation

pca_var=pca_dat_dnm$sdev^2
pve=pca_var/sum(pca_var)
cumsum(pve)[2:4] # cumulative variances explained for PCs 2:4
#----------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
# cluster PCA
#----------------------------------------------------------------------------------------#
k = 3
set.seed(89) # ensures initial step in algorithm can be replicated; k-means output fully reproducible
# nstart is the number of random sets 
km.out_dnm <- kmeans(d_pca_dnm, k, nstart=100); km.out_dnm
km.out_dnm$centers

# add x & y back based on row numbers
samp_rows <- rownames(analysis_dat_dnm)
clust_dat_dnm <- bind_cols(d_pca_dnm, dat_stk_select_dnm[samp_rows,c(1,2)])

# plot spatially 
clust_dat_sp_dnm <- SpatialPointsDataFrame(clust_dat_dnm[,c("x", "y")], clust_dat_dnm)
# resample to hectare scale & plot with color gradient
projection(clust_dat_sp_dnm) <- projection(trait_stack)
# create hectare grid (100 x 100 m) for extracting cluster data
grid_1 <- trait_stack[[1]]
grid <- grid_1; res(grid)

# convert cluster point data back to raster 
clust_dat_sp_dnm$cluster <- as.integer(km.out_dnm$cluster)
clust_rast1_dnm <- rasterize(x=clust_dat_sp_dnm, y=grid, field="cluster") 

clust_dat_sp_dnm$cluster <- as.numeric(clust_dat_sp_dnm$cluster)
clust_rast2_dnm <- rasterize(x=clust_dat_sp_dnm, y=grid, fun=mean) 

plot(clust_rast2_dnm$cluster, col=r2)


levelplot(cut(clust_rast2_dnm$cluster, 4), col.regions=r2, 
          at=0:4, margin=FALSE,
          colorkey=list(labels=list(at=0:2 + 0.5, labels=c("Class 1","Class 2","Class 3","Class 4")))) + 
  layer(sp.polygons(DNM_50, lwd=3, col='white'))

#par(mfrow=c(1,2)); plot(clust_rast1_dnm, col=r2); plot(clust_rast2_dnm$cluster, col=r2)
par(mfrow=c(1,2)); plot(clust_rast2_dnm$cluster, col=r2); plot(DNM_50, add=T, border="white"); plot(tch_resamp_DNM, col=r2)
#par(mfrow=c(1,2)); plot(clust_rast2_dnm$cluster, col=r2); plot(tch_resamp, col=r2)

cols2 <- c("#542788","#b2abd2","#fdb863","#b35806")
levelplot(cut(clust_rast2_dnm$cluster, 4), col.regions=cols2, 
          at=0:4, margin=FALSE,
          colorkey=list(labels=list(at=0:2 + 0.5, labels=c("Class 1","Class 2","Class 3","Class 4")))) + 
  layer(sp.polygons(DNM_50, lwd=3, col='black'))

# add "cluster" to datasets and plot above (e.g. biplot)
d_pca_dnm$cluster <- as.factor(km.out_dnm$cluster)
kdatS_dnm$cluster <- as.factor(km.out_dnm$cluster)

#---------------------------------------------------------------------
# PC 1 & 2 colored by CLUSTER
#---------------------------------------------------------------------
autoplot(pca_dat_dnm, data=kdatS_dnm, colour = 'black', 
         size = 2, shape = 21, fill = "cluster", #"P_mass", #"cluster"
         loadings = TRUE, loadings.colour = 'grey60',
         loadings.label = TRUE, loadings.label.size = 5, 
         loadings.label.colour= "black", alpha=0.2,
         loadings.label.label = c("max H","Cover20","LAI","H peak LAI","P:H","LMA","Foliar N","Foliar P","N:P","Vcmax")) + 
  scale_fill_manual(values=cols2) + 
  labs(fill="K-means cluster") + 
  theme(legend.position = "bottom")

#---------------------------------------------------------------------
# Trait distributions by CLUSTER
#---------------------------------------------------------------------
analysis_dat_dnm$cluster <- as.factor(km.out_dnm$cluster)

analysis_dat2_dnm <- subset(analysis_dat_dnm, select=-c(type2))
plot_dat_dnm <- analysis_dat2_dnm %>% gather("var", "value", -cluster)
plot_dat_dnm$cluster <- factor(plot_dat_dnm$cluster, levels=c("2","4","3","1"))

ggplot(plot_dat_dnm, aes(x=cluster, y=value, group=cluster)) + 
  geom_violin(aes(fill=cluster), trim=FALSE) +
  geom_boxplot(width=0.1) +
  scale_fill_manual("k-means cluster", values=c(cols[2],cols[4],cols[3],cols[1]), guide=F) + 
  facet_wrap(~var, scales="free", ncol=3) + 
  scale_x_discrete(labels=c("2" = "2: SPKa", "4" = "4: SPKa", "3" = "3: SPKs", "1" = "4: SPKh")) + 
  labs(x="", y="") 

