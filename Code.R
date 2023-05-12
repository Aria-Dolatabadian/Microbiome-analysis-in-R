#https://vaulot.github.io/tutorials/R_introduction_tutorial.html

library("dplyr")  
library("tidyr")  
library("readxl")
library(magrittr )
library(dplyr)
library("ggplot2") 
library("treemap")
library("phyloseq")
library("FactoMineR")  # For PCA
library("maps")

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("phyloseq")

tara <- read.delim("R_Tara.txt")
str(tara)

#or

#tara <- read_excel("R_Tara.xlsx", sheet = "R Tara")
#str(tara)


#Compute derived quantities and Statistics (using dplyr library)

tara <- tara %>% mutate(Baci_pct = Bacillariophyta/Photo_all * 100, Pela_pct = Pelagophyceae/Photo_all * 
    100)

tara_stat <- tara %>% group_by(fraction, depth_level) %>% summarise(Baci_pct_mean = mean(Baci_pct), 
    Baci_pct_SD = sd(Baci_pct), n = n())
tara_stat

write.csv(tara_stat, "tara_stat.csv", row.names = FALSE)

#Do simple X-Y plots (using ggplot2 library)

#X vs Y
qplot(tara_salinity, tara_temp, data = tara)

#X vs Y with variation in color of points with size fraction
qplot(Baci_pct, Pela_pct, data = tara, color = fraction)

#X vs Y with variation in color of points with size fraction and shape with depth level
qplot(Baci_pct, Pela_pct, data = tara, color = fraction, shape = depth_level)

#X vs Y with variation in color of points with size fraction and shape with depth level
qplot(Baci_pct, Pela_pct, data = tara, color = fraction, shape = depth_level)

#X vs Y with variation sampling_depth for color of points and shape with with size fraction.
qplot(Baci_pct, Pela_pct, data = tara, color = sampling_depth, shape = fraction)

#Categorical data vs y with variation in color of points with depth level
qplot(fraction, Baci_pct, data = tara, color = depth_level)

#Other types of plots

#Boxplot for the same data
qplot(fraction, Baci_pct, data = tara, color = depth_level, geom = "boxplot")

#Histogram for all the data
qplot(Baci_pct, data = tara, geom = "histogram")

#Histogram with different color for each size fraction
qplot(Baci_pct, data = tara, fill = fraction, geom = "histogram")

#Histogram with different graphs (facets) for each size fraction and depth and change bin width

qplot(Baci_pct, data = tara, facets = fraction ~ depth_level, geom = "histogram", 
    binwidth = 10)



#Tree maps (much better than Pie charts…)

tara_tree <- tara %>% select(Sample, depth_level:fraction, Strameno_all, Bacillariophyta:Raphidophyceae) %>% 
    gather(key = Class, value = n_seq, Bacillariophyta:Raphidophyceae)


treemap(tara_tree, index = "Class", vSize = "n_seq", title = "Read numbers")

#Split the tree map according to size fraction

treemap(tara_tree, index = c("fraction", "Class"), vSize = "n_seq", title = "Read numbers")

#Bar graphs

#Only keep surface samples

tara_bar <- tara_tree %>% filter((depth_level == "SUR") & (fraction == "0.8-5"))

ggplot(tara_bar, aes(x = Sample, y = n_seq, fill = Class)) + geom_bar(stat = "identity") + 
    theme_bw() + ggtitle("Tara - Surface - Fraction 0.8-5 µm") + xlab("Samples") + 
    ylab("Number of metabarcodes") + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1))



#Relative abundance

tara_bar <- tara_bar %>% mutate(n_seq_rel = n_seq/Strameno_all)

ggplot(tara_bar, aes(x = Sample, y = n_seq_rel, fill = Class)) + geom_bar(stat = "identity") + 
    theme_bw() + ggtitle("Tara - Surface - Fraction 0.8-5 µm") + xlab("Samples") + 
    ylab("Fraction of metabarcodes") + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1))


#Heat maps

tara_heat <- tara %>% filter(fraction == "0.8-5") %>% select(Bacillariophyta:Raphidophyceae)
tara_heat.matrix <- data.matrix(tara_heat)

heatmap(tara_heat.matrix, margins = c(20, 6))


#Multivariate analysis (FactoMiner package)


# Select only the 0.8-5 µm fraction and only the colums with phytplankon
# data and metadata
tara_multi <- tara %>% filter(fraction == "0.8-5")

# Define row names as 'Station_Depth level' (points with be labelled by row
# names)
row.names(tara_multi) <- paste(tara_multi$station, tara_multi$depth_level, sep = "_")

# Select only with phytoplankon data and metadata
tara_multi <- tara_multi %>% select(Bacillariophyta:Raphidophyceae, chloro_hplc:tara_salinity)

# Scale the matrix
tara_multi <- scale(tara_multi)

# Do the PCA
tara_pca <- PCA(tara_multi)

#Maps

tara_map <- tara %>% filter((fraction == "0.8-5") & (depth_level == "SUR"))
# Draw the world map
map(database = "world", fill = TRUE)

# Add stations
points(tara_map$Longitude, tara_map$Latitude, pch = 3, col = "red", cex = 1)

# Add data - circle size is proprotional to proportion of
points(tara_map$Longitude, tara_map$Latitude, pch = 19, col = "blue", cex = tara_map$Baci_pct * 
    3/100)

# Add title
title("Bacilliorophyta as % of Photosynthetic - 0.8-5 µm - surface", cex.main = 1)
















