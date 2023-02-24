################################################################################
############ POPULATION STRUCTURE ANALYSIS IN Orcinus orca ECOTYPES ############ 
################################################################################

library(vcfR) # reads vcf file
library(Rtsne) # computes t-sne
library(cluster) # clusters using kmeans
library(mclust) # adjusted rand index
library(ggplot2) # plotting
library(ggpubr) # plotting
library(maps) # plotting
library(dplyr) # data manipulation
library(caTools) # splitting the data set
library(randomForest) # create random forest
library(MLeval) # ROC curve
library(caret) # confusion matrix


##### INPUT DATA
vcf_file <- "snps_o.orca.vcf"

# READ IN THE DATA using vcfR
vcf <- read.vcfR(vcf_file)

##### Create a map plot indicating sampling locations
palette <- c("#CC0000","#EDADC7","#414040","#0080FF",
                      "#FF8000", "#B2FF66", "#EEC643",
                      "#00FFFF","#CC00CC")
                      
# Create data frame with location data
coordinates <- data.frame(
  location = c("AR", "SR", "RU", "BS", "AT", "CT", "OS", "IC", "MI"),
  lat = c(52.568369,45.188298,53.414827,50.647426,58.886533,
          40.660677,42.863116,59.782961,-46.655438),
  lon = c(-133.398362,-125.551147,164.682160,-170.526769,-149.444446,
          -127.776410,-152.372453,-22.769677,54.955717),
  palette = palette)

# Plot the world map
world_coordinates <- map_data("world")
map_plot <- ggplot() + theme_void() +
  geom_map(data = world_coordinates, 
           map = world_coordinates,
           aes(long, lat, map_id = region),
           fill = "lightgrey") +
  #geom_sf() +
  #coord_sf(ylim = c(-80, 90), xlim = c(-185,195), expand = FALSE) +
  geom_point(data=coordinates, aes(lon,lat,color=location), size=6) +
  scale_color_manual(values=coordinates$palette,
                     breaks=c("AR", "SR","RU","BS","AT","CT","OS","IC","MI"))+
  theme(legend.position="none", panel.background = element_rect(fill = "white")) +
  geom_text(data=coordinates, aes(lon, lat, label=location), hjust=-.7, vjust=.4,
            fontface="bold", size=4)




#################################################################
######################## PRE-PROCESSING #########################
#################################################################

##### Create a data frame with the extracted SNP genotypes in the VCF file
geno_df <- as.data.frame(extract.gt(vcf))
# Transform the data to binary according to the genotype (individuals are either 0 or 1 for a SNP)
geno_df <- ifelse(geno_df == "0/0", 0, 1)
# Transpose the data so samples are in rows and features (SNPs) are in columns
geno <- as.data.frame(t(geno_df))


##### Create the metadata variables
# There is no available metadata table so this data is manually reconstructed
# following the info in the Supplementary Materials of doi: 10.1111/mec.12929 (Hoelzel, 2014)
# The resulting metadata table contains both character and numeric versions of the vars

# Ecotype of every sample (5 categories)
ecotype <- rep(c("Resident","Offshore","Resident","Transient", "Unknwon","Offshore",
                 "Resident","Transient","Antartic Type B"),
               times=c(22,1,8,36,6,6,22,1,13))
num.ecotype <- as.numeric(as.factor(ecotype))

# Population codes (9 categories)
popcode <- rep(c("AR","SR","OS","SR","AT","CT","IC","OS","RU","BS","CT","MI"),
               times=c(17,5,1,8,21,15,6,6,9,13,1,13))
num.popcode <- as.numeric(as.factor(popcode))

# Sample location (9 categories)
location <- rep(c("Alaska Residents (AR)", "Washington Residents (SR)", "Pacific Offshores (OS)","Washington Residents (SR)", 
                  "Alaska Transients (AT)", "California Transients (CT)", "Iceland (IC)","Pacific Offshores (OS)", 
                  "Kamchatka Residents (RU)", "Bering Residents (BS)", "California Transients (CT)","Marion Island Antartics (MI)"), 
                times=c(17,5,1,8,21,15,6,6,9,13,1,13))
num.location <- as.numeric(as.factor(location))

# Create a dataframe with the metadata variables
metadata <- data.frame(popcode,ecotype,location, num.popcode, num.ecotype, num.location)




############################
### UNSUPERVISED METHODS ###
############################

##### Let's reduce dimensionality
# perform principal component analysis on the given data and return an object of class prcomp
pcaout <- prcomp(geno) 
summary(pcaout) # PC1 accounts for 17.7% of variance, PC2 accounts for 6.6%
# visualize that PC1 and 2 are indeed the main contributors to variance
screeplot(pcaout, main="Contribution of every PCA component", col="red") 
# create a dataframe with the loadings
pcadata <- pcaout$x # value of the centred data multiplied by the rotation matrix
pca <- as.data.frame(pcadata)

# Perform t-SNE. Value for Perplexity is the maximum allowed for this number of samples
tsneout <- Rtsne(geno, PCA=FALSE, perplexity=38, theta=0, check_duplicates=FALSE)
# create a dataframe with the principal values
tsne <- as.data.frame(tsneout$Y)


##### Plot the results

# custom palette for all pops

# pca 
pca_plot <- ggplot(pca, aes(x=PC1, y=PC2, color=metadata$location)) +
  geom_point(size=3) + ggtitle("PCA") + theme_bw() + 
  scale_color_manual(values=palette, 
                     breaks=c("Alaska Residents (AR)", "Washington Residents (SR)",
                              "Kamchatka Residents (RU)", "Bering Residents (BS)",
                              "Alaska Transients (AT)", "California Transients (CT)",
                              "Pacific Offshores (OS)", "Iceland (IC)", 
                              "Marion Island Antartics (MI)")) + 
  theme(legend.position="none",legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# t-sne 
tsne_plot <- ggplot(tsne, aes(x=V1, y=V2, color=metadata$location)) +
  geom_point(size=3) + ggtitle("t-SNE") + theme_bw() + 
  scale_color_manual(values=palette, 
                     breaks=c("Alaska Residents (AR)", "Washington Residents (SR)",
                              "Kamchatka Residents (RU)", "Bering Residents (BS)",
                              "Alaska Transients (AT)", "California Transients (CT)",
                              "Pacific Offshores (OS)", "Iceland (IC)", 
                              "Marion Island Antartics (MI)")) + 
  theme(legend.position="none",legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# combine
ggarrange(pca_plot, tsne_plot, labels= c("A","B"), ncol=2, widths=c(2,2),
          common.legend = TRUE,legend="top") 


##### Clustering
# Optimal number of clusters will not be estimated due to the research question
# Do clustering methods reveal a genetic structure based on population or ecotype?

# K-means clustering on the pca data with k= 5(ecotypes) and k=9(pops)
kmeans_k5 <- kmeans(pca, centers= 5, nstart = 100)
kmeans_k9 <- kmeans(pca, centers= 9, nstart = 100)
# add the results to the metadata info
metadata$kmeans_k5 <- as.factor(kmeans_k5$cluster)
metadata$kmeans_k9 <- as.factor(kmeans_k9$cluster)

# Hierarchical clustering on the pca data 
# calculate the distance matrix using euclidean distance
dist_matrix <- dist(pca, method = "euclidean")
# hierarchical clutering 
hc <- hclust(dist_matrix, method="ward.D2")
# cut the model in k clusters as with kmeans
hc_k5 <- cutree(hc, k=5)
hc_k9 <- cutree(hc, k=9)


##### Plot the results
# pca 
pca_plot <- ggplot(pca, aes(x=PC1, y=PC2, color=metadata$popcode)) +
  geom_point(size=3) + ggtitle("PCA") + theme_bw() + 
  scale_color_manual(values=palette, 
                     breaks=c("AR", "SR","RU","BS","AT","CT","OS","IC","MI")) +
  theme(legend.title = element_blank(), legend.justification=c(0,1), legend.position=c(0,1),
        legend.background = element_rect(color="black", fill=NA),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# kmeans k=5 plot
kmeans_k5_plot <- ggplot(pca, aes(x=PC1, y=PC2, color=as.factor(kmeans_k5$cluster))) +
  geom_point(size=3) + ggtitle("K-Means k = 5") + theme_bw() +
  geom_point(data = as.data.frame(kmeans_k5$centers), aes(x=PC1,y=PC2), color="black", shape = 5, size = 3) +
  theme(legend.title = element_blank(), legend.justification=c(0,1), legend.position=c(0,1),
        legend.background = element_rect(color="black", fill=NA),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# kmeans k=9 plot
kmeans_k9_plot <- ggplot(pca, aes(x=PC1, y=PC2, color=as.factor(kmeans_k9$cluster))) +
  geom_point(size=3) + ggtitle("K-Means k = 9") + theme_bw() + scale_color_brewer(palette="Paired") +
  geom_point(data = as.data.frame(kmeans_k9$centers), aes(x=PC1,y=PC2), color="black", shape = 5, size = 3) +
  theme(legend.title = element_blank(), legend.justification=c(0,1), legend.position=c(0,1),
        legend.background = element_rect(color="black", fill=NA),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# hc k=5 plot
hc_k5_plot <- ggplot(pca, aes(x=PC1, y=PC2, color=as.factor(hc_k5))) +
  geom_point(size=3) + ggtitle("Hierarchical Clustering k = 5") + theme_bw() +
  geom_point(data = as.data.frame(kmeans_k5$centers), aes(x=PC1,y=PC2), color="black", shape = 5, size = 3) +
  theme(legend.title = element_blank(), legend.justification=c(0,1), legend.position=c(0,1),
        legend.background = element_rect(color="black", fill=NA),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# hc k=9 plot
hc_k9_plot <- ggplot(pca, aes(x=PC1, y=PC2, color=as.factor(hc_k9))) +
  geom_point(size=3) + ggtitle("Hierarchical Clustering k = 9") + theme_bw() + scale_color_brewer(palette="Paired") +
  geom_point(data = as.data.frame(kmeans_k9$centers), aes(x=PC1,y=PC2), color="black", shape = 5, size = 3) +
  theme(legend.title = element_blank(), legend.justification=c(0,1), legend.position=c(0,1),
        legend.background = element_rect(color="black", fill=NA),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# combine
ggarrange(kmeans_k5_plot,kmeans_k9_plot, hc_k5_plot, hc_k9_plot,
                    ncol=2, nrow=2)


#### Confusion matrices
kmeans_k5_conf <- addmargins(table(metadata$ecotype, kmeans_k5$cluster))
kmeans_k9_conf <- addmargins(table(metadata$popcode, kmeans_k9$cluster))
hc_k5_conf <- addmargins(table(metadata$ecotype, hc_k5))
hc_k9_conf <- addmargins(table(metadata$popcode, hc_k9))

#### Get the Adjusted Rand Index using k=5 as reference
k5_ari <- adjustedRandIndex(kmeans_k5$cluster, hc_k5)
k9_ari <- adjustedRandIndex(kmeans_k9$cluster, hc_k9)
print(k5_ari) # 1
print(k9_ari) # 0.89





############################
#### SUPERVISED METHODS ####
############################
# Let's see how good a random forest can perform in order to
# distinguish the samples in different populations 

# add population and ecotype  data to the pca data
pca_rf <- pca
pca_rf$population <- as.factor(popcode)


# Split the data in training (0.8) and testing (0.2) datasets
split <- sample.split(pca_rf, SplitRatio = 0.8)
data_train <- na.omit(subset(pca_rf, split == TRUE))
data_test <- na.omit(subset(pca_rf, split== FALSE))

# Find the best value of random variables
pop_mtry <- tuneRF(data_train, data_train$population,
                   stepFactor = 1.2, improve = 0.01, trace=T, plot= T) 

rf <- randomForest(population~., data=data_train, ntree=20000, 
                   importance=TRUE, mtry=9)

rf # summary and confusion matrix
#importance(rf) # importance of each variable
#varImpPlot(rf) # visualize importance of each variable, PC2 and PC1 are the most relevant

# Make predictions on test data
pred_test <- predict(rf, newdata=data_test, type="class")
pred_test
pred_test_prob <- predict(rf, newdata=data_test, type="prob")[,2]


# Confusion matrix
conf_rf <- confusionMatrix(table(pred_test,data_test$population))
conf_rf

# AUC of ROC curve
roc <- multiclass.roc(response=data_test$population, 
           predictor=pred_test_prob)
auc(roc) #Multi-class area under the curve: 0.9366
