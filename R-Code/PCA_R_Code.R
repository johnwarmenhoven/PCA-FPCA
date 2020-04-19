# Code here has been taken directly from R-Markdown tutorial files supplied at: https://github.com/johnwarmenhoven/PCA-FPCA. 
# Please go to this site and download the tutorial files for a comprehensive description of this code. 
# This code is designed to provide the "fundamentals" of PCA applied to waveform data. 

# Authors: John Warmenhoven, Norma Bargary, Dominik Liebl, Andrew Harrison, Mark Robinson, Edward Gunning & Giles Hooker (2020).

# Load fda library and import data and define and then transpose the knee in preparation for PCA. 

library(fda)
data(gait)
knee = gait[,,2]
Knee_Curves = t(knee)

#### Running PCA of Waveforms

# Define time interval and subject identifiers. 

time = seq(0, 100,length.out =  20)
id = dimnames(gait)[[2]]

preplot_df = data.frame(time, knee)
colnames(preplot_df) = c("time", id)

# Re-shape in preparation to inspect data.

library(reshape2) 
preplot_df.mlt = melt(data=preplot_df, id.vars = "time", value.name = "knee", variable.name = "id")

# Plot raw data.

library(ggplot2)
theme_set(theme_classic())
theme_update(plot.title = element_text(hjust = 0.5, size=20, face="bold"),
             plot.subtitle = element_text(hjust = 0.5, size=16),
             axis.title= element_text(size=16, face="bold"),
             axis.text = element_text(size = 14))

ggplot(data=preplot_df.mlt)+ 
  geom_line(mapping=aes(x=time, y=knee, group=id, color=id))+
  theme(legend.position = "none")+
  labs(x = "Time (% Gait Cycle)",
       y="Knee (°)",
       title = "Gait Data",
       subtitle = "Knee Angle for 39 repetitions of a 20 point gait cycle")

# princomp function and applying PCA. princomp is an inbuilt R function.

ClassicPCA <- princomp(t(knee),cor=FALSE,scores=TRUE)

# Extract the mean knee curve, calculate the standard deviation for PC scores (for each component) and calculate the percentage of variation for each component. 
# Each of these will be used in plotting results. 

mean_curve = apply(Knee_Curves,2,mean)
stdevscores = apply(ClassicPCA$scores,2,sd)
PoV <- ClassicPCA$sdev^2/sum(ClassicPCA$sdev^2)

#### Visualising Results.

# Plotting PC1.

plotdf_PC1 = data.frame(time=time, mean=mean_curve, plus=mean_curve + 2*stdevscores[1]*ClassicPCA$loadings[,1],minus= mean_curve-2*stdevscores[1]*ClassicPCA$loadings[,1])
colnames(plotdf_PC1) <- c("time", "mean","plus","minus")

ggplot(data=plotdf_PC1)+
  geom_point(mapping=aes(x=time, y=plus), shape='+', size = 4)+
  geom_point(mapping=aes(x=time, y=minus), shape='-', size = 4)+
  geom_line(mapping=aes(x=time, y=mean), linetype="solid", size = 1)+
  labs(title="PC1",
       subtitle = paste("Proportion of Variance ", 100*round(PoV[1],2), "%", sep=""),
       x="Time (% Gait Cycle)",
       y="Knee (°)")

# Plotting PC2.

plotdf_PC2 = data.frame(time=time, mean=mean_curve, plus=mean_curve + 2*stdevscores[2]*ClassicPCA$loadings[,2],minus= mean_curve-2*stdevscores[2]*ClassicPCA$loadings[,2])
colnames(plotdf_PC2) <- c("time", "mean","plus","minus")

ggplot(data=plotdf_PC2)+
  geom_point(mapping=aes(x=time, y=plus), shape='+', size = 4)+
  geom_point(mapping=aes(x=time, y=minus), shape='-', size = 4)+
  geom_line(mapping=aes(x=time, y=mean), linetype="solid", size = 1)+
  labs(title="PC2",
       subtitle = paste("Proportion of Variance ", 100*round(PoV[2],2), "%", sep=""),
       x="Time (% Gait Cycle)",
       y="Knee (°)")

# Plotting PC1 & PC2 scores.

plotscores = data.frame(id=id, ClassicPCA$scores)
colnames(plotscores) = c("id","PC1", "PC2", "PC3", "PC4", "PC5")

ggplot(mapping=aes(x=PC1, y=PC2), data=plotscores)+
  geom_point(shape=1, size=3)+
  labs(title = "Scatter Plot of PC1 & PC2 Scores",
       x = "PC1 Scores",
       y = "PC2 Scores")

# Adding subject id labels to this plot.

library(ggrepel)

ggplot(aes(x=PC1,y=PC2, label = id), data=plotscores)+
  geom_point(shape=1, size=3)+
  geom_label_repel()+
  labs(title = "Scatter Plot of PC1 & PC2 Scores",
       x = "PC1 Scores",
       y = "PC2 Scores")
