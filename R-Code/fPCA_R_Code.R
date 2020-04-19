
# Code here has been taken directly from R-Markdown tutorial files supplied at: https://github.com/johnwarmenhoven/PCA-FPCA. 
# Please go to this site and download the tutorial files for a comprehensive description of this code. 
# This code is designed to provide the "fundamentals" of fPCA and as a consequence reconstruction of components using 
# gganimate is not included in this sample script. This can however be sourced from the Markdown file.

# Authors: John Warmenhoven, Norma Bargary, Dominik Liebl, Andrew Harrison, Mark Robinson, Edward Gunning & Giles Hooker (2020).

# Load fda library and import data and define the knee. 

library(fda)
data(gait)
knee = gait[,,2]

#### FDA Preliminary Steps

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

# Define parameters for smoothing.

kneebasis = create.bspline.basis(rangeval=c(0,100), nbasis=20, norder=4)
Lfdobj = int2Lfd(2)

# Run a generalized cross validation to identify lambda for smoothing and plot the results.

lam         = c(100, 10, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5)
nlam        = length(lam)
GCVsave        = rep(NA,nlam)
Dfsave         = GCVsave;
for(ilam in 1:nlam){
  print(c('lambda = ',lam[ilam]))
  lambda        = lam[ilam];
  fdParobj      = fdPar(kneebasis, Lfdobj, lambda)
  sfd           = smooth.basis(time, knee, fdParobj)
  Dfsave[ilam]  = sfd$df;
  GCVsave[ilam] = sum(sfd$gcv)
}

gcvplot_df = data.frame(lam, GCVsave)
ggplot(data = gcvplot_df, mapping = aes(x=lam, y=GCVsave))+
  geom_line()+
  geom_point(size=3)+
  labs(x=expression(lambda),
       y = expression(GCV(lambda)),
       title= "Lambda Selection")

# Define functional parameter object.

smoothing.parameter = 10;
kneefdPar = fdPar(kneebasis, Lfdobj, smoothing.parameter)

# Define the functional data object.

knee.fd = smooth.basis(time, knee, kneefdPar)

# Define the number of harmonics, set-up a functional parameter object for the fPCs and run the fPCA. 

nharm  = 5
kneePCAfdPar = fdPar(kneebasis, Lfdobj, 1e-15)
knee.pcastr = pca.fd(knee.fd$fd, nharm, kneePCAfdPar)

# Extract important features from the knee.pcastr including the the mean function, fPC functions, proportion of variance for fPCs and scores.

kneemeanfd = mean(knee.fd$fd)
kneeharmfd  = knee.pcastr$harmonics
kneemeanvec = eval.fd(time, kneemeanfd)
kneeharmmat = eval.fd(time, kneeharmfd)
kneevarprop = knee.pcastr$varprop
kneescores = knee.pcastr$scores
stdevfPCscores = apply(kneescores,2,sd)

#### Visualising Results.

# Plotting fPC1.

plotdf_PC1 = data.frame(time=time, mean=kneemeanvec, plus=kneemeanvec + 2*stdevfPCscores[1]*kneeharmmat[,1],minus= kneemeanvec-2*stdevfPCscores[1]*kneeharmmat[,1])
colnames(plotdf_PC1) <- c("time", "mean","plus","minus")

varprop1 = knee.pcastr$varprop[1]

ggplot(data=plotdf_PC1)+
  geom_point(mapping=aes(x=time, y=plus), shape='+', size = 4)+
  geom_point(mapping=aes(x=time, y=minus), shape='-', size = 4)+
  geom_line(mapping=aes(x=time, y=mean), linetype="solid", size = 1)+
  labs(title="fPC1",
       subtitle = paste("Proportion of Variance ", 100*round(varprop1,2), "%", sep=""),
       x="Time (% Gait Cycle)",
       y="Knee (°)")

# Plotting fPC2.

plotdf_PC2 = data.frame(time=time, mean=kneemeanvec, plus=kneemeanvec + 2*stdevfPCscores[2]*kneeharmmat[,2],minus= kneemeanvec-2*stdevfPCscores[2]*kneeharmmat[,2])
colnames(plotdf_PC2) <- c("time", "mean","plus","minus")

varprop2 = knee.pcastr$varprop[2]

ggplot(data=plotdf_PC2)+
  geom_point(mapping=aes(x=time, y=plus), shape='+', size = 4)+
  geom_point(mapping=aes(x=time, y=minus), shape='-', size = 4)+
  geom_line(mapping=aes(x=time, y=mean), linetype="solid", size = 1)+
  labs(title="fPC2",
       subtitle = paste("Proportion of Variance ", 100*round(varprop2,2), "%", sep=""),
       x="Time (% Gait Cycle)",
       y="Knee (°)")

# Plotting fPC1 & fPC2 scores.

plotscores = data.frame(id=id, kneescores)
colnames(plotscores) = c("id","PC1", "PC2", "PC3", "PC4", "PC5")

ggplot(mapping=aes(x=PC1, y=PC2), data=plotscores)+
  geom_point(shape=1, size=3)+
  labs(title = "Scatter Plot of fPC1 & fPC2 Scores",
       x = "fPC1 Scores",
       y = "fPC2 Scores")

# Adding subject id labels to this plot.

library(ggrepel)

ggplot(aes(x=PC1,y=PC2, label = id), data=plotscores)+
  geom_point(shape=1, size=3)+
  geom_label_repel()+
  labs(title = "Scatter Plot of fPC1 & fPC2 Scores",
       x = "fPC1 Scores",
       y = "fPC2 Scores")

#### Visualising rotated components

# Visualising rotated fPC1.

plotdf_PC1_vx = data.frame(time=time, mean=kneemeanvec, plus=kneemeanvec + 2*stdevfPCscores_vx[1]*kneeharmmat_vx[,1],minus= kneemeanvec-2*stdevfPCscores_vx[1]*kneeharmmat_vx[,1])
colnames(plotdf_PC1_vx) <- c("time", "mean","plus","minus")

varprop1_vx = knee.pcastr_vx$varprop[1]

ggplot(data=plotdf_PC1_vx)+
  geom_point(mapping=aes(x=time, y=plus), shape='+', size = 4)+
  geom_point(mapping=aes(x=time, y=minus), shape='-', size = 4)+
  geom_line(mapping=aes(x=time, y=mean), linetype="solid", size = 1)+
  labs(title="fPC1",
       subtitle = paste("Proportion of Variance ", 100*round(varprop1_vx,2), "%", sep=""),
       x="Time (% Gait Cycle)",
       y="Knee (°)")

# Visualising rotated fPC2.

plotdf_PC2_vx = data.frame(time=time, mean=kneemeanvec, plus=kneemeanvec + 2*stdevfPCscores_vx[2]*kneeharmmat_vx[,2],minus= kneemeanvec-2*stdevfPCscores_vx[2]*kneeharmmat_vx[,2])
colnames(plotdf_PC2_vx) <- c("time", "mean","plus","minus")

varprop2_vx = knee.pcastr_vx$varprop[2]

ggplot(data=plotdf_PC2_vx)+
  geom_point(mapping=aes(x=time, y=plus), shape='+', size = 4)+
  geom_point(mapping=aes(x=time, y=minus), shape='-', size = 4)+
  geom_line(mapping=aes(x=time, y=mean), linetype="solid", size = 1)+
  labs(title="fPC1",
       subtitle = paste("Proportion of Variance ", 100*round(varprop2_vx,2), "%", sep=""),
       x="Time (% Gait Cycle)",
       y="Knee (°)")

## Plotting fPC1 & fPC2 scores after rotation.

plotscores_vx = data.frame(id=id, kneescores_vx)
colnames(plotscores_vx) = c("id","PC1", "PC2", "PC3", "PC4", "PC5")

ggplot(mapping=aes(x=PC1, y=PC2), data=plotscores_vx)+
  geom_point(shape=1, size=3)+
  labs(title = "Scatter Plot of fPC1 & fPC2 Scores",
       x = "fPC1 Scores",
       y = "fPC2 Scores")

