%% A Tutorial for Applying PCA of Waveforms in Biomechanics.
%
% John Warmenhoven, Norma Bargary, Dominik Liebl, Andrew Harrison, Mark 
% Robinson, Edward Gunning & Giles Hooker

%% Load sample data
%
% Data for this tutorial can be taken from the Functional Data Analysis
% website: http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/
% (once the fdaM.zip folder is downloaded, data is located at
% 'fdaM\examples\gait')
%
% For any questions related to this tutorial (and script), please email
% john.warmenhoven@hotmail.com.
%
% The data can be loaded below (make sure this is in the same folder as
% your m-file).

%%
load('gait');

%%
% This loads open access data collected at the Motion Analysis Laboratory
% at Children’s Hospital in California. Full details on the collection of
% the data can be found at Olshen, et al. (1989).
%
% The data consists of hip and knee joint curves of 39 children.
% Each curve consists of 20 data points across the gait cycle. This
% tutorial will only focus on the knee joint.

%% Principal Components Analysis (PCA) of Waveforms
%
% Define the knee joint curves and transpose them in preparation for PCA.

%%
Knee_Curves = knee';

%%
% The knees have been normalised to 20 data points using an interpolating
% cublic spline. Here they are centered within the time interval (which in
% this case after temporal normalization is percentage of the gait cycle).

%%
time = linspace(0, 100, 20)'; 

%%
% The 'pca' function is an inbuilt Matlab function and can be applied to
% waveform data as a means of dimension reduction. This function returns
% co-efficients (x), which are principal components (PC), PC scores (y), 
% eigenvalues (z) of the covariance matrix of 'Knee_Curves,' Hotelling's 
% T-squared statistic for each observation (ii), and the percentage of the
% total variation attributed to each PC (iii). 

%%
[x,y,z,ii,iii] = pca(Knee_Curves);

%%
% Within the context of PCA of waveforms, the principal components (PC)
% themselves are time-series and can be observed graphically. See the
% figures below. Figure(1) shows PC1 and figure(2) shows PC2. 

%%
figure(1)
plot(time, x(:,1), 'k');
hold on
plot(time, zeros(20,1), 'k--')
hold off
axis([0 100 -0.3 0.5])

figure(2)
plot(time, x(:,2), 'k');
hold on
plot(time, zeros(20,1), 'k--')
hold off
axis([0 100 -0.3 0.5])

%%
% The PCs can also be observed from a more practical perspective by adding
% and subtracting them from the average/mean curve. To do this, first the
% mean curve is calculated:

%%
mean_curve = mean(Knee_Curves);

%% Visualising Results
%
% Then a constant is created to scale each of the PCs before plotting them
% relative to the mean. A constant that is commonly used is 1 or 2
% standard deviations (SD) of the PC scores. So an SD of PC scores is
% calculated for each PC. 

%%
stdevscores = std(y);

%%
% Then this is used to create the constant.

%%
con1 = stdevscores(1)*1;

%%
% And then scaled PCs are added and subtracted from the mean curve, before
% being plotted with the mean curve. 

%%
PC1_Pos = mean_curve' + con1.*(x(:,1));
PC1_Neg = mean_curve' - con1.*(x(:,1));

%%
% Figure(3) shows PC1, with positive scorers plotted using the '+'
% symbols and negative scorers plotted using the '-' symbols. 

%%
figure(3)
phdl = plot(time, mean_curve', 'k-');
set(phdl, 'LineWidth', 1)
hold on
phdl = text(time, PC1_Pos, '+');
set(phdl, 'LineWidth', 1)
hold on
phdl = text(time, PC1_Neg, '-');
set(phdl, 'LineWidth', 1)
hold off
box on
set(gca,'FontSize',10);
title('\fontsize{10} PC1');
xlabel('\fontsize{10} Gait Cycle (%)')
ylabel('\fontsize{10} Knee Angle (°)')

%%
% The same process is also carried out for PC2. Inclusive of creating a
% constant to scale the second PC.

%%
con2 = stdevscores(2)*1;

PC2_Pos = mean_curve' + con2.*(x(:,2));
PC2_Neg = mean_curve' - con2.*(x(:,2));

figure(4)
phdl = plot(time, mean_curve', 'k-');
set(phdl, 'LineWidth', 1)
hold on
phdl = text(time, PC2_Pos, '+');
set(phdl, 'LineWidth', 1)
hold on
phdl = text(time, PC2_Neg, '-');
set(phdl, 'LineWidth', 1)
hold off
box on
set(gca,'FontSize',10);
title('\fontsize{10} PC2');
xlabel('\fontsize{10} Gait Cycle (%)')
ylabel('\fontsize{10} Knee Angle (°)')

%%
% An example of a way to graph the principal component scores can also be 
% seen below, plotting scores for one PC1 relative to PC2. 

figure(5)
scatter(y(:,1),y(:,2), 'k');
box on
set(gca,'FontSize',10);
title('\fontsize{10} PC1 vs PC2');
xlabel('\fontsize{10} PC1 Scores')
ylabel('\fontsize{10} PC2 Scores')

%% References
% 
% Olshen, R. A., Biden, E. N., Wyatt, M. P. & Sutherland, D. H. (1989).
% Gait analysis and the bootstrap. _The Annals of Statistics_, 1419-1440.
