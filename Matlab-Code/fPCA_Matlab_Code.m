%% A tutorial on applying Functional Principal Components Analysis (fPCA) in Matlab
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
% All content from this tutorial can be found in two key texts.
% All theoretical information underpinning the use of these processes can
% be found at Ramsay and Silverman (2005). Additionally, a useful text
% outlining computational approaches for applying these processes in
% Matlab and R can be found at Ramsay, Hooker and Graves (2009).
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

%% Add file paths
%
% The files for using FDA can be obtained from the software/download
% sections of the www.functionaldata.org website. The title of the folder
% containing the Matlab functions is "fdaM."
%
% The software packages for Matlab are available and can be downloaded to
% a convenient location on your computer. The paths of these Matlab
% functions will need to be added before using these techniques.
%
% Make sure to name the filepath correctly (see below):

%%
addpath('E:\fdaM\')

%% Functional Data Analysis (Preliminary Steps)
%
% FDA involves a series of preliminary steps prior to applying any
% techniques. These are mainly centred around function fitting, smoothing
% processes and registration of data.
%
% An example of function fitting and smoothing is included in this script,
% but registration was not necessary for this data set so was not
% included.
%
% The knees have been normalised to 20 data points using an interpolating
% cublic spline. Here they are centered within the time interval (which in
% this case after temporal normalization is percentage of the gait
% cycle).

%%
time = linspace(0, 100, 20)'; 

%%
% The next step involves representation of functions using a suitable
% basis expansion. This involves representing the observed curves as a
% linear (or weighted) combination of known functions (basis functions),
% where the coefficients (or weights) are chosen from the data.
%
% FDA has options for a number of different bases (i.e. Fourier,
% B-splines, Wavelets, etc.). In this demonstration the functions were
% fitted using B-splines, which are commonly drawn upon for biomechanics
% data. Spline functions are the most common choice of expansion for
% non-periodic functional data, which constitutes many types of human
% movements (i.e. movements that are not necessarily cyclical).
%
% To set up for function fitting with B-splines, a fourth-order spline
% (cubic spline) is selected, the maximum number of functions is selected
% (at 20 basis functions), given that a smoothing parameter will be added
% later to the function fitting process. Smoothness of the data can also
% be manipulated by the number of basis functions selected (less basis
% functions results in a smoother fit of the functions).
% 
% A cubic spline (spline of order 4) is selected as these are considered
% sufficient to approximate nearly all sufficiently smooth functions (with
% the knee angles in this tutorial being one example of this). If we were
% interested in velocity or acceleration, we may consider selecting a
% higher order spline (5th and 6th respectively).

%%
kneebasis = create_bspline_basis([0,100], 20, 4); 

%%
% When choosing the coefficients for the basis representation we want to 
% penalise both fit to the observed data and roughness of our estimated
% function. Below we set up the penalty for roughness by defining a linear
% differential operator. This defines how we will smooth our data. If we
% want the original data to be smooth we penalise the second derivative,
% which is related directly to the curvature of our original data (and we
% do this by selecting 2 within the |int2Lfd| function). In this sense,
% |Lfdobj| is defining “what” is smooth.

%%
Lfdobj = int2Lfd(2); 

%%
% If we wanted to focus on the smoothness of our first derivative (knee
% angular velocity in this instance) we would penalize the third
% derivative. Similarly knee angular acceleration smoothness would require
% penalizing the fourth derivative.
%
% We have defined how we plan to smooth the data, but before we can finish
% defining a functional parameter object, generalized cross validation
% (GCV) can be used to trial values of lambda (a smoothing parameter), to
% determine a suitable level of smoothing.
%
% Lambda, also known as the roughness penalty parameter, controls the
% trade-off between the fit to the data and the roughness penalty we have
% defined. A large value of lambda would put more emphasis on the
% roughness penalty and result in a smoother estimated function. A smaller
% value would conversely put more emphasis on fit to the data and result
% in a less smooth estimated function.
%
% The GCV criterion (defined within this loop as |GCVsave|) is a measure
% of the predictive error for different values of lambda, with this often
% being trialled across a range of potential values. This is the most
% common way to get an estimate for an appropriate value of lambda and was
% first defined by Craven and Wahba (1979).
%
% Below is a for-loop for calculating GCV at a range of different values
% of lambda, varying between 100 and 1e-5 and being represented by |lam|.

%%
lam         = [100, 10, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5];
nlam        = length(lam);
Gcvsave        = zeros(nlam,1);
Dfsave         = Gcvsave;
for ilam = 1:nlam
     lambda        = lam(ilam);
    fdParobj      = fdPar(kneebasis, Lfdobj, lambda);
    [logknee_fd, logknee_df, logknee_gcv] = ...
                    smooth_basis(time, knee, fdParobj);
    Dfsave(ilam)  = logknee_df;
    Gcvsave(ilam) = sum(logknee_gcv);
end

%%
% Here is a plot of the GCV criterion at different values of lambda.

%%
figure(1);
title('Lambda Selection');
phdl=plot(lam, Gcvsave, 'ko-'); 
set(phdl, 'LineWidth', 1)
xlabel('\fontsize{13} \lambda')
ylabel('\fontsize{13} GCV(\lambda)')

%%
% Given that the GCV criterion relates directly to the predictive error
% for different values of lambda, a good starting point for selecting a
% value for lambda is the smallest GCV value trialled in the minimization
% routine (which in this case was 10). This may not always work perfectly
% in practice and should always be confirmed through visual inspection
% (i.e. graphing the fitted data).
%
% Now we can complete defining a functional parameter object smoothing and
% fitting the functions to the data. The results of GCV have suggested
% that a smoothing parameter of 10 may be suitable for this data.

%%
smoothing_parameter = 10;

kneefdPar = fdPar(kneebasis, Lfdobj, smoothing_parameter);

%%
% The |fdPar| function serves as a way of capturing all the information
% required for function fitting and smoothing. This is inclusive of the
% expansion being used (in this case B-splines), what is being focused on
% for smoothing (in this case penalising the second derivative) and
% roughness penalty in the form of lambda (in this case a value of 10). 

%% Functional Principal Components Analysis (fPCA)
%
% In fPCA we store all smoothed curves in the same functional data object
% using the |smooth.basis| function. 

%%
knee_fd = smooth_basis(time, knee, kneefdPar);

%%
% We have also listed the number of functional principal components to be
% retained as five.

%%
nharm  = 5;

%%
% Similar to the function fitting and smoothing processes described in the
% FDA (Preliminary Steps) section, it is also possible to smooth fPC
% functions as a part of describing them. This time a negligible
% smoothing parameter was selected (a very minimal penalty).

%%
kneePCAfdPar = fdPar(kneebasis, Lfdobj, 1e-15);
knee_pcastr = pca_fd(knee_fd, nharm, kneePCAfdPar);

%%
% If researchers wish to have more control over graphing of the fPCs
% (similar to the example provided for conventional PCA), relevant parts 
% of the FDA and fPCA processes can be extracted and called on for
% plotting. 
%
% Most of these will come from a 'struct' that has been built as a part of 
% the |pca_fd| function. In this case the struct is called |knee_pcastr|.
%
% As an example the mean function (as a functional data object) can be 
% called on using the following.

%%
kneemeanfd = mean(knee_fd);

%%
% Similarly, the fPC functions (as functional data objects) can be called
% upon using the following.

%%
kneeharmfd  = knee_pcastr.harmfd;

%%
% We can also unwrap the mean function and fPC functions to a vector of
% points that can be used for graphical observation and plotting of 
% results.

%%
kneemeanvec = squeeze(eval_fd((time), kneemeanfd));
kneeharmmat = eval_fd(time, kneeharmfd);

%%
% We can also identify the amount of variation attributed to each fPC by 
% exploring the |varprop| part of the struct. 

%%
kneevarprop = knee_pcastr.varprop;

%%
% And also derive the weights attributed to each of the individual curves
% relative to each fPC. These are also referred to as fPC scores (similar
% to PC scores in the previous example). 

%%
kneescores = knee_pcastr.harmscr; 

%% Visualising Results
%
% Similar to PCA, a constant can be created to scale fPCs before adding 
% and substracting them from the mean. A constant that is commonly used 
% is 1 or 2 standard deviations (SD) of the fPC scores, with these being
% calculated for each fPC. 

%%
stdevfPCscores = std(kneescores);

%%
% After this 1 SD is used to create the constant:

%%
con1 = stdevfPCscores(1)*1;

%%
% And then scaled these fPCs are added and subtracted from the mean curve, 
% to show positive and negative scoring deviations away from the mean for
% each fPC.

%%
fPC1_Pos = kneemeanvec + con1.*(kneeharmmat(:,1));
fPC1_Neg = kneemeanvec - con1.*(kneeharmmat(:,1));

%%
% |figure(2)| shows fPC1, with positive scorers plotted using the '+'
% symbols and negative scorers plotted using the '-' symbols.

%%
figure(2)
phdl = plot(time, kneemeanvec, 'k-');
set(phdl, 'LineWidth', 1)
hold on
phdl = text(time, fPC1_Pos, '+');
set(phdl, 'LineWidth', 1)
hold on
phdl = text(time, fPC1_Neg, '-');
set(phdl, 'LineWidth', 1)
hold off
box on
set(gca,'FontSize',10);
title('\fontsize{10} fPC1');
xlabel('\fontsize{10} Gait Cycle (%)')
ylabel('\fontsize{10} Knee Angle (°)')

%%
% The same process is also carried out for fPC2. Inclusive of creating a
% constant to scale the fPCs.

%%
con2 = stdevfPCscores(2)*1;

fPC2_Pos = kneemeanvec + con2.*(kneeharmmat(:,2));
fPC2_Neg = kneemeanvec - con2.*(kneeharmmat(:,2));

figure(3)
phdl = plot(time, kneemeanvec, 'k-');
set(phdl, 'LineWidth', 1)
hold on
phdl = text(time, fPC2_Pos, '+');
set(phdl, 'LineWidth', 1)
hold on
phdl = text(time, fPC2_Neg, '-');
set(phdl, 'LineWidth', 1)
hold off
box on
set(gca,'FontSize',10);
title('\fontsize{10} fPC2');
xlabel('\fontsize{10} Gait Cycle (%)')
ylabel('\fontsize{10} Knee Angle (°)')

%%
% An example of a way to graph the principal component scores can also be 
% seen below, plotting to scores for one fPC relative to another fPC. 

%%
figure(4)
phdl = scatter(kneescores(:,1),kneescores(:,2), 'k');
set(phdl, 'LineWidth', 1)
box on
set(gca,'FontSize',10);
title('\fontsize{10} fPC1 vs fPC2');
xlabel('\fontsize{10} fPC1 Scores')
ylabel('\fontsize{10} fPC2 Scores')

%% Varimax Rotations
%
% Varimax rotations are used to construct new components based on original
% principal components obtained from the above process of using fPCA.
% Varimax rotations maximize the variability of the squared principal
% component weights, for a selected group of fPCs. The resulting modes of
% variability tend to be concentrated on a part of the range of the
% function in question. Generally this focuses more acutely on areas of
% the curve within original fPCs, making the results sometimes more easier
% to interpret.
%
% First we use the |varmx_pca| function, entering the |knee.pcastr| struct
% (from the original fPCA) as an input. |knee.pcastr_vx| is the outputted 
% as a struct with rotated components and recalculated scores.

%%
knee_pcastr_vx = varmx_pca(knee_pcastr);

%%
% We can now follow the same convention as before for visualising rotated
% components and scores.

%%
kneeharmfd_vx  = knee_pcastr_vx.harmfd;
kneeharmmat_vx = eval_fd(time, kneeharmfd_vx);
kneevarprop_vx = knee_pcastr_vx.varprop;
kneescores_vx = knee_pcastr_vx.harmscr; 
stdevfPCscores_vx = std(kneescores_vx);

%%
% Below is the visual representation of rotated fPC1 & fPC2 in addition to
% the rotated scores for fPC1 & fPC2 plotted together. 

%%
con3 = stdevfPCscores_vx(1)*1;
fPC1_Pos_vx = kneemeanvec + con3.*(kneeharmmat_vx(:,1));
fPC1_Neg_vx = kneemeanvec - con3.*(kneeharmmat_vx(:,1));

figure(5)
phdl = plot(time, kneemeanvec, 'k-');
set(phdl, 'LineWidth', 1)
hold on
phdl = text(time, fPC1_Pos_vx, '+');
set(phdl, 'LineWidth', 1)
hold on
phdl = text(time, fPC1_Neg_vx, '-');
set(phdl, 'LineWidth', 1)
hold off
box on
set(gca,'FontSize',10);
title('\fontsize{10} fPC1 Varimax Rotated');
xlabel('\fontsize{10} Gait Cycle (%)')
ylabel('\fontsize{10} Knee Angle (°)')

con4 = stdevfPCscores_vx(2)*1;
fPC2_Pos_vx = kneemeanvec + con4.*(kneeharmmat_vx(:,2));
fPC2_Neg_vx = kneemeanvec - con4.*(kneeharmmat_vx(:,2));

figure(6)
phdl = plot(time, kneemeanvec, 'k-');
set(phdl, 'LineWidth', 1)
hold on
phdl = text(time, fPC2_Pos_vx, '+');
set(phdl, 'LineWidth', 1)
hold on
phdl = text(time, fPC2_Neg_vx, '-');
set(phdl, 'LineWidth', 1)
hold off
box on
set(gca,'FontSize',10);
title('\fontsize{10} fPC2 Varimax Rotated');
xlabel('\fontsize{10} Gait Cycle (%)')
ylabel('\fontsize{10} Knee Angle (°)')

figure(7)
phdl = scatter(kneescores_vx(:,1),kneescores_vx(:,2), 'k');
set(phdl, 'LineWidth', 1)
box on
set(gca,'FontSize',10);
title('\fontsize{10} fPC1 Varimax vs fPC2 Varimax');
xlabel('\fontsize{10} fPC1 Scores')
ylabel('\fontsize{10} fPC2 Scores')

%% References
% 
% Craven, P. & Wahba, G. (1979). Smoothing noisy data with spline
% functions: Estimating the correct degree of smoothing by the method of
% generalized crossvalidation. _Numerische Mathematik_ 31, 377–403.
%
% Olshen, R. A., Biden, E. N., Wyatt, M. P. & Sutherland, D. H. (1989).
% Gait analysis and the bootstrap. _The Annals of Statistics_, 1419-1440.
%
% Ramsay, J. O., Hooker, G. & Graves, S. (2009). _Functional data analysis
% with R and MATLAB_, Springer, New York, NY.
%
% Ramsay J. O. & Silverman B. W. _Functional data analysis_, Wiley Online
% Library; 2005.

