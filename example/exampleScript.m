%% Example Script
% William Davis, 17/01/20
% Illustrating how to use the secPR package
clearvars,clc,close all

%% Import data and true model
[observationVectorA,observationVectorB,trueModel] = importExampleData();

%% Make observation object
observationDataAB = buildObservation(trueModel.dt,...
    observationVectorA,observationVectorB);
figure,observationDataAB.plotObservation();

%% Make moment matrices
% Setting options (timeShiftSamplePoints,nEvalPoints,evalLims,bandwidth)
momentOptions = MomentOptionsClass(1:60,20,[-1,1],0.1);

% Making moments
momentDataAB = buildMoments(observationDataAB,momentOptions);
figure,momentDataAB.plotMoment1()

% Making adjustments to moment calculations
momentOptions.timeShiftSamplePoints = 1:15;
momentDataAB = buildMoments(observationDataAB,momentOptions);
figure,momentDataAB.plotMoment1()

%% Fitting to functions
% Setting options (no arguments => default)
fitOptions = FitOptionsClass();

% Estimate stochastic process model
SPmodelAB = estimateSPmodel(momentDataAB,fitOptions);

% Review autocorrelation and moments functional fit
figure,SPmodelAB.plotAutocorrelation();
figure,SPmodelAB.plotMoment1();
figure,SPmodelAB.plotMoment2();

% View estimated functions
figure
subplot(1,2,1)
h = SPmodelAB.plotDrift();
h.model = plot(SPmodelAB.momentData.evalPoints,...
    trueModel.drift(SPmodelAB.momentData.evalPoints),'r-','LineWidth',2);
legend([h.scatter,h.model],{'Data','True model'},'Location','NorthEast')
subplot(1,2,2)
h = SPmodelAB.plotNoise();
h.model = plot(SPmodelAB.momentData.evalPoints,...
    trueModel.noise(SPmodelAB.momentData.evalPoints),'r-','LineWidth',2);
legend([h.scatter,h.model],{'Data','True model'},'Location','NorthWest')

%% Bootstrapping uncertainties
% Setting options (nBootstrapSamples)
bootstrapOptions = BootstrapOptionsClass(20);

% Estimate uncertainties
SPbootstrapAB = estimateBootstrapModel(SPmodelAB,bootstrapOptions);

% View estimated uncertainties
figure
subplot(1,2,1)
h = SPbootstrapAB.plotDrift();
h.model = plot(SPmodelAB.momentData.evalPoints,...
    trueModel.drift(SPmodelAB.momentData.evalPoints),'r-','LineWidth',2);
legend([h.scatter,h.fill,h.model],'Best estimate',...
    '95% confidence interval','True model')
subplot(1,2,2)
h = SPbootstrapAB.plotNoise();
h.model = plot(SPmodelAB.momentData.evalPoints,...
    trueModel.noise(SPmodelAB.momentData.evalPoints),'r-','LineWidth',2);
legend([h.scatter,h.fill,h.model],'Best estimate',...
    '95% confidence interval','True model')
