%% Example Script 3
% William Davis, 10/02/22
% Illustrating a third way to use the secPR package
clearvars,clc,close all
addpath('../src/')
%% Import data and true model
[observationVectorA,observationVectorB,trueModel] = importExampleData();

%% Make observation object
observationDataAB = buildObservation(trueModel.dt,...
    observationVectorA,observationVectorB);
figure,observationDataAB.plotObservation();

%% Make moment matrices
% Setting options (timeShiftSamplePoints,nEvalPoints,evalLims,bandwidth)
timeShiftSamplePoints = 40:60;
momentOptions = MomentOptionsClass(timeShiftSamplePoints,20,[-1,1],0.3);
momentDataAB = buildMoments(observationDataAB,momentOptions);
figure,momentDataAB.plotMoment1()

%% Fitting to functions
% Setting options
fitOptions = FitOptionsClass();
fitOptions.fixTheta = true;
fitOptions.fixThetaValue = 0;

% Estimate stochastic process model
SPmodelAB = estimateSPmodel(momentDataAB,fitOptions);

% Review autocorrelation and moments functional fit
figure,SPmodelAB.plotAutocorrelation();
xlim([0,inf])
figure,SPmodelAB.plotMoment1();
xlim([0,inf])
figure,SPmodelAB.plotMoment2();
xlim([0,inf])

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

% %% Bootstrapping uncertainties
% % Setting options (nBootstrapSamples)
% bootstrapOptions = BootstrapOptionsClass(20);
% 
% % Estimate uncertainties
% SPbootstrapAB = estimateBootstrapModel(SPmodelAB,bootstrapOptions);
% 
% % View estimated uncertainties
% figure
% subplot(1,2,1)
% h = SPbootstrapAB.plotDrift();
% h.model = plot(SPmodelAB.momentData.evalPoints,...
%     trueModel.drift(SPmodelAB.momentData.evalPoints),'r-','LineWidth',2);
% legend([h.scatter,h.fill,h.model],'Best estimate',...
%     '95% confidence interval','True model')
% subplot(1,2,2)
% h = SPbootstrapAB.plotNoise();
% h.model = plot(SPmodelAB.momentData.evalPoints,...
%     trueModel.noise(SPmodelAB.momentData.evalPoints),'r-','LineWidth',2);
% legend([h.scatter,h.fill,h.model],'Best estimate',...
%     '95% confidence interval','True model')
