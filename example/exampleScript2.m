%% Example Script 2
% William Davis, 24/01/20
% Illustrating another way to use the secPR package
clearvars,clc,close all

%% Import data and true model
[observationVectorA,observationVectorB,trueModel] = importExampleData();

%% Make observation object
observationDataAB = buildObservation(trueModel.dt,...
    observationVectorA,observationVectorB);

%% Full estimation function
% Options
momentOptions = MomentOptionsClass(1:15,25,[-1,1],0.1);
fitOptions = FitOptionsClass();
fitOptions.fixTheta = true;
fitOptions.fixThetaValue = 0.01;
bootstrapOptions = BootstrapOptionsClass(100);

% Run full estimation
SPbootstrapAB = fullSPestimate(observationDataAB,momentOptions,...
    fitOptions,bootstrapOptions);

%% Plot results

% Functions
figure
subplot(1,2,1)
h = SPbootstrapAB.plotDrift();
h.model = plot(SPbootstrapAB.SPmodel.momentData.evalPoints,...
    trueModel.drift(SPbootstrapAB.SPmodel.momentData.evalPoints),...
    'r-','LineWidth',2);
legend([h.scatter,h.fill,h.model],'Best estimate',...
    '95% confidence interval','True model')
subplot(1,2,2)
h = SPbootstrapAB.plotNoise();
h.model = plot(SPbootstrapAB.SPmodel.momentData.evalPoints,...
    trueModel.noise(SPbootstrapAB.SPmodel.momentData.evalPoints),...
    'r-','LineWidth',2);
legend([h.scatter,h.fill,h.model],'Best estimate',...
    '95% confidence interval','True model')
