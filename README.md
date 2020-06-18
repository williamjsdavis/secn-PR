# secn-PR
Regression of stochastic processes with exponentially correleted noise. `secn-PR` is a MATLAB package designed to estimate a stochastic differential equation from multiple data sources. The regression is flexable enough to account for white or exponentially correlated noise. It employs a block bootstrap algorithm to estimate uncertainties. 

# Theory

To be added

# Using the code
This example follows the script `exampleScript.m` in the `./examples/` directory. Further information on the functions/classes/options can be found in the headers of each file.

### Import some observations
Take some time-series observations that are vectors, and use the `buildObservation()` function to form an `observationData` object. The observations **must have the same time-step**. 

This data set was calculated from a prescribed SDE, so the function estimates can be compared to the true values later.

```MATLAB
%% Make observation object
observationDataAB = buildObservation(trueModel.dt,...
    observationVectorA,observationVectorB);
figure,observationDataAB.plotObservation;
```

<img src="/figures/example_figure01.png" height="400"/>

### Calculate moments
From the observations, conditional moments can be calculated using the `buildMoments()` function. Options for the calculations are in the `MomentOptionsClass` class definition.

```MATLAB
% Setting options (timeShiftSamplePoints,nEvalPoints,evalLims,bandwidth)
momentOptions = MomentOptionsClass(1:60,20,[-1,1],0.1);

% Making moments
momentDataAB = buildMoments(observationDataAB,momentOptions);
figure,momentDataAB.plotMoment1()
```
<img src="/figures/example_figure02.png" height="400"/>

Adjustments to the options can be made and the . The choice of maximum lag should be as small as possible whilst still retaining the shape of the moments.

```MATLAB
% Making adjustments to moment calculations
momentOptions.timeShiftSamplePoints = 1:15;
momentDataAB = buildMoments(observationDataAB,momentOptions);
figure,momentDataAB.plotMoment1()
```
<img src="/figures/example_figure03.png" height="400"/>

### Estimate stochastic process
From the conditional moments, the correlation time and drift and noise functions of the stochastic process can be estimated using the `estimateSPmodel()` function. Options for the calculations are in the `FitOptionsClass` class definition. For example, `fixThetaValue` allows for setting the correlation time to a particular value rather than solving for it.

```MATLAB
% Setting options (no arguments => default)
fitOptions = FitOptionsClass();

% Estimate stochastic process model
SPmodelAB = estimateSPmodel(momentDataAB,fitOptions);

% Review autocorrelation and moments functional fit
figure,SPmodelAB.plotAutocorrelation();
figure,SPmodelAB.plotMoment1();
```

<img src="/figures/example_figure04.png" height="400"/>
<img src="/figures/example_figure05.png" height="400"/>

The drift and noise functions are plotted here with the true input functions.

```MATLAB
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
```

<img src="/figures/example_figure07.png" height="400"/>

### Bootstrap uncertainties
Uncertainties can be estimated through a block bootstrap algorithm, using the `estimateBootstrapModel()` function. Options for the calculations are in the `BootstrapOptionsClass` class definition.

```MATLAB
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
```
<img src="/figures/example_figure08.png" height="400"/>

The drift and noise functions and their uncertainties are plotted here with the true input functions.

This workflow is streamlined in the `fullSPestimate()` function, see `exampleScript2.m` for more details.

# TODOs

- [ ] Implement arbitrary `timeShiftSamplePoints` vector
- [ ] Implement arbitrary spatial evaluation points
- [x] Implement option of forcing a particular correlation time
- [ ] Add compatibility for combining observations of different time-steps
- [x] Make a plotting function for the distribution of correlation times
- [ ] Add more kernel functions to the library

# Changelog

- Version 0.9 - Introduced version I used for my research.

# References

- Lehle, B., & Peinke, J. (2018). Analyzing a stochastic process driven by Ornstein-Uhlenbeck noise. Physical Review E, 97(1), 012113.
- Lamouroux, D., & Lehnertz, K. (2009). Kernel-based regression of drift and diffusion coefficients of stochastic processes. Physics Letters A, 373(39), 3507-3512.
- Kunsch, H. R. (1989). The jackknife and the bootstrap for general stationary observations. The annals of Statistics, 1217-1241.


