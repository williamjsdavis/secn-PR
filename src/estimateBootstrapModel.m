function outSPbootstrapObject = estimateBootstrapModel(SPmodel,bootstrapOptions)
%Estimating bootstrapped stochastic process model 
%   William Davis, 21/01/20
%
%   Notes:
%   Takes estimated stochastic process model and uses a block bootstrap
%   algorithm to estimate uncertainties on correlation time and drift and
%   noise function parameters.
%
%   Inputs:
%   - "SPmodel"                 Stochastic process model object, 
%                                   SPmodelClass
%   - "bootstrapOptions"        Options, BootstrapOptionsClass
%       - "blockLength"         Length of blocks, double
%       - "nBootstrapSamples"   Number of bootstrap samples, double
%       - "parallelProcess"     Switch for parallel processing, logical
%
%   Problems:
%   - Recover fit errors?
%   - Add option to surpress print to terminal?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Printing to terminal
disp('Starting correlated process regression and bootstrap iterations')
disp(['Number of data sources: ',num2str(...
    SPmodel.momentData.observationData.dataNumber)])
disp(['Total number of data points: ',num2str(...
    SPmodel.momentData.observationData.dataTotalLength)])
disp(['Parallel processing: ',mat2str(bootstrapOptions.parallelProcess)])
%disp(['Kernel bandwidth: ',response.h])
disp(['Number of bootstrap iterations: ',...
    num2str(bootstrapOptions.nBootstrapSamples)])
disp(['Bootstrap block length: ',num2str(bootstrapOptions.blockLength)])
disp('Starting calculations')

%% Bootstrapping uncertainties
% Get settings
handleSPmodel = SPmodel;
handleSPmodel.fitOptions.keepObservations = false;
handleSPmodel.fitOptions.printOutput = false;
handleBootstrapOptions = BootstrapOptionsClass(0);

% Input function
SPhandle = @(x) fullSPestimate(x,handleSPmodel.momentData.momentOptions,...
    handleSPmodel.fitOptions,handleBootstrapOptions);

% Full data estimate
originalEstimate = {SPmodel.correlationEstimate,SPmodel.driftEstimate,...
    SPmodel.noiseEstimate}; 

% Bootstrapping
[unbiasedBootstrapDistributions,meanAbsoluteError] = blockBootstrap(...
    SPmodel.momentData.observationData,...
    SPhandle,originalEstimate,bootstrapOptions.blockLength,...
    bootstrapOptions.nBootstrapSamples,...
    bootstrapOptions.parallelProcess);

%% Bootstrap results statistics
[standardErrors,percentiles95] = calculateBootstrapStatistics(...
    unbiasedBootstrapDistributions);

% Finished
disp('Finished calculations')

%% Making object
outSPbootstrapObject = SPbootstrapClass(unbiasedBootstrapDistributions,...
    standardErrors,percentiles95,meanAbsoluteError,SPmodel,...
    bootstrapOptions);
end
function [standardErrors,percentiles95] = ...
    calculateBootstrapStatistics(distributions)
%Getting default bootstrapping settings
%   William Davis, 21/01/20
%
%   Notes:
%   Calculating statistics on bootstrap estimates. Standard error, and 95%
%   percentiles are reported.
%
%   Inputs:
%   - "distributions"           Distributions of correlation time and 
%                               functions, structure
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Standard errors
standardErrors.correlationEstimate = ...
    std(distributions.correlationEstimate);
standardErrors.driftEstimate = std(distributions.driftEstimate);
standardErrors.noiseEstimate = std(distributions.noiseEstimate);

%% Percentiles (95%)
percentileLevel = 95;
percentileLowerBound = (100-percentileLevel)/2;
percentileUpperBound = (100+percentileLevel)/2;

percentiles95.correlationEstimate(2) = ...
    prctile(distributions.correlationEstimate,percentileLowerBound);
percentiles95.correlationEstimate(1) = ...
    prctile(distributions.correlationEstimate,percentileUpperBound);
percentiles95.driftEstimate(2,:) = ...
    prctile(distributions.driftEstimate,percentileLowerBound);
percentiles95.driftEstimate(1,:) = ...
prctile(distributions.driftEstimate,percentileUpperBound);
percentiles95.noiseEstimate(2,:) = ...
prctile(distributions.noiseEstimate,percentileLowerBound);
percentiles95.noiseEstimate(1,:) = ...
    prctile(distributions.noiseEstimate,percentileUpperBound);
end
function [unbiasedBootstrapDistributions,meanAbsoluteError] = ...
    blockBootstrap(observationData,fullFunctionHandle,originalEstimate,...
    blockLength,nBootstrapSamples,parallelProcess)
%Block bootstrap uncertainties
%   William Davis, 15/07/19
%
%   Notes:
%   Block bootstrap algorithm. Based on moving block bootstrap method of 
%   Kunsch (1989):
%
%   The Jackknife and the Bootstrap for General Stationary Observations,
%   The Annals of Statistics, Vol. 17, No. 3 (Sep., 1989), pp. 1217-1241
%
%   Inputs:
%   - "observationData"         Observation data object, ObservationClass
%   - "fullFunctionHandle"      Function of data to result, function_handle
%   - "originalEstimate"        Original data estimate, cell
%   - "blockLength"             Length of blocks, double
%   - "nBootstrapSamples"       N. of bootstrap iterations, double
%   - "nEvalPoints"             N. points for spatial sampling, double
%   - "parallelProcess"         Use parallel processing? logical
%
%   Problems:
%   - Check biases
%   - Explore spmd option?
%   - Add option to surpress print to terminal?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bootstrap sampling
% Make method
drawFunction = makeBootstrapMethod(observationData.dataNumber,...
    observationData.dataLength,blockLength);

% Draw sample and pass to function (can be done parallel)
bootstrapEstimate = cell(nBootstrapSamples,1); % Preallocate
meanAbsoluteError = zeros(nBootstrapSamples,1);
if parallelProcess
    parfor n = 1:nBootstrapSamples
        disp(['Draw: ',num2str(n)])
        [bootstrapEstimate{n},meanAbsoluteError(n)] = ...
            parallelBootstrap(observationData,...
            drawFunction,fullFunctionHandle);
    end
else
    for n = 1:nBootstrapSamples
        disp(['Draw: ',num2str(n)])
        [bootstrapEstimate{n},meanAbsoluteError(n)] = ...
            parallelBootstrap(observationData,...
            drawFunction,fullFunctionHandle);
    end
end

%% Make distributions
bootstrapDistributions = makeDistributions(bootstrapEstimate,...
    nBootstrapSamples);

%% Correct for bias
unbiasedBootstrapDistributions = correctBootstrapBias(...
    bootstrapDistributions,originalEstimate,blockLength);
end
function drawFunction = makeBootstrapMethod(dataNumber,dataLength,...
    blockLength)
%Making bootstrap sampling method
%   William Davis, 23/01/20
%
%   Notes:
%   Methodology for bootstrap sampling. Creates a function handle that
%   creates random bootstrap samples for the input data observations.
%
%   Inputs:
%   - "dataNumber"              Number of observations, double
%   - "dataLength"              Observation lengths, vector
%   - "blockLength"             Length of blocks, double
%
%   Problems:
%   -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sampling functions
rdraw = @(N) randi(N,1,ceil(N/blockLength)); % Random draw
drawin = @(x) x:(x+blockLength-1); % Select indices of blocks
moddraw = @(x,N) mod(x-1,N) + 1; % Modulo
drawN = @(x,N) moddraw(x(1:N),N); % Cut after N samples
allin = @(x) cell2mat(arrayfun(drawin,x,'UniformOutput',false)); % Concat.
bsdrawi = @(N) drawN(allin(rdraw(N)),N); % Draw pseudo sample indices

%% Making function handle
drawFunction = cell(1,dataNumber);
for nd = 1:dataNumber
    drawFunction{nd} = @(x) bsdrawi(dataLength(nd));
end
end
function [bootstrapEstimate,meanAbsoluteError] = parallelBootstrap(...
    observationData,drawFunction,fullFunctionHandle)
%Parallel bootstrap function
%   William Davis, 17/07/19
%
%   Notes:
%   Parallel passing of data for bootstrap sampling and function
%   evaluation. Written in this way to avoid MATLAB thinking there are
%   broadcast variables.
%
%   Inputs:
%   - "observationData"         Observation data object, ObservationClass
%   - "drawFunction"            Bootstrap draw function, function_handle
%   - "fullFunctionHandle"      Function of data to result, function_handle
%
%   Problems:
%   -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing
% Make bootstrap sample
bootstrapData = makeBootstrapSample(observationData,drawFunction);

% Function of bootstrap data
bootstrapAllResults = fullFunctionHandle(bootstrapData);
bootstrapEstimate = {bootstrapAllResults.correlationEstimate,...
    bootstrapAllResults.driftEstimate,bootstrapAllResults.noiseEstimate};
meanAbsoluteError = bootstrapAllResults.meanAbsoluteError.bothMoments;
end
function bootstrapData = makeBootstrapSample(observationData,drawFunction)
%Making bootstrap sample data
%   William Davis, 22/01/20
%
%   Notes:
%   Takes original observations and makes a block bootstrap sample to be
%   used as surrogate data.
%
%   Inputs:
%   - "observationData"         Observation data object, ObservationClass
%   - "drawFunction"            Bootstrap draw function, function_handle
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing
cellData = arrayfun(@(n) observationData.dataCell{n}(drawFunction{n}()),...
    1:observationData.dataNumber,'UniformOutput',false);

% Making object
bootstrapData = ObservationClass(observationData.timeStep,cellData);
end
function distributions = makeDistributions(bootstrapEstimate,...
    nBootstrapSamples)
%Extract distributions from bootstrap results
%   William Davis, 23/01/20
%
%   Notes:
%   Extract results for correlation time and drift and noise functions for
%   each bootstrap result, and make into a distribution. Distributions are
%   either vectors or matricies.
%
%   Inputs:
%   - "bootstrapEstimate"       Cell array of bootstrap results, cell
%   - "bootstrapSamples"        Number of bootstrap iterations, double
%   - "nEvalPoints"             N. points for spatial sampling, double
%
%   Problems:
%   -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing correlation time
distributions.correlationEstimate = ...
    arrayfun(@(n) bootstrapEstimate{n}{1},1:nBootstrapSamples);

%% Processing functions
distributions.driftEstimate = ...
    cell2mat(cellfun(@(x) x{2},bootstrapEstimate,'UniformOutput',false));
distributions.noiseEstimate = ...
    cell2mat(cellfun(@(x) x{3},bootstrapEstimate,'UniformOutput',false));
end
function unbiasedDistributions = correctBootstrapBias(...
    bootstrapDistributions,originalEstimate,blockLength)
%Correct for block size bias 
%   William Davis, 16/07/19
%
%   Notes:
%   Scales result from block bootstrap sample; adjusting for bias. 
%
%   Inputs:
%   - "bootstrapDistributions"  Distribution of bootstrap results, matrix
%   - "originalEstimate"        Original data estimate, cell
%   - "blockLength"             Length of blocks (not used yet), double
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing correlation time
correlationEstimateBias = ...
    mean(bootstrapDistributions.correlationEstimate) - originalEstimate{1};
unbiasedDistributions.correlationEstimate = ...
    bootstrapDistributions.correlationEstimate - correlationEstimateBias;

%% Processing functions
driftEstimateBias = mean(bootstrapDistributions.driftEstimate) - ...
    originalEstimate{2};
unbiasedDistributions.driftEstimate = ...
    bootstrapDistributions.driftEstimate - driftEstimateBias;

noiseEstimateBias = mean(bootstrapDistributions.noiseEstimate) - ...
    originalEstimate{3};
unbiasedDistributions.noiseEstimate = ...
    bootstrapDistributions.noiseEstimate - noiseEstimateBias;
end