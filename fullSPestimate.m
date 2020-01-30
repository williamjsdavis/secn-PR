function outModel = fullSPestimate(observationData,momentOptions,...
    fitOptions,bootstrapOptions)
%Full function to estimate stochastic process model
%   William Davis, 22/01/20
%
%   Notes:
%   Full function for estimating stochastic process model. Chains together
%   making moments and estimating correlation time and drift and noise
%   functions. Also optionally calculates block bootstrap uncertainties.
%
%   Inputs:
%   - "observationData"         Observation data object, ObservationClass
%   - "estimationOptions"       Options for calculations, structure
%       - "timeShiftSamplePoints"   Sampling points in time, vector
%       - "nEvalPoints"             N. points for spatial sampling, double
%       - "evalLims"                Limits of spatial sampling, vector
%       - "kernelType"              Chosen kernel, string
%       - "bandwidth"               Width of kernel, double
%       - "thetaConvergence"        Iteration error limit for theta, double
%       - "functionConvergence"     Same for drift&noise functions, double
%       - "blockLength"             Length of blocks, double
%       - "bootstrapSamples"        Number of bootstrap samples, double
%       - "parallelProcess"         Switch for parallel processing, logical
%
%
%   - "momentOptions"           Options, MomentOptionsClass
%       - "timeShiftSamplePoints"   Sampling points in time, vector
%       - "nEvalPoints"             N. points for spatial sampling, double
%       - "evalLims"                Limits of spatial sampling, vector
%       - "kernelType"              Chosen kernel, string
%       - "bandwidth"               Width of kernel, double
%   - "fitOptions"              Options, FitOptionsClass
%       - "thetaConvergence"        Iteration error limit for theta, double
%       - "functionConvergence"     Same for drift and noise funs., double
%       - "keepObservations"        Keep observationData? logical
%       - "printOutput"             Print to command window? logical
%   - "bootstrapOptions"        Options, BootstrapOptionsClass
%       - "blockLength"             Length of blocks, double
%       - "nBootstrapSamples"       Number of bootstrap samples, double
%       - "parallelProcess"         Switch for parallel processing, logical
%
%   Problems:
%   - Build in default settings for moments?
%   - Investigate possible recursion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate stochastic process
% Build moments
momentData = buildMoments(observationData,momentOptions);

% Estimate stochastic process model
outModel = estimateSPmodel(momentData,fitOptions);

%% Bootstrap estimate uncertainties if required
if bootstrapOptions.nBootstrapSamples ~= 0
    outModel = estimateBootstrapModel(outModel,bootstrapOptions);
end
end