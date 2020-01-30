classdef MomentOptionsClass
%Moment option class
%   William Davis, 27/01/20
%
%   Notes:
%   Defining class for calculation options for conditional moments from 
%   observation data. The arguements of the constructor of the class are
%   the necessary properties that have no default value.
%
%   Inputs:
%   - "timeShiftSamplePoints"   Sampling points in time, vector
%   - "nEvalPoints"             N. points for spatial sampling, double
%   - "evalLims"                Limits of spatial sampling, vector
%   - "kernelType"              Chosen kernel (default), string
%   - "bandwidth"               Width of kernel, double
%
%   Problems:
%   - Only epanechnikov kernel implemented
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        timeShiftSamplePoints(1,:) double {mustBeNumeric,mustBeFinite,...
            mustBeReal,mustBeInteger}
        nEvalPoints(1,1) double {mustBeNumeric,mustBeFinite,...
            mustBeNonempty,mustBeReal,mustBeInteger}
        evalLims(1,2) double {mustBeNumeric,mustBeFinite,mustBeNonempty,...
            mustBeNonempty,mustBeReal}
        kernelType(1,:) char {mustBeMember(kernelType,...
            {'epanechnikov'})} = 'epanechnikov'
        bandwidth(1,1) double {mustBeNumeric,mustBeFinite,...
            mustBeNonempty,mustBeReal}
    end
    
    methods
        function obj = MomentOptionsClass(timeShiftSamplePoints,...
                nEvalPoints,evalLims,bandwidth)
            %%%Constructor method, initialise necessary properties
            % Validate emptyness of vector
            if isempty(timeShiftSamplePoints)
                error(['Error setting property ''nEvalPoints'' of ',...
                    'class ''MomentOptionsClass'':\n%s'],...
                    'Size of value must be scalar.')
            else
                obj.timeShiftSamplePoints = timeShiftSamplePoints;
            end
            obj.nEvalPoints = nEvalPoints;
            obj.evalLims = evalLims;
            obj.bandwidth = bandwidth;
        end
        
        function optimumBandwidth = suggestBanwidth(obj,observationData)
            %%%Normal reference bandwidth selector
            % Theory from:
            % Fan, Jianqing. 
            % Local polynomial modelling and its applications: 
            % Monographs on statistics and applied probability 66. 
            % Routledge, 2018.
            % Section 2.7.1
            
            % Multiplicative constant
            switch obj.kernelType
                case 'epanechnikov'
                    scaleConstant = 2.34;
                otherwise
                    error('Undefined kernel type')
            end
            
            % Correlation-time estimate
            maximumLag = max(obj.timeShiftSamplePoints);
            maximumTheta = maximumLag*observationData.timeStep;
            [correlationTime,~] = theta_search(observationData.dataCell,...
                observationData.timeStep,maximumLag,maximumTheta,1E-2);
            
            % Number of time-steps per correlation-time
            observationsPerCorrelation = ...
                ceil(correlationTime/observationData.timeStep);
            
            % Scaled number of observations (correction for correlation)
            scaledN = ...
                observationData.dataTotalLength/observationsPerCorrelation;
            
            % Sum standard deviation function
            weightedStandardDeviation = ...
                @(s,n) sqrt((s.^2*(n-1)')/sum(n-1)); 
            
            % Sample standard deviation
            sampleSdev = cellfun(@(x) std(x,0),observationData.dataCell);
            scaledS = weightedStandardDeviation(...
                sampleSdev,observationData.dataLength);
            
            % Normal reference bandwidth selector
            optimumBandwidth = scaleConstant*scaledS*scaledN^(-1/5); 
            
        end
    end
end
function [thetaStarNew,thetaProperties] = theta_search(X,dt,nuMax,...
    thetaMax,betaConv)
%Theta search
%   William Davis, 13/01/19
%
%   Notes:
%   Non-linear searching for theta, the correlation time of a stochastic
%   process. Fitting the autocorrelation increments to basis functions.
%   Two searches are made.
%   
%   Tau from 1 to tau_max will be used in the first search, and once theta
%   is estimated, theta becomes the tau_max for the second search. The
%   second theta is used, due to scaling arguments.
%
%   Altered to incorporate multiple data sources.
%   Grouped outputs into "thetaProperties".
%
%   Inputs:
%   - "X"                       Observed variable, cell array of data
%   - "t"                       Time
%   - "dt"                      Time-step
%   - "nu_max"                  Index of maximum time-shift
%   - "theta_max"               Maximum correlation time in search
%   - "beta_conv"               Convergence parameter
%
%   Problems:
%   - See note on nu correction.
%   - Tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing
% Autocorrelation
dA = nAutocorrIncrement(X,nuMax); % Multiple data

% First search
[thetaStarInitial,rNuMatrix,~] = thetaBasisFunctionFit(dA,dt,nuMax,...
    thetaMax,betaConv);

% New maximum tau
newNuMax = ceil(sqrt(thetaStarInitial)/dt);

if newNuMax > nuMax
    % I think nu needs to be as small as possible.
    % The case of taking a new nu max might only be appropriate for time-
    % series data that has bigger time steps. For our purposes, the time-
    % step is quite small, relative to the characteristic time-scale of the
    % process. If we do use the new nu maximum, what effectively happens is
    % the estimates from lambda get dominated by the first r function (the
    % linear term). This overpowers the other two terms, and means that the
    % algorithm cannot sense the correlation. I.e. the theta_bf_fit
    % function for the new nu maximum has a minimum very close to 0, as
    % opposed to the original theta_bf_fit run that has a finite time
    % minimum. This may need further investigation.
    
    newNuMax = nuMax;
    
    dA = nAutocorrIncrement(X,newNuMax);
    
end

% Second search
newThetaMax = newNuMax*dt; % New theta maximum
[thetaStarNew,~,lambdaStar] = thetaBasisFunctionFit(dA(1:newNuMax),dt,...
    newNuMax,newThetaMax,betaConv);

% R matrix
rMatrix = rNuMatrix(thetaStarNew); % Output r(tau,theta) matrix

%% Outputs
thetaProperties.newNuMax = newNuMax;
thetaProperties.rMatrix = rMatrix;
thetaProperties.dA = dA;
thetaProperties.lambdaStar = lambdaStar;
end
function [dA_full] = nAutocorrIncrement(X,tauMax)
%Autocorrelation increments, n-data
%   William Davis, 09/07/19
%
%   Notes:
%   Calculates the (non-normalised) autocorrelation increments for a given
%   maximum time-shift. Calculates increments from n sources of data.
%
%   Inputs:
%   - "X"                       Observed variable, cell array of data
%   - "nu_max"                  Index of maximum time-shift
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing
ndata = numel(X); % Number of data-sets

[acf,dA] = deal(cell(1,ndata)); % Preallocate
[nX,varX] = deal(zeros(1,ndata));
for nd = 1:ndata
    % Built in function (faster)
    [acf{nd},~,~] = autocorr(X{nd},'NumLags',tauMax); 
    nX(nd) = numel(X{nd}); % Number of points
    varX(nd) = var(X{nd}); % Variance
    dA{nd} = (acf{nd}(2:end)-acf{nd}(1))*var(X{nd}); % Tidy up array
end

%% Merge data
dAmat = cell2mat(dA); % Make into matrix
dA_full = (dAmat*nX')/sum(nX);

% Old code
%[acf,~,~] = autocorr(X,'NumLags',tau_max); % Built in function (faster)
%dA = (acf(2:end)-acf(1))*var(X); % Tidy up array
end
function [thetaStar,rNuMatrix,lambdaStar] = thetaBasisFunctionFit(dA,dt,...
    nuMax,thetaMax,betaConv)
%Theta basis function fit
%   William Davis, 12/01/19
%
%   Notes:
%   Non-linear searching for theta, the correlation time of a stochastic
%   process. Fits autocorrelation to basis functions to find optimal theta.
%   Assumes all tau from 1 to tau_max will be used in search.
%   
%
%   Inputs:
%   - "dA"                      Autocorrelation increments
%   - "dt"                      Time-step
%   - "nu_max"                  Index of maximum time-shift
%   - "theta_max"               Maximum correlation time in search
%   - "beta_conv"               Convergence parameter (not used yet)
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing
nu = 1:nuMax; % Indexes of time-shifts
tau_nu = nu'*dt; % Array of time-shifts

% Basis functions
r1 = @(tau,theta) tau - theta*(1 - exp(-tau/theta));
r2 = @(tau,theta) tau.^2/2 - theta*r1(tau,theta);
r3 = @(tau,theta) tau.^3/6 - theta*r2(tau,theta);
rArray = @(tau,theta) [r1(tau,theta),r2(tau,theta),r3(tau,theta)];

% Functions r matrix (reducing dependencies, new method)
rNuMatrix = @(theta) rArray(tau_nu,theta);

%% Searching for theta
% Objective function
lambdaFunc = @(theta_c) rNuMatrix(theta_c)\dA;
funcValue = @(theta_c) sum((dA-rNuMatrix(theta_c)*lambdaFunc(theta_c)).^2);

% Line search using golden section search and parabolic interpolation
thetaStar = fminbnd(funcValue,0,thetaMax);

% Best fit lambda vector
lambdaStar = lambdaFunc(thetaStar);
end