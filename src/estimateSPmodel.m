function outSPmodelObject = estimateSPmodel(momentData,fitOptions)
%Estimate stochastic process model
%   William Davis, 20/01/20
%
%   Notes:
%   Takes moments object and fits for correlation time, theta, and 
%   drift, f(x), and noise, g(x), functions. Theory from:
%
%   Lehle & Peinke, 2018
%   Analyzing a stochastic process driven by Ornstein-Uhlenbeck noise
%   PHYSICAL REVIEW E 97, 012113 (2018)
%
%   Inputs:
%   - "momentData"              Moment data object, MomentClass
%   - "fitOptions"              Options, FitOptionsClass
%       - "thetaConvergence"    Iteration error limit for theta, double
%       - "functionConvergence" Same for drift and noise functions, double
%       - "fixTheta"            Set theta value? logical
%       - "fixThetaValue"       Value of fixed theta, double
%       - "keepObservations"    Keep observationData? logical
%       - "printOutput"         Print to command window? logical
%
%   Problems:
%   - Alter theta_search.m to take arbitrary timeShiftSamplePoints vector.
%   - Alter theta_search.m to give rMatrix for vector above...
%   - Removed possible change in lambda_tau_choice. Investigate further.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get or estimate correlation time (theta)
[thetaStar,thetaProperties] = getTheta(momentData,fitOptions);

%% Estimate drift and noise functional form
% Estimation of lambda vector (linearised)
[lambdaProperties.lambda1Star,lambdaProperties.lambda2Star] = ...
    lambdaSearchLinear(momentData.moment1Matrix,...
    momentData.moment2Matrix,thetaProperties.rMatrix);

% Estimate functions f and g (drift and noise)
lambda1_1 = lambdaProperties.lambda1Star(1,:); % Array of lambda^(1)_1
lambda2_1 = lambdaProperties.lambda2Star(1,:); % Array of lambda^(2)_1
[fEstimate,gEstimate,fInitial,gInitial] = fgIter(lambda1_1,lambda2_1,...
    thetaStar,momentData.evalPoints,fitOptions.functionConvergence);

% Mean fit error
meanAbsoluteError = meanFitError(momentData,thetaProperties,...
    lambdaProperties,fitOptions.printOutput);

%% Making object
outSPmodelObject = SPmodelClass(thetaStar,fEstimate,gEstimate,fInitial,...
    gInitial,meanAbsoluteError,momentData,thetaProperties,...
    lambdaProperties,fitOptions);
end
function meanAbsoluteError = meanFitError(momentData,thetaProperties,...
    lambdaProperties,printOutput)
%Calculating mean fit error
%   William Davis, 22/01/20
%
%   Notes:
%   Calculating mean error between data and fit
%
%   Inputs:
%   - "momentData"              Moment data object, MomentClass
%   - "thetaProperties"         Correlation time properties, structure
%   - "lambdaProperties"        Function fit properties, structure
%   - "printOutput"             Print to command window? logical
%
%   Problems:
%   - Add option to surpress print to terminal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
moment1Difference = momentData.moment1Matrix-...
    thetaProperties.rMatrix*lambdaProperties.lambda1Star;
moment2Difference = momentData.moment2Matrix-...
    thetaProperties.rMatrix*lambdaProperties.lambda2Star;
meanAbsoluteError.moment1 = mean(abs(moment1Difference(:)));
meanAbsoluteError.moment2 = mean(abs(moment2Difference(:)));
meanAbsoluteError.bothMoments = mean([meanAbsoluteError.moment1,...
    meanAbsoluteError.moment2]);
if printOutput
    disp('Absolute fit error on moments')
    disp(['Mean M^(1) fit error: ',...
        sprintf('%0.4e',meanAbsoluteError.moment1)])
    disp(['Mean M^(2) fit error: ',...
        sprintf('%0.4e',meanAbsoluteError.moment2)])
    disp(['Mean fit error: ',...
        sprintf('%0.4e',meanAbsoluteError.bothMoments)])
end
end
function [thetaStar,thetaProperties] = getTheta(momentData,fitOptions)
%Get theta value
%   William Davis, 18/06/20
%
%   Notes:
%   Determine whether to search for theta or use set value.
%
%   Inputs:
%   - "momentData"              Moment data object, MomentClass
%   - "fitOptions"              Options, FitOptionsClass
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum lag is set for fixed theta case
maximumLag = max(momentData.momentOptions.timeShiftSamplePoints);

if fitOptions.fixTheta % Get fixed theta
    [thetaStar,thetaProperties] = theta_chosen(...
        momentData.observationData.dataCell,...
        momentData.observationData.timeStep,maximumLag,fitOptions);
    
else % Search for theta
    maximumTheta = maximumLag*momentData.observationData.timeStep;
    [thetaStar,thetaProperties] = theta_search(...
        momentData.observationData.dataCell,...
        momentData.observationData.timeStep,maximumLag,maximumTheta,...
        fitOptions.thetaConvergence);
end
end
function [thetaStar,thetaProperties] = theta_chosen(X,dt,nuMax,fitOptions)
%
%   William Davis, 18/06/20
%
%   Notes:
%   Get set theta value and relevant properties.
%
%   Inputs:
%   - "X"                       Observed variable, cell array of data
%   - "dt"                      Time-step
%   - "nuMax"                   Index of maximum time-shift
%   - "thetaMax"                Maximum correlation time in search
%   - "fitOptions"              Options, FitOptionsClass
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get set theta value
thetaStar = fitOptions.fixThetaValue;

% Autocorrelation
dA = nAutocorrIncrement(X,nuMax); % Multiple data

% R matrix and lambda vector
rNuMatrix = formRmatrix(dt,nuMax);
rMatrix = rNuMatrix(thetaStar); % Output r(tau,theta) matrix
lambdaStar = rMatrix\dA; 

%% Outputs
thetaProperties.newNuMax = nuMax;
thetaProperties.rMatrix = rMatrix;
thetaProperties.dA = dA;
thetaProperties.lambdaStar = lambdaStar;
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
%   - "dt"                      Time-step
%   - "nuMax"                   Index of maximum time-shift
%   - "thetaMax"                Maximum correlation time in search
%   - "betaConv"                Convergence parameter (not used yet)
%
%   Problems:
%   - See note on nu correction.
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
    acf{nd} = myAutocorr(X{nd},tauMax);
    nX(nd) = numel(X{nd}); % Number of points
    varX(nd) = var(X{nd}); % Variance
    dA{nd} = (acf{nd}(2:end)-acf{nd}(1))*var(X{nd}); % Tidy up array
end

%% Merge data
dAmat = cell2mat(dA); % Make into matrix
dA_full = (dAmat*nX')/sum(nX);
end
function acf = myAutocorr(x,lags)
%Autocorrelation function
%   William Davis, 30/01/20
%
%   Notes:
%   Reduces dependency on Signal Processing Toolbox.
%   
%   Inputs:
%
%   Problems:
%   - No checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demean and calculate
xDemean = x - mean(x);
nFFT = 2^(nextpow2(length(xDemean))+1);
F = fft(xDemean,nFFT);
F = F.*conj(F);
acf = ifft(F);
acf = acf(1:(lags+1));
acf = real(acf);
acf = acf./acf(1);
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
%   Inputs:
%   - "dA"                      Autocorrelation increments
%   - "dt"                      Time-step
%   - "nuMax"                   Index of maximum time-shift
%   - "thetaMax"                Maximum correlation time in search
%   - "betaConv"                Convergence parameter (not used)
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Searching for theta
% Objective function from matrix
rNuMatrix = formRmatrix(dt,nuMax);
lambdaFunc = @(theta_c) rNuMatrix(theta_c)\dA;
funcValue = @(theta_c) sum((dA-rNuMatrix(theta_c)*lambdaFunc(theta_c)).^2);

% Line search using golden section search and parabolic interpolation
thetaStar = fminbnd(funcValue,0,thetaMax);

% Best fit lambda vector
lambdaStar = lambdaFunc(thetaStar);
end
function rNuMatrix = formRmatrix(dt,nuMax)
%Form R matrix
%   William Davis, 18/06/20
%
%   Notes:
%   Recursive formation.
%
%   Inputs:
%   - "dt"                      Time-step
%   - "nuMax"                   Index of maximum time-shift
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basis functions
r1 = @(tau,theta) tau - theta*(1 - exp(-tau/theta));
r2 = @(tau,theta) tau.^2/2 - theta*r1(tau,theta);
r3 = @(tau,theta) tau.^3/6 - theta*r2(tau,theta);
rArray = @(tau,theta) [r1(tau,theta),r2(tau,theta),r3(tau,theta)];

% Functions r matrix (reducing dependencies, new method)
nu = 1:nuMax; % Indexes of time-shifts
tau_nu = nu'*dt; % Array of time-shifts
rNuMatrix = @(theta) rArray(tau_nu,theta);
end
function [fNew,gNew,f0,g0] = fgIter(lambda1_1,lambda2_1,theta,Xcentre,betaConv)
%Iterate functions f and g
%   William Davis, 14/01/19
%
%   Notes:
%   Iterate and converge for functions f and g, based on estimates for
%   lambda vectors. Input x positions are needed since spatial derivatives
%   are calculated in fixed point iterations.
%
%   Inputs:
%   - "lambda1_1"               Array of lambda^(1)_1, function of x
%   - "lambda2_1"               Array of lambda^(2)_1, function of x
%   - "theta"                   Estimate of correlation time
%   - "Xcentre"                 Position of x estimates
%   - "beta_conv"               Convergence parameter
%
%   Problems:
%   - Tidy up
%   - Convergence parameter too small?
%   - Write better convergence algorithm?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings
count_max = 1; % Iteration limits

%% Processing
% Starting values
f0 = lambda1_1; % Drift function
g0 = sqrt(abs(lambda2_1)); % Noise function
fNew = f0;
gNew = g0;

% Plotting for for debugging purposes
%figure(1),hold on
%plot(Xcentre,f_0,'k-','LineWidth',2),xlabel('X'),ylabel('f')
%figure(2),hold on
%plot(Xcentre,g_0,'k-','LineWidth',2),xlabel('X'),ylabel('f')

% Iterate until converged
count = 0; % Initial count
int_err = inf; % Initial error
while (int_err > betaConv) && (count < count_max)
    % Update counter
    count = count + 1; 
    
    % Update functions
    f_old = fNew; 
    g_old = gNew;
    
    % Find new values
    [fNew,gNew] = fixedPointIter(lambda1_1,lambda2_1,f_old,g_old,theta,...
        Xcentre);
    
    % Integration error
    int_err = trapz(Xcentre,(fNew-f_old).^2) + ...
        trapz(Xcentre,(gNew-g_old).^2);
    
    % Plotting for for debugging purposes
    %figure(1),plot(Xcentre,f_new),
    %figure(2),plot(Xcentre,g_new)
    
end

end
function [lambda1,lambda2] = lambdaSearchLinear(M1,M2,r_matrix)
%Lambda search
%   William Davis, 29/05/19
%
%   Notes:
%   Linear searching for lambda, the correlation time of a stochastic
%   process. Lambdas are calculated as a function of x position, although
%   the specific values of the positions are not needed. Only the size of
%   the x position vector is needed. Linearised, faster.
%
%   Inputs:
%   - "M1"                      First moments
%   - "M2"                      Second moments
%   - "r_matrix"                Function r(tau,theta) matrix
%
%   Problems:
%   - Use original tau_max or second tau_max?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over requested x points
% nx_points = size(M1,2); % Number of x points

% lambda1 = zeros(3,nx_points); % Preallocating lambda^(1) array
% lambda2 = zeros(3,nx_points); % Preallocating lambda^(2) array
% for ii = 1:nx_points
%     % Minimising normal matrix
%     lambda1(:,ii) = r_matrix\M1(:,ii); 
%     lambda2(:,ii) = r_matrix\M2(:,ii);
% end

%% Solve linear equations
% inv(r_matrix'*r_matrix)*r_matrix'*M1
lambda1 = r_matrix\M1;
lambda2 = r_matrix\M2;

end
function [fNew,gNew] = fixedPointIter(lambda1_1,lambda2_1,f,g,theta,...
    Xcentre)
%Fixed point iteration
%   William Davis, 14/01/19
%
%   Notes:
%   Fixed point iteration to update functions f and g. Assumes a small
%   theta.
%   
%
%   Inputs:
%   - "lambda1_1"               Array of lambda^(1)_1, function of x
%   - "lambda2_1"               Array of lambda^(2)_1, function of x
%   - "f"                       Array, function f estimate
%   - "g"                       Array, function g estimate
%   - "theta"                   Estimate of correlation time
%   - "Xcentre"                 Position of x estimates
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing
% Spatial f and g derivates
%f_div = gradient(f,Xcentre);
%g_div = gradient(g,Xcentre);
f_div = fdiffNU(Xcentre,f);
g_div = fdiffNU(Xcentre,g);

% Updated functions
fNew = lambda1_1 - 0.5*g.*g_div - ...
    0.5*theta*(f_div.*g.*g_div - f.*g_div.^2);
g_arg = lambda2_1 - theta*(f_div.*g.^2 - f.*g.*g_div);
gNew = sqrt(abs(g_arg));

end
function dFdx = fdiffNU(x,F)
%Irregular grid derivative
%   William Davis, 01/05/19
%
%   Notes:
%   Based off code from Bruce.
%   Polynomial interpolation and finite difference differentiation for 
%   irregular grid spacing. Notes at:
%
%   http://cfd.mace.manchester.ac.uk/twiki/pub/Main/...
%                   TimCraftNotes_All_Access/cfd1-findiffs.pdf
%
%   Second order accuracy. Three point centered stencil for body and three
%   point offset stencil for boundary points. Vectorised for efficiency.
%   Checked by hand and numerically.
%
%   Input parameters:
%   - "x"                           Abscissa
%   - "F"                           Function values
%
%   Problems:
%   - Write second order code?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing
n = length(x); % Length
m = n-2; % End point for body

h0 = x(2:n-1) - x(1:n-2); % Backwards steps
h1 = x(3:n) - x(2:n-1); % Forward steps
den = h0.*h1.*(h0+h1); % Denominator
h0_2 = h0.^2; % Backwards steps squared
h1_2 = h1.^2; % Forward steps squared

% Boundaries
dFdx(n) = ((h0_2(m)+2*h0(m)*h1(m))*F(n) - ...
    (h0(m)+h1(m))^2*F(n-1) + h1_2(m)*F(n-2))/den(m);
dFdx(1) = (-(h1_2(1) + 2*h0(1)*h1(1))* F(1) + ...
    (h0(1)+h1(1))^2*F(2) - h0_2(1)*F(3))/den(1);

% Main Body
dFdx(2:n-1)= (-h1_2.*F(1:n-2) + (h1_2-h0_2).*F(2:n-1) + h0_2.*F(3:n))./den;
end