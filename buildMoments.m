function outMomentObject = buildMoments(observationData,momentOptions)
%Build conditional moment data
%   William Davis, 17/01/20
%
%   Notes:
%   Takes multiple observation data-sets and creates an observation object
%   for calculating conditional moments. Optimal bandwidth can be estimated
%   through suggestBanwidth() function in MomentOptionsClass function.
%
%   Inputs:
%   - "observationData"         Observation data object, ObservationClass
%   - "momentOptions"           Options, MomentOptionsClass
%       - "timeShiftSamplePoints"   Sampling points in time, vector
%       - "nEvalPoints"             N. points for spatial sampling, double
%       - "evalLims"                Limits of spatial sampling, vector
%       - "kernelType"              Chosen kernel, string
%       - "bandwidth"               Width of kernel, double
%
%   Problems:
%   - Alter to take arbitrary timeShiftSamplePoints vector
%   - Alter to take arbitrary sampled state space
%   - Only epanechnikov kernel implemented
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating moments, kernel method
[evalPoints,Moment1Matrix,Moment2Matrix] = nKBR_moments_fast(...
    observationData.dataCell,...
    momentOptions.timeShiftSamplePoints,momentOptions.nEvalPoints,...
    momentOptions.evalLims,momentOptions.bandwidth);

%% Making object
outMomentObject = MomentClass(Moment1Matrix,Moment2Matrix,...
                observationData,evalPoints,momentOptions);
end
function [Xcentre,M1,M2] = nKBR_moments_fast(X,tauIn,Npoints,xLims,bandwidth)
%Kernel based moments, n-data
%   William Davis, 29/09/19
%
%   Notes:
%   Calculates kernel based moments for a given stochastic time-series.
%   Uses Epanechnikov kernel with built in computational advantages. Uses
%   Nadaraya-Watson estimator. Theory from:
%
%   Lamouroux & Lehnertz, 2009
%   Kernel-based regression of drift and diffusion coefficients of 
%   stochastic processes
%   PHYSICS LETTERS A 373, 3507-3512 (2009)
%
%   Calculates moments from n sources of data.
%   Fast algorithm. Uses running sum of moments, including the naive
%   algorithm for the second moment (error is on the order of 1E-12).
%   Results in a 20-30% speedup.
%   
%
%   Inputs:
%   - "X"                       Observed variables, cell array of data
%   - "tau_in"                  Time-shift indexes
%   - "Npoints"                 Number of evaluation points
%   - "xLims"                   Limits in upper and lower evaluation points
%   - "h"                       Bandwidth
%
%   Problems:
%   - Tidy up
%   - Implement other kernels
%   - Remove yinc{nd} nd indexing?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing
dX = (xLims(2)-xLims(1))/(Npoints-1); % Bins increment
Xcentre = xLims(1):dX:xLims(2); % Grid
heff = bandwidth*sqrt(5); % Effective bandwidth, for setting up bins
eta = floor(heff/dX+0.5); % Bandwidth for bins optimizing

% Epanechnikov kernel
K= @(u) 0*(u.^2>1)+3/4*(1-u.^2).*(u.^2<=1);
Ks = @(u) K(u/sqrt(5))/sqrt(5); % Silverman's definition of the kernel (Silverman, 1986)
Kh = @(u) Ks(u/bandwidth)/bandwidth; % Changing bandwidth

% Sort all data into bins
Bextend = dX*(eta+0.5); % Extend bins
edges = xLims(1)-Bextend:dX:xLims(2)+Bextend; % Edges
ndata = numel(X); % Number of data-sets
Xloc = cell(1,ndata); % Preallocate histogram location data
nXdata = cellfun(@numel,X); % Number of x data
for nd = 1:ndata
    [~,~,Xloc{nd}] = histcounts(X{nd},edges); % Sort
end
Xbinloc = eta+(1:Npoints); % Bin locations
BinBeg = Xbinloc-eta; % Bin beginnings
BinEnd = Xbinloc+eta; % Bin beginnings

% Preallocate
Ntau = numel(tauIn); % Number of time-steps
[M1,M2] = deal(zeros(Ntau,Npoints)); % Moments
[iXkey,Khj,yinc] = deal(cell(1,ndata)); % Preallocate increment data

% Pre calculate increments
inc = cell(Ntau,ndata);
for nd = 1:ndata
    poss_tau_ind = 1:nXdata(nd); % Possible time-shifts
    for tt = 1:Ntau
        tau_c = tauIn(tt); % Chosen shift
        tau_ind = poss_tau_ind(1+tau_c:end); % Chosen indices
        inc{tt,nd} = X{nd}(tau_ind) - X{nd}(tau_ind - tau_c);
    end
end

% Loop over evaluation points
for ii = 1:Npoints
    
    % Start and end bins
    kBinBeg = BinBeg(ii);
    kBinEnd = BinEnd(ii);
    
    % Data and weights
    for nd = 1:ndata
        iXkey{nd} = find((kBinBeg<=Xloc{nd}) & (Xloc{nd}<=kBinEnd)); % Data key
        Khj{nd} = Kh(Xcentre(ii) - X{nd}(iXkey{nd})); % Weights
    end
    
    % For each shift
    for tt = 1:Ntau
        tau_c = tauIn(tt); % Chosen shift
        
        % Initialise moment sums
        sumKhjtt = 0;
        sumKhjtt_ytt = 0;
        sumKhjtt_y2tt = 0;
        
        % Get data
        for nd = 1:ndata            
            XUin = iXkey{nd}; % Unshifted data indices
            XUin = XUin(XUin <= nXdata(nd)-tau_c); % Clip overflow
            yinc{nd} = inc{tt,nd}(XUin); % Increments
            Khjt = Khj{nd}(1:numel(yinc{nd})); % Clipped weight vector
            sumKhjtt = sumKhjtt + sum(Khjt); % Sum of weights
            sumKhjtt_ytt = sumKhjtt_ytt + sum(Khjt.*yinc{nd}); % Weighted increments
            sumKhjtt_y2tt = sumKhjtt_y2tt + sum(Khjt.*yinc{nd}.^2); % Weighted increments^2
        end
        
        % First moment sum
        M1(tt,ii) = sumKhjtt_ytt/sumKhjtt;
        
        % Second moment sum (naive algorithm)
        M2(tt,ii) = (sumKhjtt_y2tt - sumKhjtt_ytt^2/sumKhjtt)/sumKhjtt;
        
    end
end
end