function [observationA,observationB,trueModel] = importExampleData
%Build conditional moment data
%   William Davis, 27/01/20
%
%   Notes:
%   Import example data for illustration purposes.
%
%   Inputs:
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Observations
observationA = importdata('./exampleDataA.dat'); % Observation A
observationB = importdata('./exampleDataB.dat'); % Observation B

%% True parameters used in calculations
trueModel.dt = 0.005; % Time-step
trueModel.theta = 0.01; % Correlation time
trueModel.drift = @(x) -2.5*x; % Drift function
trueModel.noise = @(x) 1+x.^2/4.*(x>0); % Noise function
end