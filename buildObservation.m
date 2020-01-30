function outObservationObject = buildObservation(timeStep,varargin)
%Build observation data object
%   William Davis, 17/01/20
%
%   Notes:
%   Takes multiple observation data-sets and creates an observation object
%   for calculating conditional moments. Currently, must have the same
%   time-step for all data.
%
%   Inputs:
%   - "timeStep"                Time-step of data, double
%   - "X1"                      Observation 1, vector
%   - "X2"                      Observation 2, vector
%   - ...                       ...
%   - "Xn"                      Observation n, vector
%
%   Problems:
%   - Implement differing time-steps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-processing
N_observations = numel(varargin);

% Form into cell array
CellData = cell(1,N_observations);
for n = 1:N_observations
    CellData{n} = varargin{n}(:);
end

%% Making object
outObservationObject = ObservationClass(timeStep,CellData);
end