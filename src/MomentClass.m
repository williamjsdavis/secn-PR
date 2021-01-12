classdef MomentClass
%Moment data class
%   William Davis, 17/01/20
%
%   Notes:
%   Defining class for conditional moments from observation data.
%
%   Inputs:
%   - "moment1Matrix"           First moments, matrix
%   - "moment2Matrix"           First moments, matrix
%   - "observationData"         Observation data object, ObservationClass
%   - "evalPoints"              Spatial sampling points, vector
%   - "momentOptions"           Options, MomentOptionsClass
%       - "timeShiftSamplePoints"   Sampling points in time, vector
%       - "nEvalPoints"             N. points for spatial sampling, double
%       - "evalLims"                Limits of spatial sampling, vector
%       - "kernelType"              Chosen kernel, string
%       - "bandwidth"               Width of kernel, double
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        moment1Matrix
        moment2Matrix
        observationData
        evalPoints
        momentOptions
    end
    
    methods
        function obj = MomentClass(moment1Matrix,moment2Matrix,...
                observationData,evalPoints,momentOptions)
            %%%Constructor method, initialise properties
            obj.moment1Matrix = moment1Matrix;
            obj.moment2Matrix = moment2Matrix;
            obj.observationData = observationData;
            obj.evalPoints = evalPoints;
            obj.momentOptions = momentOptions;
        end
        
        function plotMoment1(obj)
            %%%Plot moment 1 for visual inspection
            % Moment grids
            timeShiftVector = obj.observationData.timeStep*...
                obj.momentOptions.timeShiftSamplePoints;
            [~,timeShiftGrid] = meshgrid(obj.evalPoints,timeShiftVector);
            moment1ByTimeShift = obj.moment1Matrix./timeShiftGrid;
            
            % Make figure
            hold on,box on
            p1 = plot(timeShiftVector,moment1ByTimeShift,'k:');
            p2 = scatter(timeShiftGrid(:),moment1ByTimeShift(:),20,...
                'filled','k');
            title('Data for moment M^{(1)}')
            xlabel('Time-shift, \tau')
            ylabel('M^{(1)}/\tau')
            legend([p1(1),p2],...
                {'Same evaluation points','Individual fit data'},...
                'Location','NorthWest')
        end
        function plotMoment2(obj)
            %%%Plot moment 2 for visual inspection
            % Moment grids
            timeShiftVector = obj.observationData.timeStep*...
                obj.momentOptions.timeShiftSamplePoints;
            [~,timeShiftGrid] = meshgrid(obj.evalPoints,timeShiftVector);
            moment2ByTimeShift = obj.moment2Matrix./timeShiftGrid;
            
            % Make figure
            hold on,box on
            p1 = plot(timeShiftVector,moment2ByTimeShift,'k-');
            p2 = scatter(timeShiftGrid(:),moment2ByTimeShift(:),20,...
                'filled','k');
            title('Data for moment M^{(2)}')
            xlabel('Time-shift, \tau')
            ylabel('M^{(2)}/\tau')
            legend([p1(1),p2],...
                {'Same evaluation points','Individual fit data'},...
                'Location','NorthWest')
        end
    end
end