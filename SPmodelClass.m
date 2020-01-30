classdef SPmodelClass
%Fitted stochastic process model data class
%   William Davis, 20/01/20
%
%   Notes:
%   Defining class for exponentially correlated stochastic process model
%   fit to given data.
%
%   Inputs:
%   - "correlationEstimate"     Correlation time estimate, double
%   - "driftEstimate"           Drift function estimate, vector
%   - "noiseEstimate"           Noise function estimate, vector
%   - "driftInitial"            Drift function initial estimate, vector
%   - "noiseInitial"            Noise function initial estimate, vector
%   - "meanAbsoluteError"       Error between moment data and fit, double
%   - "momentData"              Moment data object, MomentClass
%   - "thetaProperties"         Correlation time properties, structure
%   - "lambdaProperties"        Function fit properties, structure
%   - "fitOptions"              Options, FitOptionsClass
%       - "thetaConvergence"    Iteration error limit for theta, double
%       - "functionConvergence" Same for drift and noise functions, double
%       - "keepObservations"    Keep observationData? logical
%       - "printOutput"         Print to command window? logical
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        correlationEstimate
        driftEstimate
        noiseEstimate
        driftInitial
        noiseInitial
        meanAbsoluteError
        momentData
        thetaProperties
        lambdaProperties
        fitOptions
    end
    
    methods
        function obj = SPmodelClass(correlationEstimate,driftEstimate,...
                noiseEstimate,driftInitial,noiseInitial,...
                meanAbsoluteError,momentData,thetaProperties,...
                lambdaProperties,fitOptions)
            %%%Constructor method, initialise properties
            obj.correlationEstimate = correlationEstimate;
            obj.driftEstimate = driftEstimate;
            obj.noiseEstimate = noiseEstimate;
            obj.driftInitial = driftInitial;
            obj.noiseInitial = noiseInitial;
            obj.meanAbsoluteError = meanAbsoluteError;
            obj.momentData = momentData;
            obj.thetaProperties = thetaProperties;
            obj.lambdaProperties = lambdaProperties;
            obj.fitOptions = fitOptions;
            
            % Keep observations?
            if ~obj.fitOptions.keepObservations
                obj.momentData.observationData = [];
            end
        end
        
        function h = plotDrift(obj)
            %%%Plot drift function for visual inspection
            % Make figure
            hold on,box on
            h.scatter = scatter(obj.momentData.evalPoints,...
                obj.driftEstimate,30,'filled','b');
            title('Estimated drift function')
            xlabel('Position')
            ylabel('Drift function')
        end
        function h = plotNoise(obj)
            %%%Plot noise function for visual inspection
            % Make figure
            hold on,box on
            h.scatter = scatter(obj.momentData.evalPoints,...
                obj.noiseEstimate,30,'filled','g');
            title('Estimated noise function')
            xlabel('Position')
            ylabel('Noise function')
            ylim([0,inf])
        end
        function [l,h] = plotAutocorrelation(obj)
            %%%Plot autocorrelation function and fit for visual inspection
            % Make vectors
            timeShiftVector = obj.momentData.observationData.timeStep*...
                obj.momentData.momentOptions.timeShiftSamplePoints;
            
            timeShiftVectorIncludingZero = [0,timeShiftVector];
            dAfitIncludingZero = [0;obj.thetaProperties.rMatrix*...
                obj.thetaProperties.lambdaStar];
            
            % Make figure
            hold on,box on
            h.scatter = scatter(timeShiftVector,obj.thetaProperties.dA,...
                30,'filled','k');
            h.fit = plot(timeShiftVectorIncludingZero,...
                dAfitIncludingZero,'r-');
            title('Autocorrelation increments')
            xlabel('Time-shift')
            ylabel('Autocorrelation increments')
            l = legend([h.scatter,h.fit],{'Data','Fit function'},...
                'Location','NorthEast');
        end
        function [l,h] = plotMoment1(obj)
            %%%Plot moment 1 and fit for visual inspection
            % Make vectors and grids
            timeShiftVector = obj.momentData.observationData.timeStep*...
                obj.momentData.momentOptions.timeShiftSamplePoints;
            
            [~,timeShiftGrid] = meshgrid(obj.momentData.evalPoints,...
                timeShiftVector);
            
            moment1ByTimeShift = obj.momentData.moment1Matrix./...
                timeShiftGrid;
            
            M1fit = obj.thetaProperties.rMatrix*...
                obj.lambdaProperties.lambda1Star./timeShiftGrid;
            
            % Make figure
            hold on,box on
            h.scatter = scatter(timeShiftGrid(:),moment1ByTimeShift(:),...
                20,'filled','k');
            h.fit = plot(timeShiftGrid,M1fit,'r-');
            title('Data and fit for moment M^{(1)}')
            xlabel('Time-shift')
            ylabel('M^{(1)}/\tau')
            l = legend([h.scatter,h.fit(1)],{'Data','Fit function'},...
                'Location','NorthWest');
        end
        function [l,h] = plotMoment2(obj)
            %%%Plot moment 2 and fit for visual inspection
            % Make vectors and grids
            timeShiftVector = obj.momentData.observationData.timeStep*...
                obj.momentData.momentOptions.timeShiftSamplePoints;
            
            [~,timeShiftGrid] = meshgrid(obj.momentData.evalPoints,...
                timeShiftVector);
            
            moment2ByTimeShift = obj.momentData.moment2Matrix./...
                timeShiftGrid;
            
            M2fit = obj.thetaProperties.rMatrix*...
                obj.lambdaProperties.lambda2Star./timeShiftGrid;
            
            % Make figure
            hold on,box on
            h.scatter = scatter(timeShiftGrid(:),moment2ByTimeShift(:),...
                20,'filled','k');
            h.fit = plot(timeShiftGrid,M2fit,'r-');
            title('Data and fit for moment M^{(2)}')
            xlabel('Time-shift')
            ylabel('M^{(2)}/\tau')
            l = legend([h.scatter,h.fit(1)],{'Data','Fit function'},...
                'Location','NorthWest');
        end
    end
end