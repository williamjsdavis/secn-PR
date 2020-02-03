classdef SPbootstrapClass
%Bootstrapped uncertainties stochastic process model data class
%   William Davis, 22/01/20
%
%   Notes:
%   Defining class for exponentially correlated stochastic process model
%   fit to given data.
%
%   Inputs:
%   - "distributions"           Distribution of bootstrap results,
%                                   structure
%   - "standardErrors"          Standard errors on bootstraps, structure
%   - "percentiles95"           Confidence intervals on bootstraps,
%                                   structure
%   - "meanAbsoluteError"       Fit errors for bootstrap moments, vector
%   - "SPmodel"                 Stochastic process model object, 
%                                   SPmodelClass
%   - "bootstrapOptions"        Options for calculations, structure
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        distributions
        standardErrors
        percentiles95
        meanAbsoluteError
        SPmodel
        bootstrapOptions
    end
    
    methods
        function obj = SPbootstrapClass(distributions,standardErrors,...
                percentiles95,meanAbsoluteError,SPmodel,bootstrapOptions)
            %%%Constructor method, initialise properties
            obj.distributions = distributions;
            obj.standardErrors = standardErrors;
            obj.percentiles95 = percentiles95;
            obj.meanAbsoluteError = meanAbsoluteError;
            obj.SPmodel = SPmodel;
            obj.bootstrapOptions = bootstrapOptions;
        end
        
        function h = plotDrift(obj)
            %%%Plot drift function for visual inspection
            % Make figure
            hold on,box on
            h.fill = fillErrorArea(obj.SPmodel.momentData.evalPoints,...
                obj.percentiles95.driftEstimate,[0.9,0.9,0.9]);
            h.scatter = scatter(obj.SPmodel.momentData.evalPoints,...
                obj.SPmodel.driftEstimate,30,'filled','b');
            title('Estimated drift function')
            xlabel('Position')
            ylabel('Drift function')
        end
        function h = plotNoise(obj)
            %%%Plot noise function for visual inspection
            % Make figure
            hold on,box on
            h.fill = fillErrorArea(obj.SPmodel.momentData.evalPoints,...
                obj.percentiles95.noiseEstimate,[0.9,0.9,0.9]);
            h.scatter = scatter(obj.SPmodel.momentData.evalPoints,...
                obj.SPmodel.noiseEstimate,30,'filled','g');
            title('Estimated noise function')
            xlabel('Position')
            ylabel('Noise function')
            ylim([0,inf])
        end
        function [l,h] = plotCorrelationTime(obj)
            %%%Plot correlation time and bootstraps for visual inspection
            % Make figure
            hold on,box on
            h.histogram = histogram(obj.distributions.correlationEstimate);
            h.plot = plot(obj.SPmodel.correlationEstimate*ones(2,1),...
                ylim,'r-','LineWidth',2);
            title('Correlation time and bootstrap results')
            xlabel('Correlation time')
            ylabel('Counts')
            xlim([0,inf])
            l = legend([h.histogram,h.plot],'Bootstrap','Original',...
                'Location','NorthWest');
        end
        function [l,h] = plotMeanFitError(obj)
            %%%Plot mean absolute fit errors for bootetrap results
            % Make figure
            hold on,box on
            h.histogram = histogram(obj.meanAbsoluteError);
            h.plot = plot(obj.SPmodel.meanAbsoluteError.bothMoments*...
                ones(2,1),ylim,'r-','LineWidth',2);
            title('Bootstrap moment fit error')
            xlabel('Mean absolute fit, moments 1 & 2')
            ylabel('Counts')
            xlim([0,inf])
            l = legend([h.histogram,h.plot],'Bootstrap','Original',...
                'Location','NorthWest');
        end
    end
end
function h = fillErrorArea(x,y,color)
% Quick function for filling error areas
h = fill([x,fliplr(x)],[y(1,:),fliplr(y(2,:))],color);
end