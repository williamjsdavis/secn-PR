classdef FitOptionsClass
%Fit option class
%   William Davis, 27/01/20
%
%   Notes:
%   Defining class for calculation options for estimating correlation time
%   and drift and noise functions from conditional moments. This class has
%   no non-default value necessary properties.
%
%   Inputs:
%   - "thetaConvergence"        Iteration error limit for theta, double
%   - "functionConvergence"     Same for drift and noise functions, double
%   - "fixTheta"                Set theta value? logical
%   - "fixThetaValue"           Value of fixed theta, double
%   - "keepObservations"        Keep observationData? logical
%   - "printOutput"             Print to command window? logical
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        thetaConvergence(1,1) double {mustBeNumeric,mustBeFinite,...
            mustBeNonempty,mustBeReal,mustBePositive} = 1E-2
        functionConvergence(1,1) double {mustBeNumeric,mustBeFinite,...
            mustBeNonempty,mustBeReal,mustBePositive} = 0.2
        fixTheta(1,1) logical = false
        fixThetaValue(1,1) double {mustBeNumeric,mustBeFinite,...
            mustBeNonempty,mustBeReal,mustBeNonnegative} = 0.0
        keepObservations(1,1) logical = true
        printOutput(1,1) logical = true
    end
    
    methods
        function obj = FitOptionsClass()
            %%%Constructor method, initialise necessary properties
            
        end
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end