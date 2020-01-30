classdef BootstrapOptionsClass
%Bootstrap option class
%   William Davis, 27/01/20
%
%   Notes:
%   Defining class for calculation options for bootstrapping uncertainties
%   on correlation time and drift and noise functions. This class has no
%   non-default value necessary properties.
%
%   Inputs:
%   - "blockLength"             Length of blocks, double
%   - "nBootstrapSamples"       Number of bootstrap samples, double
%   - "parallelProcess"         Switch for parallel processing, logical
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        blockLength(1,1) double {mustBeNumeric,mustBeFinite,...
            mustBeNonempty,mustBeReal,mustBePositive,mustBeInteger} = 500
        nBootstrapSamples(1,1) double {mustBeNumeric,mustBeFinite,...
            mustBeNonempty,mustBeReal,mustBeNonnegative,...
            mustBeInteger} = 0
        parallelProcess(1,1) logical = false
    end
    
    methods
        function obj = BootstrapOptionsClass(nBootstrapSamples)
            %%%Constructor method, initialise necessary properties
            % Validate emptyness of vector
            obj.nBootstrapSamples = nBootstrapSamples;
        end
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end