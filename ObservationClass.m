classdef ObservationClass
%Observation data class
%   William Davis, 17/01/20
%
%   Notes:
%   Defining class for observations from multiple data-sets.
%
%   Inputs:
%   - "timeStep"                Time-step of data, double
%   - "CellData"                Observations, cell
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        dataCell(1,:) cell 
        dataNumber
        dataLength(1,:)
        dataTotalLength
        timeStep(1,1) double {mustBeNumeric,mustBeFinite,mustBeNonempty,...
            mustBeReal}
    end
    
    methods
        function obj = ObservationClass(timeStep,CellData)
            %%%Constructor method, initialise properties
            % Cells of data
            obj.dataCell = CellData;
            
            % Number of observations
            obj.dataNumber = numel(CellData);
            
            % Vector of observation lengths
            obj.dataLength = cellfun(@numel,CellData);
            
            % Total number of points
            obj.dataTotalLength = sum(obj.dataLength);
            
            % Time-step
            obj.timeStep = timeStep;
            
            % Check cell properties
            checkCellProperties(CellData,obj.dataNumber)
        end
        
        function NewObject = horzcat(Obj,AddObject)
            %%%Append new observations to object
            % Update observations cell
            NewObject.dataCell = horzcat(...
                Obj.dataCell,AddObject.dataCell);
            
            % Update number of observations
            NewObject.dataNumber = Obj.dataNumber + ...
                AddObject.dataNumber;
            
            % Update vector of observation lengths
            NewObject.dataLength = [Obj.dataLength,...
                AddObject.dataLength];
        end
        function NewObject = vertcat(obj,AddObject)
            %%%Append new observations to object (same as horzcat)
            NewObject = horzcat(obj,AddObject);
        end
        function plotObservation(obj)
            %%%Plot observations for visual inspection
            
            % Make figure
            for n = 1:obj.dataNumber
                subplot(obj.dataNumber,1,n)
                hold on,box on
                plot(obj.timeStep*(1:obj.dataLength(n)),...
                    obj.dataCell{n})
                title(['Data: observation no. ',num2str(n)])
                xlabel('Time')
                ylabel('Amplitude')
            end
        end
    end
end
function checkCellProperties(CellData,dataNumber)
% Checking cell properties
for n = 1:dataNumber
    VectorData = CellData{n};
    mustBeNumeric(VectorData)
    mustBeFinite(VectorData)
    mustBeNonempty(VectorData)
    mustBeReal(VectorData)
end
end