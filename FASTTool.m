%clear all
clear
close all
addpath(genpath('subfunctions'));

% Create the inputfiles folder if it does not exist
if not(exist(['.' filesep 'subfunctions' filesep 'inputfiles'], 'dir'))
    eval(['mkdir .' filesep 'subfunctions' filesep 'inputfiles'])
end

% Remove Deny-permissions from the inputfiles-folder, to avoid problems
% with BModes not being able to write its output file
[~,~] = system(['icacls subfunctions' filesep 'inputfiles /inheritance:d']);
[~,~] = system(['icacls subfunctions' filesep 'inputfiles /remove:d Everyone']);

% FASTTool for the course AE4W09 - Wind Turbine Design
WindTurbineDesign_App