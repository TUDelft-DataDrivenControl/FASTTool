%clear all
clear
close all
addpath(genpath('subfunctions'));

% Create the inputfiles folder if it does not exist
if not(exist(['.' filesep 'subfunctions' filesep 'inputfiles'], 'dir'))
    eval(['mkdir .' filesep 'subfunctions' filesep 'inputfiles'])
end

% FASTTool for the course AE4W09 - Wind Turbine Design
run('WindTurbineDesign.m');