clear all
close all
addpath(genpath('subfunctions'));

% Create the inputfiles folder if it does not exist
if not(exist('.\subfunctions\inputfiles', 'dir'))
    mkdir '.\subfunctions\inputfiles'
end

% FAST GUI for the course AE4W09 - Wind Turbine Design
run('WindTurbineDesign.m');