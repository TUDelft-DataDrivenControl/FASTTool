function [y11_shape, y11_coeff, y11_freq, ...
    y12_shape, y12_coeff, y12_freq, ...
    y21_shape, y21_coeff, y21_freq, ...
    y22_shape, y22_coeff, y22_freq] = bModesOutput(data)
% bModesOutput extracts mode shape and natural-frequency data from a string
% array that has been extracted from an output file of bModes. The
% extraction is separated by this function from reading the file in
% FileIOClass.m, to separate agnostic fileIO from the interpretation of the
% content of the file.
%
%   Syntax
%     [y11_shape, y11_coeff, y11_freq, ...
%         y12_shape, y12_coeff, y12_freq, ...
%         y21_shape, y21_coeff, y21_freq, ...
%         y22_shape, y22_coeff, y22_freq] = bModesOutput(data)
%
%   Input arguments
%     data - String array with content read from bModes output file by the
%     function readBModesOutput in fileIOClass.m
%
%   Output
%     y11_shape - Mode shape of first side-to-side or flapwise mode
%     y11_coeff - Coefficients of polynomial fit to first side-to-side or flapwise mode
%     y11_freq - Natural frequency of first side-to-side or flapwise mode in Hz
%     y12_shape - Mode shape of second side-to-side or flapwise mode
%     y12_coeff - Coefficients of polynomial fit to second side-to-side or flapwise mode
%     y12_freq - Natural frequency of second side-to-side or flapwise mode in Hz
%     y21_shape - Mode shape of first fore-aft or edgewise mode
%     y21_coeff - Coefficients of polynomial fit to first fore-aft or edgewise mode
%     y21_freq - Natural frequency of first fore-aft or edgewise mode in Hz
%     y22_shape - Mode shape of second fore-aft or edgewise mode
%     y22_coeff - Coefficients of polynomial fit to second fore-aft or edgewise mode
%     y22_freq - Natural frequency of second fore-aft or edgewise mode in Hz
%
%   Behaviour
%     The function steps through the data of the 20 mode shapes that are
%     stored in the output file. The natural frequency of the mode is
%     extracted from the header text of the mode. The type of mode shape
%     (side-to-side/flapwise, fore-aft/edgewise or twist) is determined
%     by detecting in which dimension the deformation in the mode shape is
%     largest. For each type of mode, the mode with the lowest natural
%     frequency is identified as the first mode and the one with the next
%     lowest frequency as the second mode. Results for higher modes and for
%     torsional modes are ignored.
%
%   Called by
%     WindTurbineClass.performModalAnalysis
%     WindTurbineClass.linearise
%     WindTurbineClass.simulate

% The extraction of data uses the knowledge that information about 20 modes
% are available in the bModes output file. This number (20) is hardcoded in
% the function bModesInput.m, which generates the bModes input data file.
% The data contains 7 cell arrays for the 6 columns of the bModes output
% file. (The 7th array contains void data.) Each column is divided over the
% 20 modes using the parameter n, which contains the number of elements per
% mode.
% Initialisation of parametes is only done to give them the correct
% dimensions. The values of the initialisation (zero) are meaningless.
n = (length(data{1}) - 1)/20;
mode = zeros(20,1);
freq = zeros(20,1);
y1 = zeros(20,n-2); % Blade flapwise/tower side-to-side deformations of modes
y2 = zeros(20,n-2); % Blade edgewise/tower fore-aft deformations of modes
y3 = zeros(20,n-2); % Torsional deformations of modes
for i = 1:20
    % Extract the natural frequency from the first character cell of the
    % mode.
    header = cell2mat(data{1}((i-1)*n+1));
    freq(i) = str2double(header(end-15:end-4));
    % Extract the deformations of the mode
    for j = 1:n-2
        y1(i,j) = str2double(data{2}((i-1)*n+2+j));
        y2(i,j) = str2double(data{4}((i-1)*n+2+j));
        y3(i,j) = pi/180*str2double(data{6}((i-1)*n+2+j));
    end
    % Establish the type of mode from the dimension for which the
    % deformation in the mode shape is largest. Store the type in the array
    % mode. Later, only data for modes of type 1 and type 2 will be copied
    % to the output parameters.
    % mode = 1 means blade flapwise/tower side-to-side
    % mode = 2 means blade edgewise/tower fore-aft
    % mode = 3 means blade twist/ tower torsion
    if max(abs(y1(i,:))) > max(abs(y2(i,:))) && max(abs(y1(i,:))) > max(abs(y3(i,:)))
        mode(i) = 1;
    elseif max(abs(y2(i,:))) > max(abs(y1(i,:))) && max(abs(y2(i,:))) > max(abs(y3(i,:)))
        mode(i) = 2;
    else
        mode(i) = 3;
    end
end
y1 = y1(:,1:2:end);
y2 = y2(:,1:2:end);
y3 = y3(:,1:2:end);
    
% Copy data of blade flapwise/tower side-to-side modes. The natural
% frequency and shape information is directly copied. The 'coeff' data are
% coefficients of a 6th order polynomial fit.
i = find(mode == 1);
x = linspace(0,1,(n-3)/2+1);
y11_shape = [y1(i(1),:); y2(i(1),:); y3(i(1),:)];
y11_coeff = polyfit(x,y1(i(1),:),6);
y11_coeff = y11_coeff / sum(y11_coeff(1:5));
y11_freq = freq(i(1));
y12_shape = [y1(i(2),:); y2(i(2),:); y3(i(2),:)];
y12_coeff = polyfit(x,y1(i(2),:),6);
y12_coeff = y12_coeff / sum(y12_coeff(1:5));
y12_freq = freq(i(2));

% Copy data of blade edgewise/tower fore-aft modes. The natural
% frequency and shape information is directly copied. The 'coeff' data are
% coefficients of a 6th order polynomial fit.
i = find(mode == 2);
x = linspace(0,1,(n-3)/2+1);
y21_shape = [y1(i(1),:); y2(i(1),:); y3(i(1),:)];
y21_coeff = polyfit(x,y2(i(1),:),6);
y21_coeff = y21_coeff / sum(y21_coeff(1:5));
y21_freq = freq(i(1));
y22_shape = [y1(i(2),:); y2(i(2),:); y3(i(2),:)];
y22_coeff = polyfit(x,y2(i(2),:),6);
y22_coeff = y22_coeff / sum(y22_coeff(1:5));
y22_freq = freq(i(2));
