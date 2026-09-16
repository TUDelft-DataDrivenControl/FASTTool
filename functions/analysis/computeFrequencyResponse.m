function [MagResponse, PhaseResponse] = computeFrequencyResponse(transfer_function, omega)
% computeFrequencyResponse determines magnitude and phase of the response
% from the provided transfer functions at the frequencies in array omega.
%
%   Syntax
%     [MagResponse, PhaseResponse] = computeFrequencyResponse(transfer_function, omega)
%
%   Input arguments
%     transfer_function - Array of transfer functions at different pitch
%     angles
%     omega - Array of frequencies at which the response should be returned
%
%   Output
%     MagResponse - Magnitude of the response at the frequencies in omega
%     PhaseResponse - Phase of the response at the frequencies in omega
%
%   Called by
%     WindTurbineClass.analyseControl

    % The next code determines the magnitude and phase from the frequency
    % response. The function squeeze is used, because the matrix
    % transfer_function contains unnecessary dimensions. These dimensions
    % have been introduced by the function tf, that has produced the
    % transfer functions in computeController.m
    frf = freqresp(transfer_function, omega);
    MagResponse = mag2db(squeeze(abs(frf)))';
    PhaseResponse = angle(squeeze(frf))'*180/pi;
    
