function text = elastoDynTower(Tower)
% elastoDynTower writes data in a character array, which can be written to
% file. The writing of the character array is separated by this function
% from writing the file in FileIOClass.m, to separate agnostic fileIO from
% the generation of the content of the file.
%
%   Syntax
%     text = elastoDynTower(Tower)
%
%   Input arguments
%     Tower - Structure with tower specifications
%
%   Output
%     text - Character array formatted as an ELASTODYN V1.00.* tower input
%     file that can be written to a plain-text file named
%     ElastoDyn_tower.dat.
%
%   Behaviour
%     Relevant data from the WindTurbineClass properties is printed to a
%     character array according to the specified format. Several parameter
%     values are not taken from the input data, but are hard coded. This is
%     particularly the case for most model parameters or settings. Some
%     parameters are given a constant or calculated value.
%
%   Called by
%     WindTurbineClass.linearise
%     WindTurbineClass.simulate

text = '';
text = append(text, sprintf('------- ELASTODYN V1.00.* TOWER INPUT FILE -------------------------------------\n'));
text = append(text, sprintf('Created %s.\n', datetime));
text = append(text, sprintf('---------------------- TOWER PARAMETERS ----------------------------------------\n'));
text = append(text, sprintf('  %i        NTwInpSt    - Number of input stations to specify tower geometry\n', length(Tower.Height)));
text = append(text, sprintf('          1   TwrFADmp(1) - Tower 1st fore-aft mode structural damping ratio (%%)\n'));
text = append(text, sprintf('          1   TwrFADmp(2) - Tower 2nd fore-aft mode structural damping ratio (%%)\n'));
text = append(text, sprintf('          1   TwrSSDmp(1) - Tower 1st side-to-side mode structural damping ratio (%%)\n'));
text = append(text, sprintf('          1   TwrSSDmp(2) - Tower 2nd side-to-side mode structural damping ratio (%%)\n'));
text = append(text, sprintf('---------------------- TOWER ADJUSTMUNT FACTORS --------------------------------\n'));
text = append(text, sprintf('          1   FAStTunr(1) - Tower fore-aft modal stiffness tuner, 1st mode (-)\n'));
text = append(text, sprintf('          1   FAStTunr(2) - Tower fore-aft modal stiffness tuner, 2nd mode (-)\n'));
text = append(text, sprintf('          1   SSStTunr(1) - Tower side-to-side stiffness tuner, 1st mode (-)\n'));
text = append(text, sprintf('          1   SSStTunr(2) - Tower side-to-side stiffness tuner, 2nd mode (-)\n'));
text = append(text, sprintf('          1   AdjTwMa     - Factor to adjust tower mass density (-)\n'));
text = append(text, sprintf('          1   AdjFASt     - Factor to adjust tower fore-aft stiffness (-)\n'));
text = append(text, sprintf('          1   AdjSSSt     - Factor to adjust tower side-to-side stiffness (-)\n'));
text = append(text, sprintf('---------------------- DISTRIBUTED TOWER PROPERTIES ----------------------------\n'));
text = append(text, sprintf('  HtFract       TMassDen         TwFAStif       TwSSStif\n'));
text = append(text, sprintf('   (-)           (kg/m)           (Nm^2)         (Nm^2)\n'));

% Ensure that tower top < hub height (to prevent negative Twr2Sft)
if Tower.Height(end) >= Tower.HubHeight
    Tower.Height(end) = Tower.HubHeight - 0.1;
end
for i = 1:length(Tower.Height)
    text = append(text, sprintf('%7.7E\t%7.7E\t%7.7E\t%7.7E\n', ...
        Tower.Height(i)/Tower.Height(end), ...
        Tower.Mass(i), ...
        Tower.EI(i), ...
        Tower.EI(i)));
end

text = append(text, sprintf('---------------------- TOWER FORE-AFT MODE SHAPES ------------------------------\n'));
text = append(text, sprintf('    %9.4f   TwFAM1Sh(2) - Mode 1, coefficient of x^2 term\n', Tower.ForeAft1_coeff(5)));
text = append(text, sprintf('    %9.4f   TwFAM1Sh(3) -       , coefficient of x^3 term\n', Tower.ForeAft1_coeff(4)));
text = append(text, sprintf('    %9.4f   TwFAM1Sh(4) -       , coefficient of x^4 term\n', Tower.ForeAft1_coeff(3)));
text = append(text, sprintf('    %9.4f   TwFAM1Sh(5) -       , coefficient of x^5 term\n', Tower.ForeAft1_coeff(2)));
text = append(text, sprintf('    %9.4f   TwFAM1Sh(6) -       , coefficient of x^6 term\n', Tower.ForeAft1_coeff(1)));
text = append(text, sprintf('    %9.4f   TwFAM2Sh(2) - Mode 2, coefficient of x^2 term\n', Tower.ForeAft2_coeff(5)));
text = append(text, sprintf('    %9.4f   TwFAM2Sh(3) -       , coefficient of x^3 term\n', Tower.ForeAft2_coeff(4)));
text = append(text, sprintf('    %9.4f   TwFAM2Sh(4) -       , coefficient of x^4 term\n', Tower.ForeAft2_coeff(3)));
text = append(text, sprintf('    %9.4f   TwFAM2Sh(5) -       , coefficient of x^5 term\n', Tower.ForeAft2_coeff(2)));
text = append(text, sprintf('    %9.4f   TwFAM2Sh(6) -       , coefficient of x^6 term\n', Tower.ForeAft2_coeff(1)));
text = append(text, sprintf('---------------------- TOWER SIDE-TO-SIDE MODE SHAPES --------------------------\n'));
text = append(text, sprintf('    %9.4f   TwSSM1Sh(2) - Mode 1, coefficient of x^2 term\n', Tower.SideSide1_coeff(5)));
text = append(text, sprintf('    %9.4f   TwSSM1Sh(3) -       , coefficient of x^3 term\n', Tower.SideSide1_coeff(4)));
text = append(text, sprintf('    %9.4f   TwSSM1Sh(4) -       , coefficient of x^4 term\n', Tower.SideSide1_coeff(3)));
text = append(text, sprintf('    %9.4f   TwSSM1Sh(5) -       , coefficient of x^5 term\n', Tower.SideSide1_coeff(2)));
text = append(text, sprintf('    %9.4f   TwSSM1Sh(6) -       , coefficient of x^6 term\n', Tower.SideSide1_coeff(1)));
text = append(text, sprintf('    %9.4f   TwSSM2Sh(2) - Mode 2, coefficient of x^2 term\n', Tower.SideSide2_coeff(5)));
text = append(text, sprintf('    %9.4f   TwSSM2Sh(3) -       , coefficient of x^3 term\n', Tower.SideSide2_coeff(4)));
text = append(text, sprintf('    %9.4f   TwSSM2Sh(4) -       , coefficient of x^4 term\n', Tower.SideSide2_coeff(3)));
text = append(text, sprintf('    %9.4f   TwSSM2Sh(5) -       , coefficient of x^5 term\n', Tower.SideSide2_coeff(2)));
text = append(text, sprintf('    %9.4f   TwSSM2Sh(6) -       , coefficient of x^6 term\n', Tower.SideSide2_coeff(1)));
