function text = bModesBlade(Blade)
% bModesBlade writes data in a character array, which can be written to
% file. The writing of the character array is separated by this function
% from writing the file in FileIOClass.m, to separate agnostic fileIO from
% the generation of the content of the file.
%
%   Syntax
%     text = bModesBlade(Blade)
%
%   Input arguments
%     Blade - Structure with blade specifications
%
%   Output
%     text - Character array formatted as a blade input file for bModes
%     that can be written to a plain-text file named BModes_blade.dat.
%
%   Behaviour
%     Relevant data from the WindTurbineClass properties is printed to a
%     character array according to the specified format. Several parameter
%     values are not taken from the input data, but are hard coded. This is
%     particularly the case for most model parameters or settings. Some
%     parameters are given a constant or calculated value.
%
%   Called by
%     WindTurbineClass.performModalAnalysis
%     WindTurbineClass.linearise
%     WindTurbineClass.simulate

    n_secs = length(Blade.Radius);
    text = sprintf('Blade section properties\n');
    text = append(text, sprintf('%0.0f\tn_secs:\tnumber of blade or tower sections at which properties are specified (-)\n\n', n_secs));
    text = append(text, sprintf('sec_loc\tstr_tw\ttw_iner\tmass_den\tflp_iner\tedge_iner\tflp_stff\tedge_stff\ttor_stff\taxial_stff\tcg_offst\tsc_offst\ttc_offst\n'));
    text = append(text, sprintf('(-)\t(deg)\t(deg)\t(kg/m)\t\t(kg-m)\t\t(kg-m)\t\t(Nm^2)\t\t(Nm^2)\t\t(Nm^2)\t\t(N)\t\t(m)\t\t(m)\t\t(m)\n'));
    for i = 1:n_secs
        text = append(text, sprintf('%0.5f\t%7.3f\t%7.3f\t%7.2f\t\t%8.2f\t%8.2f\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\n', ...
            (Blade.Radius(i)-Blade.Radius(1))/(Blade.Radius(end)-Blade.Radius(1)), ...
            Blade.Twist(i), ...
            Blade.Twist(i), ...
            Blade.Mass(i), ...
            Blade.FlapIner(i), ...
            Blade.EdgeIner(i), ...
            Blade.EIflap(i), ...
            Blade.EIedge(i), ...
            Blade.GJ(i), ...
            Blade.EA(i), ...
            Blade.cg(i), ...
            Blade.sc(i), ...
            0));
    end
