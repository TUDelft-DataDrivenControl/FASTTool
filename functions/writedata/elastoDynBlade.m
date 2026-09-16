function text = elastoDynBlade(Blade)
% elastoDynBlade writes data in a character array, which can be written to
% file. The writing of the character array is separated by this function
% from writing the file in FileIOClass.m, to separate agnostic fileIO from
% the generation of the content of the file.
%
%   Syntax
%     text = elastoDynBlade(Blade)
%
%   Input arguments
%     Blade - Structure with blade specifications
%
%   Output
%     text - Character array formatted as an ELASTODYN V1.00.* blade input
%     file that can be written to a plain-text file named
%     ElastoDyn_blade.dat.
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
text = append(text, sprintf('------- ELASTODYN V1.00.* INDIVIDUAL BLADE INPUT FILE --------------------------\n'));
text = append(text, sprintf('Created %s.\n', datetime));
text = append(text, sprintf('---------------------- BLADE PARAMETERS ----------------------------------------\n'));
text = append(text, sprintf('  %i        NBlInpSt    - Number of blade input stations (-)\n', length(Blade.Radius)));
text = append(text, sprintf('   0.477465   BldFlDmp(1) - Blade flap mode #1 structural damping in percent of critical (%%)\n'));
text = append(text, sprintf('   0.477465   BldFlDmp(2) - Blade flap mode #2 structural damping in percent of critical (%%)\n'));
text = append(text, sprintf('   0.477465   BldEdDmp(1) - Blade edge mode #1 structural damping in percent of critical (%%)\n'));
text = append(text, sprintf('---------------------- BLADE ADJUSTMENT FACTORS --------------------------------\n'));
text = append(text, sprintf('          1   FlStTunr(1) - Blade flapwise modal stiffness tuner, 1st mode (-)\n'));
text = append(text, sprintf('          1   FlStTunr(2) - Blade flapwise modal stiffness tuner, 2nd mode (-)\n'));
text = append(text, sprintf('   1.057344   AdjBlMs     - Factor to adjust blade mass density (-)  !bjj: value for AD14=1.04536; value for AD15=1.057344 (it would be nice to enter the requested blade mass instead of a factor here)\n'));
text = append(text, sprintf('          1   AdjFlSt     - Factor to adjust blade flap stiffness (-)\n'));
text = append(text, sprintf('          1   AdjEdSt     - Factor to adjust blade edge stiffness (-)\n'));
text = append(text, sprintf('---------------------- DISTRIBUTED BLADE PROPERTIES ----------------------------\n'));
text = append(text, sprintf('    BlFract      PitchAxis      StrcTwst       BMassDen        FlpStff        EdgStff\n'));
text = append(text, sprintf('      (-)           (-)          (deg)          (kg/m)         (Nm^2)         (Nm^2)\n'));
for i = 1:length(Blade.Radius)
    text = append(text, sprintf('%7.7E\t%7.7E\t%7.7E\t%7.7E\t%7.7E\t%7.7E\n', ...
        (Blade.Radius(i)-Blade.Radius(1))/(Blade.Radius(end)-Blade.Radius(1)), ...
        Blade.PitchAxis(i), ...
        Blade.Twist(i), ...
        Blade.Mass(i), ...
        Blade.EIflap(i), ...
        Blade.EIedge(i)));
end
text = append(text, sprintf('---------------------- BLADE MODE SHAPES ---------------------------------------\n'));
text = append(text, sprintf('    %9.4f   BldFl1Sh(2) - Flap mode 1, coeff of x^2\n', Blade.Flap1_coeff(5)));
text = append(text, sprintf('    %9.4f   BldFl1Sh(3) -            , coeff of x^3\n', Blade.Flap1_coeff(4)));
text = append(text, sprintf('    %9.4f   BldFl1Sh(4) -            , coeff of x^4\n', Blade.Flap1_coeff(3)));
text = append(text, sprintf('    %9.4f   BldFl1Sh(5) -            , coeff of x^5\n', Blade.Flap1_coeff(2)));
text = append(text, sprintf('    %9.4f   BldFl1Sh(6) -            , coeff of x^6\n', Blade.Flap1_coeff(1)));
text = append(text, sprintf('    %9.4f   BldFl2Sh(2) - Flap mode 2, coeff of x^2\n', Blade.Flap2_coeff(5)));
text = append(text, sprintf('    %9.4f   BldFl2Sh(3) -            , coeff of x^3\n', Blade.Flap2_coeff(4)));
text = append(text, sprintf('    %9.4f   BldFl2Sh(4) -            , coeff of x^4\n', Blade.Flap2_coeff(3)));
text = append(text, sprintf('    %9.4f   BldFl2Sh(5) -            , coeff of x^5\n', Blade.Flap2_coeff(2)));
text = append(text, sprintf('    %9.4f   BldFl2Sh(6) -            , coeff of x^6\n', Blade.Flap2_coeff(1)));
text = append(text, sprintf('    %9.4f   BldEdgSh(2) - Edge mode 1, coeff of x^2\n', Blade.Edge1_coeff(5)));
text = append(text, sprintf('    %9.4f   BldEdgSh(3) -            , coeff of x^3\n', Blade.Edge1_coeff(4)));
text = append(text, sprintf('    %9.4f   BldEdgSh(4) -            , coeff of x^4\n', Blade.Edge1_coeff(3)));
text = append(text, sprintf('    %9.4f   BldEdgSh(5) -            , coeff of x^5\n', Blade.Edge1_coeff(2)));
text = append(text, sprintf('    %9.4f   BldEdgSh(6) -            , coeff of x^6\n', Blade.Edge1_coeff(1)));
