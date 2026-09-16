function text = aeroDynBlade(Blade)
% aeroDynBlade writes data in a character array, which can be written to
% file. The writing of the character array is separated by this function
% from writing the file in FileIOClass.m, to separate agnostic fileIO from
% the generation of the content of the file.
%
%   Syntax
%     text = aeroDynBlade(Blade)
%
%   Input arguments
%     Blade - Structure with blade specifications
%
%   Output
%     text - Character array formatted as an AERODYN v15.00.* blade input
%     file that can be written to a plain-text file named
%     AeroDyn_blade.dat.
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

% Blade input file
text = '';
text = append(text, sprintf('------- AERODYN v15.00.* BLADE DEFINITION INPUT FILE -------------------------------------\n'));
text = append(text, sprintf('Created %s.\n', datetime));
text = append(text, sprintf('======  Blade Properties =================================================================\n'));
text = append(text, sprintf('         %i   NumBlNds           - Number of blade nodes used in the analysis (-)\n', length(Blade.Radius)));
text = append(text, sprintf('  BlSpn        BlCrvAC        BlSwpAC        BlCrvAng       BlTwist        BlChord          BlAFID\n'));
text = append(text, sprintf('   (m)           (m)            (m)            (deg)         (deg)           (m)              (-)\n'));
Blade.Radius(end) = Blade.Radius(end) - 1e-4;
for i = 1:length(Blade.Radius)
    text = append(text, sprintf('%7.7E\t%7.7E\t%7.7E\t%7.7E\t%7.7E\t%7.7E\t%i\n', ...
        Blade.Radius(i)-Blade.Radius(1), ...
        0, ...
        0, ...
        0, ...
        Blade.Twist(i), ...
        Blade.Chord(i), ...
        Blade.NFoil(i)));
end
