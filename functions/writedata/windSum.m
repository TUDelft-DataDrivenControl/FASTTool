function text = windSum(U, H)
% windSum writes data in a character array, which can be written to file.
% The writing of the character array is separated by this function from
% writing the file in FileIOClass.m, to separate agnostic fileIO from the
% generation of the content of the file.
%
%   Syntax
%     text = windSum(U, H)
%
%   Input arguments
%     U - Average wind speed
%     H = Hub height
%
%   Output
%     text - Character array formatted as an input file that can be written
%     to a plain-text file named wind.sum. This type of file is used when a
%     Bladed style wind file is used in a simulation.
%
%   Behaviour
%     Relevant data from the WindTurbineClass properties is printed to a
%     character array according to the specified format. Several parameter
%     values are not taken from the input data, but are hard coded. This is
%     particularly the case for most model parameters or settings. Some
%     parameters are given a constant or calculated value.
%
%   Called by
%     WindTurbineClass.simulate

text = '';
text = append(text, sprintf('T\tCLOCKWISE\n'));
text = append(text, sprintf('%0.0f\tHUB HEIGHT\n\n', H));
text = append(text, sprintf('%0.3f\tUBAR\n', U));
text = append(text, sprintf('%0.3f\tTI(u)\n', 100));
text = append(text, sprintf('%0.3f\tTI(v)\n', 100));
text = append(text, sprintf('%0.3f\tTI(w)\n\n', 100));
text = append(text, sprintf('0\tHEIGHT OFFSET'));
