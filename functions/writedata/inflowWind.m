function text = inflowWind(type, U, H)
% inflowWind writes data in a character array, which can be written to file. The writing of the character array is separated by this function from writing the file in FileIOClass.m, to separate agnostic fileIO from the generation of the content of the file.
%
%   Syntax
%     text = inflowWind(type, U, H)
%
%   Input arguments
%     type - Type of wind conditions
%     U - 10-Minute average wind speed
%     H - Hub height
%
%   Output
%     text - Character array formatted as an InflowWind v3.01.* input file
%     that can be written to a plain-text file named InflowWind.dat.
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

% ToDo: Comment why these settings have been chosen
if type == 1
    WindType = 1;
    PLexp = 0;
elseif type == 3
    WindType = 1;
    PLexp = 0.2;
else
    WindType = 4;
    PLexp = 0.2;
end

text = '';
text = append(text, sprintf('------- InflowWind v3.01.* INPUT FILE -------------------------------------------------------------------------\n'));
text = append(text, sprintf('12 m/s turbulent winds on 31x31 FF grid and tower for FAST CertTests #18, #19, #21, #22, #23, and #24\n'));
text = append(text, sprintf('---------------------------------------------------------------------------------------------------------------\n'));
text = append(text, sprintf('False         Echo           - Echo input data to <RootName>.ech (flag)\n'));
text = append(text, sprintf('          %i  WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined)\n', WindType));
text = append(text, sprintf('          0   PropagationDir - Direction of wind propagation (meteoroligical rotation from aligned with X (positive rotates towards -Y) -- degrees)\n'));
text = append(text, sprintf('          1   NWindVel       - Number of points to output the wind velocity    (0 to 9)\n'));
text = append(text, sprintf('          0   WindVxiList    - List of coordinates in the inertial X direction (m)\n'));
text = append(text, sprintf('          0   WindVyiList    - List of coordinates in the inertial Y direction (m)\n'));
text = append(text, sprintf('      %5.4f   WindVziList    - List of coordinates in the inertial Z direction (m)\n', H));
text = append(text, sprintf('================== Parameters for Steady Wind Conditions [used only for WindType = 1] =========================\n'));
text = append(text, sprintf('      %5.4f   HWindSpeed     - Horizontal windspeed                            (m/s)\n', U));
text = append(text, sprintf('      %5.4f   RefHt          - Reference height for horizontal wind speed      (m)\n', H));
text = append(text, sprintf('      %5.4f   PLexp          - Power law exponent                              (-)\n', PLexp));
text = append(text, sprintf('================== Parameters for Uniform wind file   [used only for WindType = 2] ============================\n'));
text = append(text, sprintf('"unused"    Filename       - Filename of time series data for uniform wind field.      (-)\n'));
text = append(text, sprintf('         90   RefHt          - Reference height for horizontal wind speed                (m)\n'));
text = append(text, sprintf('     125.88   RefLength      - Reference length for linear horizontal and vertical sheer (-)\n'));
text = append(text, sprintf('================== Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3] ==============\n'));
text = append(text, sprintf('"unused"    Filename       - Name of the Full field wind file to use (.bts)\n'));
text = append(text, sprintf('================== Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4] =========\n'));
text = append(text, sprintf('"wind"    FilenameRoot   - Rootname of the full-field wind file to use (.wnd, .sum)\n'));
text = append(text, sprintf('False         TowerFile      - Have tower file (.twr) (flag)\n'));
text = append(text, sprintf('================== Parameters for HAWC-format binary files  [Only used with WindType = 5] =====================\n'));
text = append(text, sprintf('"unused"    FileName_u     - name of the file containing the u-component fluctuating wind (.bin)\n'));
text = append(text, sprintf('"unused"    FileName_v     - name of the file containing the v-component fluctuating wind (.bin)\n'));
text = append(text, sprintf('"unused"    FileName_w     - name of the file containing the w-component fluctuating wind (.bin)\n'));
text = append(text, sprintf('         64   nx             - number of grids in the x direction (in the 3 files above) (-)\n'));
text = append(text, sprintf('         32   ny             - number of grids in the y direction (in the 3 files above) (-)\n'));
text = append(text, sprintf('         32   nz             - number of grids in the z direction (in the 3 files above) (-)\n'));
text = append(text, sprintf('         16   dx             - distance (in meters) between points in the x direction    (m)\n'));
text = append(text, sprintf('          3   dy             - distance (in meters) between points in the y direction    (m)\n'));
text = append(text, sprintf('          3   dz             - distance (in meters) between points in the z direction    (m)\n'));
text = append(text, sprintf('         90   RefHt          - reference height; the height (in meters) of the vertical center of the grid (m)\n'));
text = append(text, sprintf('  -------------   Scaling parameters for turbulence   ---------------------------------------------------------\n'));
text = append(text, sprintf('          1   ScaleMethod    - Turbulence scaling method   [0 = none, 1 = direct scaling, 2 = calculate scaling factor based on a desired standard deviation]\n'));
text = append(text, sprintf('          1   SFx            - Turbulence scaling factor for the x direction (-)   [ScaleMethod=1]\n'));
text = append(text, sprintf('          1   SFy            - Turbulence scaling factor for the y direction (-)   [ScaleMethod=1]\n'));
text = append(text, sprintf('          1   SFz            - Turbulence scaling factor for the z direction (-)   [ScaleMethod=1]\n'));
text = append(text, sprintf('         12   SigmaFx        - Turbulence standard deviation to calculate scaling from in x direction (m/s)    [ScaleMethod=2]\n'));
text = append(text, sprintf('          8   SigmaFy        - Turbulence standard deviation to calculate scaling from in y direction (m/s)    [ScaleMethod=2]\n'));
text = append(text, sprintf('          2   SigmaFz        - Turbulence standard deviation to calculate scaling from in z direction (m/s)    [ScaleMethod=2]\n'));
text = append(text, sprintf('  -------------   Mean wind profile parameters (added to HAWC-format files)   ---------------------------------\n'));
text = append(text, sprintf('          5   URef           - Mean u-component wind speed at the reference height (m/s)\n'));
text = append(text, sprintf('          2   WindProfile    - Wind profile type (0=constant;1=logarithmic,2=power law)\n'));
text = append(text, sprintf('        0.2   PLExp          - Power law exponent (-) (used for PL wind profile type only)\n'));
text = append(text, sprintf('       0.03   Z0             - Surface roughness length (m) (used for LG wind profile type only)\n'));
text = append(text, sprintf('====================== OUTPUT ==================================================\n'));
text = append(text, sprintf('False         SumPrint     - Print summary data to <RootName>.IfW.sum (flag)\n'));
text = append(text, sprintf('              OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n'));
text = append(text, sprintf('"Wind1VelX"               X-direction wind velocity at point WindList(1)\n'));
text = append(text, sprintf('"Wind1VelY"               Y-direction wind velocity at point WindList(1)\n'));
text = append(text, sprintf('"Wind1VelZ"               Z-direction wind velocity at point WindList(1)\n'));
text = append(text, sprintf('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n'));
text = append(text, sprintf('---------------------------------------------------------------------------------------\n'));
