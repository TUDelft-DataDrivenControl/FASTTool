function text = aeroDyn(Blade, Airfoil, Tower, mode, AirDensity)
% aeroDyn writes data in a character array, which can be written to file.
% The writing of the character array is separated by this function from
% writing the file in FileIOClass.m, to separate agnostic fileIO from the
% generation of the content of the file.
%
%   Syntax
%     text = aeroDyn(Blade, Airfoil, Tower, mode, AirDensity)
%
%   Input arguments
%     Blade - Structure with blade specifications
%     Airfoil - Data structure with aerofoil data
%     Tower - Structure with tower specifications
%     mode - Mode of the type of simulation (linearisation or simulation,
%     with simulation as default)
%     AirDensity - Density of air
%
%   Output
%     text - Character array formatted as an AERODYN v15.03.* main input
%     file that can be written to a plain-text file named AeroDyn.dat.
%
%   Behaviour
%     Relevant data from the WindTurbineClass properties is printed to a
%     character array according to the specified format. Several parameter
%     values are not taken from the input data, but are hard coded. This is
%     particularly the case for most model parameters or settings. Some
%     parameters are given a constant or calculated value.
%     If the mode is 'Linearize', the blade aerofoil model is set to a
%     steady model. Otherwise, it defaults to the Beddous-Leishman unsteady
%     model.
%
%   Called by
%     WindTurbineClass.linearise
%     WindTurbineClass.simulate

% ToDo: Comment why these settings have been chosen
if contains(mode,'Linearize')
    AFAeroMod = 1;
else
    AFAeroMod = 2;
end

% AeroDyn input file
text = '';
text = append(text, sprintf('------- AERODYN v15.03.* INPUT FILE ------------------------------------------------\n'));
text = append(text, sprintf('Created %s.\n', datetime));
text = append(text, sprintf('======  General Options  ============================================================================\n'));
text = append(text, sprintf('False         Echo               - Echo the input to "<rootname>.AD.ech"?  (flag)\n'));
text = append(text, sprintf('"default"     DTAero             - Time interval for aerodynamic calculations {or "default"} (s)\n'));
text = append(text, sprintf('          1   WakeMod            - Type of wake/induction model (switch) {0=none, 1=BEMT}\n'));
text = append(text, sprintf('          %i   AFAeroMod          - Type of blade airfoil aerodynamics model (switch) {1=steady model, 2=Beddoes-Leishman unsteady model}\n', AFAeroMod));
text = append(text, sprintf('          0  TwrPotent          - Type tower influence on wind based on potential flow around the tower (switch) {0=none, 1=baseline potential flow, 2=potential flow with Bak correction}\n'));
text = append(text, sprintf('False          TwrShadow          – Calculate tower influence on wind based on downstream tower shadow? (flag)\n'));
text = append(text, sprintf('False           TwrAero            - Calculate tower aerodynamic loads? (flag)\n'));
text = append(text, sprintf('False          FrozenWake         - Assume frozen wake during linearization? (flag) [used only when WakeMod=1 and when linearizing]\n'));
text = append(text, sprintf('======  Environmental Conditions  ===================================================================\n'));
text = append(text, sprintf('      %4.3f   AirDens            - Air density (kg/m^3)\n', AirDensity));
text = append(text, sprintf('  1.464E-05   KinVisc            - Kinematic air viscosity (m^2/s)\n'));
text = append(text, sprintf('        335   SpdSound           - Speed of sound (m/s)\n'));
text = append(text, sprintf('======  Blade-Element/Momentum Theory Options  ====================================================== [used only when WakeMod=1]\n'));
text = append(text, sprintf('          2   SkewMod            - Type of skewed-wake correction model (switch) {1=uncoupled, 2=Pitt/Peters, 3=coupled} [used only when WakeMod=1]\n'));
text = append(text, sprintf('True          TipLoss            - Use the Prandtl tip-loss model? (flag) [used only when WakeMod=1]\n'));
text = append(text, sprintf('True          HubLoss            - Use the Prandtl hub-loss model? (flag) [used only when WakeMod=1]\n'));
text = append(text, sprintf('True          TanInd             - Include tangential induction in BEMT calculations? (flag) [used only when WakeMod=1]\n'));
text = append(text, sprintf('False         AIDrag             - Include the drag term in the axial-induction calculation? (flag) [used only when WakeMod=1]\n'));
text = append(text, sprintf('False         TIDrag             - Include the drag term in the tangential-induction calculation? (flag) [used only when WakeMod=1 and TanInd=TRUE]\n'));
text = append(text, sprintf('"Default"     IndToler           - Convergence tolerance for BEMT nonlinear solve residual equation {or "default"} (-) [used only when WakeMod=1]\n'));
text = append(text, sprintf('        100   MaxIter            - Maximum number of iteration steps (-) [used only when WakeMod=1]\n'));
text = append(text, sprintf('======  Beddoes-Leishman Unsteady Airfoil Aerodynamics Options  ===================================== [used only when AFAeroMod=2]\n'));
text = append(text, sprintf('          3   UAMod              - Unsteady Aero Model Switch (switch) {1=Baseline model (Original), 2=Gonzalez’s variant (changes in Cn,Cc,Cm), 3=Minemma/Pierce variant (changes in Cc and Cm)} [used only when AFAeroMod=2]\n'));
text = append(text, sprintf('True          FLookup            - Flag to indicate whether a lookup for f will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE)); if FALSE S1-S4 must be provided in airfoil input files (flag) [used only when AFAeroMod=2]\n'));
text = append(text, sprintf('======  Airfoil Information =========================================================================\n'));
text = append(text, sprintf('          1   InCol_Alfa         - The column in the airfoil tables that contains the angle of attack (-)\n'));
text = append(text, sprintf('          2   InCol_Cl           - The column in the airfoil tables that contains the lift coefficient (-)\n'));
text = append(text, sprintf('          3   InCol_Cd           - The column in the airfoil tables that contains the drag coefficient (-)\n'));
text = append(text, sprintf('          4   InCol_Cm           - The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-)\n'));
text = append(text, sprintf('          0   InCol_Cpmin        - The column in the airfoil tables that contains the Cpmin coefficient; use zero if there is no Cpmin column (-)\n'));
text = append(text, sprintf('          %i   NumAFfiles         - Number of airfoil files used (-)\n', length(Blade.IFoil)));
text = append(text, sprintf('"%s"                         AFNames            - Airfoil file names (NumAFfiles lines) (quoted strings)\n', 'AeroDyn_Cylinder 1.dat'));
for i = 2:length(Blade.IFoil)
    text = append(text, sprintf('"%s"\n', ['AeroDyn_', Airfoil.Name{Blade.IFoil(i)}, '.dat']));
end
text = append(text, sprintf('======  Rotor/Blade Properties  =====================================================================\n'));
text = append(text, sprintf('True          UseBlCm            - Include aerodynamic pitching moment in calculations?  (flag)\n'));
text = append(text, sprintf('"AeroDyn_blade.dat"    ADBlFile(1)        - Name of file containing distributed aerodynamic properties for Blade #1 (-)\n'));
text = append(text, sprintf('"AeroDyn_blade.dat"    ADBlFile(2)        - Name of file containing distributed aerodynamic properties for Blade #2 (-) [unused if NumBl < 2]\n'));
text = append(text, sprintf('"AeroDyn_blade.dat"    ADBlFile(3)        - Name of file containing distributed aerodynamic properties for Blade #3 (-) [unused if NumBl < 3]\n'));
text = append(text, sprintf('======  Tower Influence and Aerodynamics ============================================================= [used only when TwrPotent/=0, TwrShadow=True, or TwrAero=True]\n'));
text = append(text, sprintf('          %i   NumTwrNds         - Number of tower nodes used in the analysis  (-) [used only when TwrPotent/=0, TwrShadow=True, or TwrAero=True]\n', length(Tower.Height)));
text = append(text, sprintf('TwrElev        TwrDiam        TwrCd\n'));
text = append(text, sprintf('(m)              (m)           (-)\n'));
for i = 1:length(Tower.Height)
    text = append(text, sprintf('%5.4f    %5.4f    %5.4f\n', Tower.Height(i), Tower.Diameter(i), 0.6));
end
text = append(text, sprintf('======  Outputs  ====================================================================================\n'));
text = append(text, sprintf('False         SumPrint            - Generate a summary file listing input options and interpolated properties to "<rootname>.AD.sum"?  (flag)\n'));
text = append(text, sprintf('          0   NBlOuts             - Number of blade node outputs [0 - 9] (-)\n'));
text = append(text, sprintf(' 1, 9, 19     BlOutNd             - Blade nodes whose values will be output  (-)\n'));
text = append(text, sprintf('          0   NTwOuts             - Number of tower node outputs [0 - 9]  (-)\n'));
text = append(text, sprintf(' 1, 2, 3, 4, 5     TwOutNd             - Tower nodes whose values will be output  (-)\n'));
text = append(text, sprintf('                   OutList             - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n'));
text = append(text, sprintf('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n'));
text = append(text, sprintf('---------------------------------------------------------------------------------------\n'));
