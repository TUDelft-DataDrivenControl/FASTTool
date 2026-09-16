function text = elastoDyn(Blade, Tower, Nacelle, Drivetrain, Control, mode, varargin)
% elastoDyn writes data in a character array, which can be written to file.
% The writing of the character array is separated by this function from
% writing the file in FileIOClass.m, to separate agnostic fileIO from the
% generation of the content of the file.
%
%   Syntax
%     text = elastoDyn(Blade, Tower, Nacelle, Drivetrain, Control, mode, varargin)
%
%   Input arguments
%     Blade - Structure with blade specifications
%     Tower - Structure with tower specifications
%     Nacelle - Structure with nacelle specifications, used to determine
%     the tower top mass
%     Drivetrain - Structure with drivetrain data
%     Control - Structure with control specifications
%     mode - Mode of the type of simulation (linearisation or simulation,
%     with simulation as default)
%     type - Beam-type specification: blade or tower
%     varargin{1} - Rotor speed used to initialise the simulation
%     varargin{2} - Blade pitch angle used to initialise the simulation
%
%   Output
%     text - Character array formatted as an ELASTODYN v1.03.* main input
%     file that can be written to a plain-text file named ElastoDyn.dat.
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
RotSpeed = 0;
BlPitch = 0;
YawDOF = 'True';

if length(varargin) >= 1
    RotSpeed = varargin{1};
end
if length(varargin) >= 2
    BlPitch = varargin{2};
end
if contains(mode, 'Linearize') % Modes 'Linearize' and 'LinearizeWithNRELSettings'
    YawDOF = 'False';
end

GenDOF = 'True';
DrTrDOF = 'False'; % Drivetrain DOF disabled to avoid problems that might be caused by unrealistic rotor inertia (e.g. not updated when changing from geared to direct drive), in combination with the fixed values for DTTorSpr and DTTorDmp set below
%DrTrDOF = 'True'; % Drivetrain DOF enabled. This was used in older versions of FASTTool. The controller should be designed with a linearisation file made with the same setting for DrTrDOF as used for the simulation, to ensure compatibility
DTTorSpr = '8.67637E+08'; % Value from NREL 5 MW reference turbine. (Probably) Not used with DrTrDOF set to 'False'
DTTorDmp = '  6.215E+06'; % Value from NREL 5 MW reference turbine. (Probably) Not used with DrTrDOF set to 'False'

LengthScale = sqrt(Control.Torque.SpeedC*(2*pi/60) *  Control.Torque.Demanded * Drivetrain.Generator.Efficiency/5000000);

text = '';
text = append(text, sprintf('------- ELASTODYN v1.03.* INPUT FILE -------------------------------------------\n'));
text = append(text, sprintf('Created %s.\n', datetime));
text = append(text, sprintf('---------------------- SIMULATION CONTROL --------------------------------------\n'));
text = append(text, sprintf('False         Echo        - Echo input data to "<RootName>.ech" (flag)\n'));
text = append(text, sprintf('          3   Method      - Integration method: {1: RK4, 2: AB4, or 3: ABM4} (-)\n'));
text = append(text, sprintf('"DEFAULT"     DT          - Integration time step (s)\n'));
text = append(text, sprintf('---------------------- ENVIRONMENTAL CONDITION ---------------------------------\n'));
text = append(text, sprintf('    9.80665   Gravity     - Gravitational acceleration (m/s^2)\n'));
text = append(text, sprintf('---------------------- DEGREES OF FREEDOM --------------------------------------\n'));
text = append(text, sprintf('True          FlapDOF1    - First flapwise blade mode DOF (flag)\n'));
text = append(text, sprintf('True          FlapDOF2    - Second flapwise blade mode DOF (flag)\n'));
text = append(text, sprintf('True          EdgeDOF     - First edgewise blade mode DOF (flag)\n'));
text = append(text, sprintf('False         TeetDOF     - Rotor-teeter DOF (flag) [unused for 3 blades]\n'));
text = append(text, sprintf('%s          DrTrDOF     - Drivetrain rotational-flexibility DOF (flag)\n', DrTrDOF));
text = append(text, sprintf('%s          GenDOF      - Generator DOF (flag)\n', GenDOF));
text = append(text, sprintf('%s          YawDOF      - Yaw DOF (flag)\n', YawDOF));
text = append(text, sprintf('True          TwFADOF1    - First fore-aft tower bending-mode DOF (flag)\n'));
text = append(text, sprintf('True          TwFADOF2    - Second fore-aft tower bending-mode DOF (flag)\n'));
text = append(text, sprintf('True          TwSSDOF1    - First side-to-side tower bending-mode DOF (flag)\n'));
text = append(text, sprintf('True          TwSSDOF2    - Second side-to-side tower bending-mode DOF (flag)\n'));
text = append(text, sprintf('False         PtfmSgDOF   - Platform horizontal surge translation DOF (flag)\n'));
text = append(text, sprintf('False         PtfmSwDOF   - Platform horizontal sway translation DOF (flag)\n'));
text = append(text, sprintf('False         PtfmHvDOF   - Platform vertical heave translation DOF (flag)\n'));
text = append(text, sprintf('False         PtfmRDOF    - Platform roll tilt rotation DOF (flag)\n'));
text = append(text, sprintf('False         PtfmPDOF    - Platform pitch tilt rotation DOF (flag)\n'));
text = append(text, sprintf('False         PtfmYDOF    - Platform yaw rotation DOF (flag)\n'));
text = append(text, sprintf('---------------------- INITIAL CONDITIONS --------------------------------------\n'));
text = append(text, sprintf('          0   OoPDefl     - Initial out-of-plane blade-tip displacement (meters)\n'));
text = append(text, sprintf('          0   IPDefl      - Initial in-plane blade-tip deflection (meters)\n'));
text = append(text, sprintf('      %5.4f   BlPitch(1)  - Blade 1 initial pitch (degrees)\n', BlPitch));
text = append(text, sprintf('      %5.4f   BlPitch(2)  - Blade 2 initial pitch (degrees)\n', BlPitch));
text = append(text, sprintf('      %5.4f   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n', BlPitch));
text = append(text, sprintf('          0   TeetDefl    - Initial or fixed teeter angle (degrees) [unused for 3 blades]\n'));
text = append(text, sprintf('          0   Azimuth     - Initial azimuth angle for blade 1 (degrees)\n'));
text = append(text, sprintf('      %5.4f   RotSpeed    - Initial or fixed rotor speed (rpm)\n', RotSpeed));
text = append(text, sprintf('          0   NacYaw      - Initial or fixed nacelle-yaw angle (degrees)\n'));
text = append(text, sprintf('          0   TTDspFA     - Initial fore-aft tower-top displacement (meters)\n'));
text = append(text, sprintf('          0   TTDspSS     - Initial side-to-side tower-top displacement (meters)\n'));
text = append(text, sprintf('          0   PtfmSurge   - Initial or fixed horizontal surge translational displacement of platform (meters)\n'));
text = append(text, sprintf('          0   PtfmSway    - Initial or fixed horizontal sway translational displacement of platform (meters)\n'));
text = append(text, sprintf('          0   PtfmHeave   - Initial or fixed vertical heave translational displacement of platform (meters)\n'));
text = append(text, sprintf('          0   PtfmRoll    - Initial or fixed roll tilt rotational displacement of platform (degrees)\n'));
text = append(text, sprintf('          0   PtfmPitch   - Initial or fixed pitch tilt rotational displacement of platform (degrees)\n'));
text = append(text, sprintf('          0   PtfmYaw     - Initial or fixed yaw rotational displacement of platform (degrees)\n'));
text = append(text, sprintf('---------------------- TURBINE CONFIGURATION -----------------------------------\n'));
text = append(text, sprintf('          %i   NumBl       - Number of blades (-)\n', Blade.Number));
text = append(text, sprintf(' %5.4f      TipRad      - The distance from the rotor apex to the blade tip (meters)\n', Blade.Radius(end)));
text = append(text, sprintf(' %5.4f      HubRad      - The distance from the rotor apex to the blade root (meters)\n', Blade.Radius(1)));
text = append(text, sprintf(' %5.1f      PreCone(1)  - Blade 1 cone angle (degrees)\n', -Blade.Cone));
text = append(text, sprintf(' %5.1f      PreCone(2)  - Blade 2 cone angle (degrees)\n', -Blade.Cone));
text = append(text, sprintf(' %5.1f      PreCone(3)  - Blade 3 cone angle (degrees) [unused for 2 blades]\n', -Blade.Cone));
text = append(text, sprintf('          0   HubCM       - Distance from rotor apex to hub mass [positive downwind] (meters)\n'));
text = append(text, sprintf('          0   UndSling    - Undersling length [distance from teeter pin to the rotor apex] (meters) [unused for 3 blades]\n'));
text = append(text, sprintf('          0   Delta3      - Delta-3 angle for teetering rotors (degrees) [unused for 3 blades]\n'));
text = append(text, sprintf('          0   AzimB1Up    - Azimuth value to use for I/O when blade 1 points up (degrees)\n'));
text = append(text, sprintf(' %5.5f      OverHang    - Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)\n', -Nacelle.Hub.Overhang));
text = append(text, sprintf('      1.912   ShftGagL    - Distance from rotor apex [3 blades] or teeter pin [2 blades] to shaft strain gages [positive for upwind rotors] (meters)\n'));
text = append(text, sprintf(' %5.1f      ShftTilt    - Rotor shaft tilt angle (degrees)\n', -Nacelle.Hub.ShaftTilt));
text = append(text, sprintf(' %5.5f      NacCMxn     - Downwind distance from the tower-top to the nacelle CM (meters)\n', LengthScale*1.9)); 
text = append(text, sprintf('          0   NacCMyn     - Lateral  distance from the tower-top to the nacelle CM (meters)\n'));
text = append(text, sprintf(' %5.2f      NacCMzn     - Vertical distance from the tower-top to the nacelle CM (meters)\n', 0.35*Nacelle.Housing.Diameter));
text = append(text, sprintf('   -3.09528   NcIMUxn     - Downwind distance from the tower-top to the nacelle IMU (meters)\n'));
text = append(text, sprintf('          0   NcIMUyn     - Lateral  distance from the tower-top to the nacelle IMU (meters)\n'));
text = append(text, sprintf('    2.23336   NcIMUzn     - Vertical distance from the tower-top to the nacelle IMU (meters)\n'));
text = append(text, sprintf(' %5.5f      Twr2Shft    - Vertical distance from the tower-top to the rotor shaft (meters)\n', Tower.HubHeight-Tower.Height(end)));
text = append(text, sprintf(' %5.2f      TowerHt     - Height of tower above ground level [onshore] or MSL [offshore] (meters)\n', Tower.Height(end)));
text = append(text, sprintf('          0   TowerBsHt   - Height of tower base above ground level [onshore] or MSL [offshore] (meters)\n'));
text = append(text, sprintf('          0   PtfmCMxt    - Downwind distance from the ground level [onshore] or MSL [offshore] to the platform CM (meters)\n'));
text = append(text, sprintf('          0   PtfmCMyt    - Lateral distance from the ground level [onshore] or MSL [offshore] to the platform CM (meters)\n'));
text = append(text, sprintf('          0   PtfmCMzt    - Vertical distance from the ground level [onshore] or MSL [offshore] to the platform CM (meters)\n'));
text = append(text, sprintf('          0   PtfmRefzt   - Vertical distance from the ground level [onshore] or MSL [offshore] to the platform reference point (meters)\n'));
text = append(text, sprintf('---------------------- MASS AND INERTIA ----------------------------------------\n'));
text = append(text, sprintf('          0   TipMass(1)  - Tip-brake mass, blade 1 (kg)\n'));
text = append(text, sprintf('          0   TipMass(2)  - Tip-brake mass, blade 2 (kg)\n'));
text = append(text, sprintf('          0   TipMass(3)  - Tip-brake mass, blade 3 (kg) [unused for 2 blades]\n'));
text = append(text, sprintf(' %5.3E      HubMass     - Hub mass (kg)\n', Nacelle.Hub.Mass));
text = append(text, sprintf(' %5.3E      HubIner     - Hub inertia about rotor axis [3 blades] or teeter axis [2 blades] (kg m^2)\n', (Nacelle.Hub.Mass/56780.0)*(Blade.Radius(1)/1.5)^2*115926.0));
text = append(text, sprintf(' %5.3f      GenIner     - Generator inertia about HSS (kg m^2)\n', Drivetrain.Generator.HSSInertia));
text = append(text, sprintf(' %5.3E      NacMass     - Nacelle mass (kg)\n', Nacelle.Housing.Mass));
text = append(text, sprintf(' %5.3E      NacYIner    - Nacelle inertia about yaw axis (kg m^2)\n', 3.0*Nacelle.Housing.Mass*(LengthScale*1.9)^2)); % Factor 3 comes from NacYIner/(NacMass*NacCMxn)^2 of NREL 5 MW turbine
text = append(text, sprintf('          0   YawBrMass   - Yaw bearing mass (kg)\n'));
text = append(text, sprintf('          0   PtfmMass    - Platform mass (kg)\n'));
text = append(text, sprintf('          0   PtfmRIner   - Platform inertia for roll tilt rotation about the platform CM (kg m^2)\n'));
text = append(text, sprintf('          0   PtfmPIner   - Platform inertia for pitch tilt rotation about the platform CM (kg m^2)\n'));
text = append(text, sprintf('          0   PtfmYIner   - Platform inertia for yaw rotation about the platform CM (kg m^2)\n'));
text = append(text, sprintf('---------------------- BLADE ---------------------------------------------------\n'));
text = append(text, sprintf('         17   BldNodes    - Number of blade nodes (per blade) used for analysis (-)\n'));
text = append(text, sprintf('"ElastoDyn_blade.dat"    BldFile(1)  - Name of file containing properties for blade 1 (quoted string)\n'));
text = append(text, sprintf('"ElastoDyn_blade.dat"    BldFile(2)  - Name of file containing properties for blade 2 (quoted string)\n'));
text = append(text, sprintf('"ElastoDyn_blade.dat"    BldFile(3)  - Name of file containing properties for blade 3 (quoted string) [unused for 2 blades]\n'));
text = append(text, sprintf('---------------------- ROTOR-TEETER --------------------------------------------\n'));
text = append(text, sprintf('          0   TeetMod     - Rotor-teeter spring/damper model {0: none, 1: standard, 2: user-defined from routine UserTeet} (switch) [unused for 3 blades]\n'));
text = append(text, sprintf('          0   TeetDmpP    - Rotor-teeter damper position (degrees) [used only for 2 blades and when TeetMod=1]\n'));
text = append(text, sprintf('          0   TeetDmp     - Rotor-teeter damping constant (N-m/(rad/s)) [used only for 2 blades and when TeetMod=1]\n'));
text = append(text, sprintf('          0   TeetCDmp    - Rotor-teeter rate-independent Coulomb-damping moment (N-m) [used only for 2 blades and when TeetMod=1]\n'));
text = append(text, sprintf('          0   TeetSStP    - Rotor-teeter soft-stop position (degrees) [used only for 2 blades and when TeetMod=1]\n'));
text = append(text, sprintf('          0   TeetHStP    - Rotor-teeter hard-stop position (degrees) [used only for 2 blades and when TeetMod=1]\n'));
text = append(text, sprintf('          0   TeetSSSp    - Rotor-teeter soft-stop linear-spring constant (N-m/rad) [used only for 2 blades and when TeetMod=1]\n'));
text = append(text, sprintf('          0   TeetHSSp    - Rotor-teeter hard-stop linear-spring constant (N-m/rad) [used only for 2 blades and when TeetMod=1]\n'));
text = append(text, sprintf('---------------------- DRIVETRAIN ----------------------------------------------\n'));
text = append(text, sprintf(' %5.1f      GBoxEff     - Gearbox efficiency (%%)\n', 100*Drivetrain.Gearbox.Efficiency));
text = append(text, sprintf(' %5.1f      GBRatio     - Gearbox ratio (-)\n', Drivetrain.Gearbox.Ratio));
text = append(text, sprintf('%s   DTTorSpr    - Drivetrain torsional spring (N-m/rad)\n', DTTorSpr));
text = append(text, sprintf('%s   DTTorDmp    - Drivetrain torsional damper (N-m/(rad/s))\n', DTTorDmp));
text = append(text, sprintf('---------------------- FURLING -------------------------------------------------\n'));
text = append(text, sprintf('False         Furling     - Read in additional model properties for furling turbine (flag) [must currently be FALSE)\n'));
text = append(text, sprintf('"unused"      FurlFile    - Name of file containing furling properties (quoted string) [unused when Furling=False]\n'));
text = append(text, sprintf('---------------------- TOWER ---------------------------------------------------\n'));
text = append(text, sprintf('  %i        TwrNodes    - Number of tower nodes used for analysis (-)\n', length(Tower.Height)));
text = append(text, sprintf('"ElastoDyn_tower.dat"    TwrFile     - Name of file containing tower properties (quoted string)\n'));
text = append(text, sprintf('---------------------- OUTPUT --------------------------------------------------\n'));
text = append(text, sprintf('True          SumPrint    - Print summary data to "<RootName>.sum" (flag)\n'));
text = append(text, sprintf('          1   OutFile     - Switch to determine where output will be placed: {1: in module output file only; 2: in glue code output file only; 3: both} (currently unused)\n'));
text = append(text, sprintf('True          TabDelim    - Use tab delimiters in text tabular output file? (flag) (currently unused)\n'));
text = append(text, sprintf('"ES10.3E2"    OutFmt      - Format used for text tabular output (except time).  Resulting field should be 10 characters. (quoted string) (currently unused)\n'));
text = append(text, sprintf('          0   TStart      - Time to begin tabular output (s) (currently unused)\n'));
text = append(text, sprintf('          1   DecFact     - Decimation factor for tabular output {1: output every time step} (-) (currently unused)\n'));
text = append(text, sprintf('          0   NTwGages    - Number of tower nodes that have strain gages for output [0 to 9] (-)\n'));
text = append(text, sprintf('         10,         19,         28    TwrGagNd    - List of tower nodes that have strain gages [1 to TwrNodes] (-) [unused if NTwGages=0]\n'));
text = append(text, sprintf('          0   NBlGages    - Number of blade nodes that have strain gages for output [0 to 9] (-)\n'));
text = append(text, sprintf('          5,          9,         13    BldGagNd    - List of blade nodes that have strain gages [1 to BldNodes] (-) [unused if NBlGages=0]\n'));
text = append(text, sprintf('              OutList     - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n'));
% text = append(text, sprintf('"uWind"\n'));
% text = append(text, sprintf('"vWind"\n'));
% text = append(text, sprintf('"wWind"\n'));
text = append(text, sprintf('"OoPDefl1"\n'));
text = append(text, sprintf('"OoPDefl2"\n'));
text = append(text, sprintf('"OoPDefl3"\n'));
text = append(text, sprintf('"IPDefl1"\n'));
text = append(text, sprintf('"IPDefl2"\n'));
text = append(text, sprintf('"IPDefl3"\n'));
text = append(text, sprintf('"NcIMUTAxs"\n'));
text = append(text, sprintf('"NcIMUTAys"\n'));
text = append(text, sprintf('"NcIMUTAzs"\n'));
text = append(text, sprintf('"RootMOoP1"\n'));
text = append(text, sprintf('"RootMOoP2"\n'));
text = append(text, sprintf('"RootMOoP3"\n'));
text = append(text, sprintf('"RootMIP1"\n'));
text = append(text, sprintf('"RootMIP2"\n'));
text = append(text, sprintf('"RootMIP3"\n'));
text = append(text, sprintf('"RootMFlp1"\n'));
text = append(text, sprintf('"RootMFlp2"\n'));
text = append(text, sprintf('"RootMFlp3"\n'));
text = append(text, sprintf('"RootMEdg1"\n'));
text = append(text, sprintf('"RootMEdg2"\n'));
text = append(text, sprintf('"RootMEdg3"\n'));
text = append(text, sprintf('"TwrBsMxt"\n'));
text = append(text, sprintf('"TwrBsMyt"\n'));
text = append(text, sprintf('"TwrBsMzt"\n'));
text = append(text, sprintf('"LSSTipMya"\n'));
text = append(text, sprintf('"LSSTipMza"\n'));
text = append(text, sprintf('"LSSTipVxa"\n'));
% text = append(text, sprintf('"GenPwr"\n'));
% text = append(text, sprintf('"GenTq"\n'));
text = append(text, sprintf('"GenSpeed"\n'));
text = append(text, sprintf('"BlPitch1"\n'));
text = append(text, sprintf('"Azimuth"\n'));
text = append(text, sprintf('"RotPwr"\n'));
text = append(text, sprintf('"RotThrust"\n'));
text = append(text, sprintf('"RotSpeed"\n'));
text = append(text, sprintf('"HSShftTq"\n'));
% text = append(text, sprintf('"OoPDefl1"                - Blade 1 out-of-plane and in-plane deflections and tip twist\n'));
% text = append(text, sprintf('"OoPDefl2"                - Blade 1 out-of-plane and in-plane deflections and tip twist\n'));
% text = append(text, sprintf('"OoPDefl3"                - Blade 1 out-of-plane and in-plane deflections and tip twist\n'));
% text = append(text, sprintf('"IPDefl1"                 - Blade 1 out-of-plane and in-plane deflections and tip twist\n'));
% text = append(text, sprintf('"IPDefl2"                 - Blade 1 out-of-plane and in-plane deflections and tip twist\n'));
% text = append(text, sprintf('"IPDefl3"                 - Blade 1 out-of-plane and in-plane deflections and tip twist\n'));
% text = append(text, sprintf('"TwstDefl1"               - Blade 1 out-of-plane and in-plane deflections and tip twist\n'));
% text = append(text, sprintf('"BldPitch1"               - Blade 1 pitch angle\n'));
% text = append(text, sprintf('"Azimuth"                 - Blade 1 azimuth angle\n'));
% text = append(text, sprintf('"RotSpeed"                - Low-speed shaft and high-speed shaft speeds\n'));
% text = append(text, sprintf('"GenSpeed"                - Low-speed shaft and high-speed shaft speeds\n'));
% text = append(text, sprintf('"TTDspFA"                 - Tower fore-aft and side-to-side displacements and top twist\n'));
% text = append(text, sprintf('"TTDspSS"                 - Tower fore-aft and side-to-side displacements and top twist\n'));
% text = append(text, sprintf('"TTDspTwst"               - Tower fore-aft and side-to-side displacements and top twist\n'));
% text = append(text, sprintf('"Spn2MLxb1"               - Blade 1 local edgewise and flapwise bending moments at span station 2 (approx. 50%% span)\n'));
% text = append(text, sprintf('"Spn2MLyb1"               - Blade 1 local edgewise and flapwise bending moments at span station 2 (approx. 50%% span)\n'));
% text = append(text, sprintf('"RootFxb1"                - Out-of-plane shear, in-plane shear, and axial forces at the root of blade 1\n'));
% text = append(text, sprintf('"RootFyb1"                - Out-of-plane shear, in-plane shear, and axial forces at the root of blade 1\n'));
% text = append(text, sprintf('"RootFzb1"                - Out-of-plane shear, in-plane shear, and axial forces at the root of blade 1\n'));
% text = append(text, sprintf('"RootMxb1"                - In-plane bending, out-of-plane bending, and pitching moments at the root of blade 1\n'));
% text = append(text, sprintf('"RootMyb1"                - In-plane bending, out-of-plane bending, and pitching moments at the root of blade 1\n'));
% text = append(text, sprintf('"RootMzb1"                - In-plane bending, out-of-plane bending, and pitching moments at the root of blade 1\n'));
% text = append(text, sprintf('"RotTorq"                 - Rotor torque and low-speed shaft 0- and 90-bending moments at the main bearing\n'));
% text = append(text, sprintf('"LSSGagMya"               - Rotor torque and low-speed shaft 0- and 90-bending moments at the main bearing\n'));
% text = append(text, sprintf('"LSSGagMza"               - Rotor torque and low-speed shaft 0- and 90-bending moments at the main bearing\n'));
% text = append(text, sprintf('"YawBrFxp"                - Fore-aft shear, side-to-side shear, and vertical forces at the top of the tower (not rotating with nacelle yaw)\n'));
% text = append(text, sprintf('"YawBrFyp"                - Fore-aft shear, side-to-side shear, and vertical forces at the top of the tower (not rotating with nacelle yaw)\n'));
% text = append(text, sprintf('"YawBrFzp"                - Fore-aft shear, side-to-side shear, and vertical forces at the top of the tower (not rotating with nacelle yaw)\n'));
% text = append(text, sprintf('"YawBrMxp"                - Side-to-side bending, fore-aft bending, and yaw moments at the top of the tower (not rotating with nacelle yaw)\n'));
% text = append(text, sprintf('"YawBrMyp"                - Side-to-side bending, fore-aft bending, and yaw moments at the top of the tower (not rotating with nacelle yaw)\n'));
% text = append(text, sprintf('"YawBrMzp"                - Side-to-side bending, fore-aft bending, and yaw moments at the top of the tower (not rotating with nacelle yaw)\n'));
% text = append(text, sprintf('"TwrBsFxt"                - Fore-aft shear, side-to-side shear, and vertical forces at the base of the tower (mudline)\n'));
% text = append(text, sprintf('"TwrBsFyt"                - Fore-aft shear, side-to-side shear, and vertical forces at the base of the tower (mudline)\n'));
% text = append(text, sprintf('"TwrBsFzt"                - Fore-aft shear, side-to-side shear, and vertical forces at the base of the tower (mudline)\n'));
% text = append(text, sprintf('"TwrBsMxt"                - Side-to-side bending, fore-aft bending, and yaw moments at the base of the tower (mudline)\n'));
% text = append(text, sprintf('"TwrBsMyt"                - Side-to-side bending, fore-aft bending, and yaw moments at the base of the tower (mudline)\n'));
% text = append(text, sprintf('"TwrBsMzt"                - Side-to-side bending, fore-aft bending, and yaw moments at the base of the tower (mudline)\n'));
text = append(text, sprintf('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n'));
text = append(text, sprintf('---------------------------------------------------------------------------------------\n'));
