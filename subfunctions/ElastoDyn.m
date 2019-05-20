function ElastoDyn(Blade,Tower,Nacelle,Drivetrain,Control,mode,varargin)

%% Operating conditions
RotSpeed = 0;
BlPitch = 0;
GenDOF = 'True';
YawDOF = 'True';

if length(varargin) >= 1
    RotSpeed = varargin{1};
end
if length(varargin) >= 2
    BlPitch = varargin{2};
end
if contains(mode, 'Linearize')
    GenDOF = 'True';
    YawDOF = 'False';
end

LengthScale = sqrt(Control.Torque.SpeedC*(2*pi/60) *  Control.Torque.Demanded * Drivetrain.Generator.Efficiency/5000000);

%% ElastoDyn input file
fid = fopen([pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'ElastoDyn.dat'], 'wt');
fprintf(fid, '------- ELASTODYN v1.03.* INPUT FILE -------------------------------------------\n');
fprintf(fid, 'Created %s.\n', datestr(now));
fprintf(fid, '---------------------- SIMULATION CONTROL --------------------------------------\n');
fprintf(fid, 'False         Echo        - Echo input data to "<RootName>.ech" (flag)\n');
fprintf(fid, '          3   Method      - Integration method: {1: RK4, 2: AB4, or 3: ABM4} (-)\n');
fprintf(fid, '"DEFAULT"     DT          - Integration time step (s)\n');
fprintf(fid, '---------------------- ENVIRONMENTAL CONDITION ---------------------------------\n');
fprintf(fid, '    9.80665   Gravity     - Gravitational acceleration (m/s^2)\n');
fprintf(fid, '---------------------- DEGREES OF FREEDOM --------------------------------------\n');
fprintf(fid, 'True          FlapDOF1    - First flapwise blade mode DOF (flag)\n');
fprintf(fid, 'True          FlapDOF2    - Second flapwise blade mode DOF (flag)\n');
fprintf(fid, 'True          EdgeDOF     - First edgewise blade mode DOF (flag)\n');
fprintf(fid, 'False         TeetDOF     - Rotor-teeter DOF (flag) [unused for 3 blades]\n');
fprintf(fid, 'True          DrTrDOF     - Drivetrain rotational-flexibility DOF (flag)\n');
fprintf(fid, '%s          GenDOF      - Generator DOF (flag)\n', GenDOF);
fprintf(fid, '%s          YawDOF      - Yaw DOF (flag)\n', YawDOF);
fprintf(fid, 'True          TwFADOF1    - First fore-aft tower bending-mode DOF (flag)\n');
fprintf(fid, 'True          TwFADOF2    - Second fore-aft tower bending-mode DOF (flag)\n');
fprintf(fid, 'True          TwSSDOF1    - First side-to-side tower bending-mode DOF (flag)\n');
fprintf(fid, 'True          TwSSDOF2    - Second side-to-side tower bending-mode DOF (flag)\n');
fprintf(fid, 'False         PtfmSgDOF   - Platform horizontal surge translation DOF (flag)\n');
fprintf(fid, 'False         PtfmSwDOF   - Platform horizontal sway translation DOF (flag)\n');
fprintf(fid, 'False         PtfmHvDOF   - Platform vertical heave translation DOF (flag)\n');
fprintf(fid, 'False         PtfmRDOF    - Platform roll tilt rotation DOF (flag)\n');
fprintf(fid, 'False         PtfmPDOF    - Platform pitch tilt rotation DOF (flag)\n');
fprintf(fid, 'False         PtfmYDOF    - Platform yaw rotation DOF (flag)\n');
fprintf(fid, '---------------------- INITIAL CONDITIONS --------------------------------------\n');
fprintf(fid, '          0   OoPDefl     - Initial out-of-plane blade-tip displacement (meters)\n');
fprintf(fid, '          0   IPDefl      - Initial in-plane blade-tip deflection (meters)\n');
fprintf(fid, '      %5.4f   BlPitch(1)  - Blade 1 initial pitch (degrees)\n', BlPitch);
fprintf(fid, '      %5.4f   BlPitch(2)  - Blade 2 initial pitch (degrees)\n', BlPitch);
fprintf(fid, '      %5.4f   BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n', BlPitch);
fprintf(fid, '          0   TeetDefl    - Initial or fixed teeter angle (degrees) [unused for 3 blades]\n');
fprintf(fid, '          0   Azimuth     - Initial azimuth angle for blade 1 (degrees)\n');
fprintf(fid, '      %5.4f   RotSpeed    - Initial or fixed rotor speed (rpm)\n', RotSpeed);
fprintf(fid, '          0   NacYaw      - Initial or fixed nacelle-yaw angle (degrees)\n');
fprintf(fid, '          0   TTDspFA     - Initial fore-aft tower-top displacement (meters)\n');
fprintf(fid, '          0   TTDspSS     - Initial side-to-side tower-top displacement (meters)\n');
fprintf(fid, '          0   PtfmSurge   - Initial or fixed horizontal surge translational displacement of platform (meters)\n');
fprintf(fid, '          0   PtfmSway    - Initial or fixed horizontal sway translational displacement of platform (meters)\n');
fprintf(fid, '          0   PtfmHeave   - Initial or fixed vertical heave translational displacement of platform (meters)\n');
fprintf(fid, '          0   PtfmRoll    - Initial or fixed roll tilt rotational displacement of platform (degrees)\n');
fprintf(fid, '          0   PtfmPitch   - Initial or fixed pitch tilt rotational displacement of platform (degrees)\n');
fprintf(fid, '          0   PtfmYaw     - Initial or fixed yaw rotational displacement of platform (degrees)\n');
fprintf(fid, '---------------------- TURBINE CONFIGURATION -----------------------------------\n');
fprintf(fid, '          %i   NumBl       - Number of blades (-)\n', Blade.Number);
fprintf(fid, ' %5.4f      TipRad      - The distance from the rotor apex to the blade tip (meters)\n', Blade.Radius(end));
fprintf(fid, ' %5.4f      HubRad      - The distance from the rotor apex to the blade root (meters)\n', Blade.Radius(1));
fprintf(fid, ' %5.1f      PreCone(1)  - Blade 1 cone angle (degrees)\n', -Blade.Cone);
fprintf(fid, ' %5.1f      PreCone(2)  - Blade 2 cone angle (degrees)\n', -Blade.Cone);
fprintf(fid, ' %5.1f      PreCone(3)  - Blade 3 cone angle (degrees) [unused for 2 blades]\n', -Blade.Cone);
fprintf(fid, '          0   HubCM       - Distance from rotor apex to hub mass [positive downwind] (meters)\n');
fprintf(fid, '          0   UndSling    - Undersling length [distance from teeter pin to the rotor apex] (meters) [unused for 3 blades]\n');
fprintf(fid, '          0   Delta3      - Delta-3 angle for teetering rotors (degrees) [unused for 3 blades]\n');
fprintf(fid, '          0   AzimB1Up    - Azimuth value to use for I/O when blade 1 points up (degrees)\n');
fprintf(fid, ' %5.5f      OverHang    - Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)\n', -Nacelle.Hub.Overhang);
fprintf(fid, '      1.912   ShftGagL    - Distance from rotor apex [3 blades] or teeter pin [2 blades] to shaft strain gages [positive for upwind rotors] (meters)\n');
fprintf(fid, ' %5.1f      ShftTilt    - Rotor shaft tilt angle (degrees)\n', -Nacelle.Hub.ShaftTilt);
fprintf(fid, ' %5.5f      NacCMxn     - Downwind distance from the tower-top to the nacelle CM (meters)\n', LengthScale*1.9); 
fprintf(fid, '          0   NacCMyn     - Lateral  distance from the tower-top to the nacelle CM (meters)\n');
fprintf(fid, ' %5.2f      NacCMzn     - Vertical distance from the tower-top to the nacelle CM (meters)\n', 0.35*Nacelle.Housing.Diameter);
fprintf(fid, '   -3.09528   NcIMUxn     - Downwind distance from the tower-top to the nacelle IMU (meters)\n');
fprintf(fid, '          0   NcIMUyn     - Lateral  distance from the tower-top to the nacelle IMU (meters)\n');
fprintf(fid, '    2.23336   NcIMUzn     - Vertical distance from the tower-top to the nacelle IMU (meters)\n');
fprintf(fid, ' %5.5f      Twr2Shft    - Vertical distance from the tower-top to the rotor shaft (meters)\n', Tower.HubHeight-Tower.Height(end));
fprintf(fid, ' %5.2f      TowerHt     - Height of tower above ground level [onshore] or MSL [offshore] (meters)\n', Tower.Height(end));
fprintf(fid, '          0   TowerBsHt   - Height of tower base above ground level [onshore] or MSL [offshore] (meters)\n');
fprintf(fid, '          0   PtfmCMxt    - Downwind distance from the ground level [onshore] or MSL [offshore] to the platform CM (meters)\n');
fprintf(fid, '          0   PtfmCMyt    - Lateral distance from the ground level [onshore] or MSL [offshore] to the platform CM (meters)\n');
fprintf(fid, '          0   PtfmCMzt    - Vertical distance from the ground level [onshore] or MSL [offshore] to the platform CM (meters)\n');
fprintf(fid, '          0   PtfmRefzt   - Vertical distance from the ground level [onshore] or MSL [offshore] to the platform reference point (meters)\n');
fprintf(fid, '---------------------- MASS AND INERTIA ----------------------------------------\n');
fprintf(fid, '          0   TipMass(1)  - Tip-brake mass, blade 1 (kg)\n');
fprintf(fid, '          0   TipMass(2)  - Tip-brake mass, blade 2 (kg)\n');
fprintf(fid, '          0   TipMass(3)  - Tip-brake mass, blade 3 (kg) [unused for 2 blades]\n');
fprintf(fid, ' %5.3E      HubMass     - Hub mass (kg)\n', Nacelle.Hub.Mass);
fprintf(fid, ' %5.3E      HubIner     - Hub inertia about rotor axis [3 blades] or teeter axis [2 blades] (kg m^2)\n', (Nacelle.Hub.Mass/56780.0)*(Blade.Radius(1)/1.5)^2*115926.0);
fprintf(fid, ' %5.3f      GenIner     - Generator inertia about HSS (kg m^2)\n', Drivetrain.Generator.HSSInertia);
fprintf(fid, ' %5.3E      NacMass     - Nacelle mass (kg)\n', Nacelle.Housing.Mass);
fprintf(fid, ' %5.3E      NacYIner    - Nacelle inertia about yaw axis (kg m^2)\n', 3.0*Nacelle.Housing.Mass*(LengthScale*1.9)^2); % Factor 3 comes from NacYIner/(NacMass*NacCMxn)^2 of NREL 5 MW turbine
fprintf(fid, '          0   YawBrMass   - Yaw bearing mass (kg)\n');
fprintf(fid, '          0   PtfmMass    - Platform mass (kg)\n');
fprintf(fid, '          0   PtfmRIner   - Platform inertia for roll tilt rotation about the platform CM (kg m^2)\n');
fprintf(fid, '          0   PtfmPIner   - Platform inertia for pitch tilt rotation about the platform CM (kg m^2)\n');
fprintf(fid, '          0   PtfmYIner   - Platform inertia for yaw rotation about the platform CM (kg m^2)\n');
fprintf(fid, '---------------------- BLADE ---------------------------------------------------\n');
fprintf(fid, '         17   BldNodes    - Number of blade nodes (per blade) used for analysis (-)\n');
fprintf(fid, '"ElastoDyn_blade.dat"    BldFile(1)  - Name of file containing properties for blade 1 (quoted string)\n');
fprintf(fid, '"ElastoDyn_blade.dat"    BldFile(2)  - Name of file containing properties for blade 2 (quoted string)\n');
fprintf(fid, '"ElastoDyn_blade.dat"    BldFile(3)  - Name of file containing properties for blade 3 (quoted string) [unused for 2 blades]\n');
fprintf(fid, '---------------------- ROTOR-TEETER --------------------------------------------\n');
fprintf(fid, '          0   TeetMod     - Rotor-teeter spring/damper model {0: none, 1: standard, 2: user-defined from routine UserTeet} (switch) [unused for 3 blades]\n');
fprintf(fid, '          0   TeetDmpP    - Rotor-teeter damper position (degrees) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '          0   TeetDmp     - Rotor-teeter damping constant (N-m/(rad/s)) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '          0   TeetCDmp    - Rotor-teeter rate-independent Coulomb-damping moment (N-m) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '          0   TeetSStP    - Rotor-teeter soft-stop position (degrees) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '          0   TeetHStP    - Rotor-teeter hard-stop position (degrees) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '          0   TeetSSSp    - Rotor-teeter soft-stop linear-spring constant (N-m/rad) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '          0   TeetHSSp    - Rotor-teeter hard-stop linear-spring constant (N-m/rad) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '---------------------- DRIVETRAIN ----------------------------------------------\n');
fprintf(fid, ' %5.1f      GBoxEff     - Gearbox efficiency (%%)\n', 100*Drivetrain.Gearbox.Efficiency);
fprintf(fid, ' %5.1f      GBRatio     - Gearbox ratio (-)\n', Drivetrain.Gearbox.Ratio);
fprintf(fid, '8.67637E+08   DTTorSpr    - Drivetrain torsional spring (N-m/rad)\n');
fprintf(fid, '  6.215E+06   DTTorDmp    - Drivetrain torsional damper (N-m/(rad/s))\n');
fprintf(fid, '---------------------- FURLING -------------------------------------------------\n');
fprintf(fid, 'False         Furling     - Read in additional model properties for furling turbine (flag) [must currently be FALSE)\n');
fprintf(fid, '"unused"      FurlFile    - Name of file containing furling properties (quoted string) [unused when Furling=False]\n');
fprintf(fid, '---------------------- TOWER ---------------------------------------------------\n');
fprintf(fid, '  %i        TwrNodes    - Number of tower nodes used for analysis (-)\n', length(Tower.Height));
fprintf(fid, '"ElastoDyn_tower.dat"    TwrFile     - Name of file containing tower properties (quoted string)\n');
fprintf(fid, '---------------------- OUTPUT --------------------------------------------------\n');
fprintf(fid, 'True          SumPrint    - Print summary data to "<RootName>.sum" (flag)\n');
fprintf(fid, '          1   OutFile     - Switch to determine where output will be placed: {1: in module output file only; 2: in glue code output file only; 3: both} (currently unused)\n');
fprintf(fid, 'True          TabDelim    - Use tab delimiters in text tabular output file? (flag) (currently unused)\n');
fprintf(fid, '"ES10.3E2"    OutFmt      - Format used for text tabular output (except time).  Resulting field should be 10 characters. (quoted string) (currently unused)\n');
fprintf(fid, '          0   TStart      - Time to begin tabular output (s) (currently unused)\n');
fprintf(fid, '          1   DecFact     - Decimation factor for tabular output {1: output every time step} (-) (currently unused)\n');
fprintf(fid, '          0   NTwGages    - Number of tower nodes that have strain gages for output [0 to 9] (-)\n');
fprintf(fid, '         10,         19,         28    TwrGagNd    - List of tower nodes that have strain gages [1 to TwrNodes] (-) [unused if NTwGages=0]\n');
fprintf(fid, '          0   NBlGages    - Number of blade nodes that have strain gages for output [0 to 9] (-)\n');
fprintf(fid, '          5,          9,         13    BldGagNd    - List of blade nodes that have strain gages [1 to BldNodes] (-) [unused if NBlGages=0]\n');
fprintf(fid, '              OutList     - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n');
% fprintf(fid, '"uWind"\n');
% fprintf(fid, '"vWind"\n');
% fprintf(fid, '"wWind"\n');
fprintf(fid, '"OoPDefl1"\n');
fprintf(fid, '"OoPDefl2"\n');
fprintf(fid, '"OoPDefl3"\n');
fprintf(fid, '"IPDefl1"\n');
fprintf(fid, '"IPDefl2"\n');
fprintf(fid, '"IPDefl3"\n');
fprintf(fid, '"NcIMUTAxs"\n');
fprintf(fid, '"NcIMUTAys"\n');
fprintf(fid, '"NcIMUTAzs"\n');
fprintf(fid, '"RootMOoP1"\n');
fprintf(fid, '"RootMOoP2"\n');
fprintf(fid, '"RootMOoP3"\n');
fprintf(fid, '"RootMIP1"\n');
fprintf(fid, '"RootMIP2"\n');
fprintf(fid, '"RootMIP3"\n');
fprintf(fid, '"RootMFlp1"\n');
fprintf(fid, '"RootMFlp2"\n');
fprintf(fid, '"RootMFlp3"\n');
fprintf(fid, '"RootMEdg1"\n');
fprintf(fid, '"RootMEdg2"\n');
fprintf(fid, '"RootMEdg3"\n');
fprintf(fid, '"TwrBsMxt"\n');
fprintf(fid, '"TwrBsMyt"\n');
fprintf(fid, '"TwrBsMzt"\n');
fprintf(fid, '"LSSTipMya"\n');
fprintf(fid, '"LSSTipMza"\n');
fprintf(fid, '"LSSTipVxa"\n');
% fprintf(fid, '"GenPwr"\n');
% fprintf(fid, '"GenTq"\n');
fprintf(fid, '"GenSpeed"\n');
fprintf(fid, '"BlPitch1"\n');
fprintf(fid, '"Azimuth"\n');
fprintf(fid, '"RotPwr"\n');
fprintf(fid, '"RotThrust"\n');
fprintf(fid, '"RotSpeed"\n');
fprintf(fid, '"HSShftTq"\n');
% fprintf(fid, '"OoPDefl1"                - Blade 1 out-of-plane and in-plane deflections and tip twist\n');
% fprintf(fid, '"OoPDefl2"                - Blade 1 out-of-plane and in-plane deflections and tip twist\n');
% fprintf(fid, '"OoPDefl3"                - Blade 1 out-of-plane and in-plane deflections and tip twist\n');
% fprintf(fid, '"IPDefl1"                 - Blade 1 out-of-plane and in-plane deflections and tip twist\n');
% fprintf(fid, '"IPDefl2"                 - Blade 1 out-of-plane and in-plane deflections and tip twist\n');
% fprintf(fid, '"IPDefl3"                 - Blade 1 out-of-plane and in-plane deflections and tip twist\n');
% fprintf(fid, '"TwstDefl1"               - Blade 1 out-of-plane and in-plane deflections and tip twist\n');
% fprintf(fid, '"BldPitch1"               - Blade 1 pitch angle\n');
% fprintf(fid, '"Azimuth"                 - Blade 1 azimuth angle\n');
% fprintf(fid, '"RotSpeed"                - Low-speed shaft and high-speed shaft speeds\n');
% fprintf(fid, '"GenSpeed"                - Low-speed shaft and high-speed shaft speeds\n');
% fprintf(fid, '"TTDspFA"                 - Tower fore-aft and side-to-side displacements and top twist\n');
% fprintf(fid, '"TTDspSS"                 - Tower fore-aft and side-to-side displacements and top twist\n');
% fprintf(fid, '"TTDspTwst"               - Tower fore-aft and side-to-side displacements and top twist\n');
% fprintf(fid, '"Spn2MLxb1"               - Blade 1 local edgewise and flapwise bending moments at span station 2 (approx. 50%% span)\n');
% fprintf(fid, '"Spn2MLyb1"               - Blade 1 local edgewise and flapwise bending moments at span station 2 (approx. 50%% span)\n');
% fprintf(fid, '"RootFxb1"                - Out-of-plane shear, in-plane shear, and axial forces at the root of blade 1\n');
% fprintf(fid, '"RootFyb1"                - Out-of-plane shear, in-plane shear, and axial forces at the root of blade 1\n');
% fprintf(fid, '"RootFzb1"                - Out-of-plane shear, in-plane shear, and axial forces at the root of blade 1\n');
% fprintf(fid, '"RootMxb1"                - In-plane bending, out-of-plane bending, and pitching moments at the root of blade 1\n');
% fprintf(fid, '"RootMyb1"                - In-plane bending, out-of-plane bending, and pitching moments at the root of blade 1\n');
% fprintf(fid, '"RootMzb1"                - In-plane bending, out-of-plane bending, and pitching moments at the root of blade 1\n');
% fprintf(fid, '"RotTorq"                 - Rotor torque and low-speed shaft 0- and 90-bending moments at the main bearing\n');
% fprintf(fid, '"LSSGagMya"               - Rotor torque and low-speed shaft 0- and 90-bending moments at the main bearing\n');
% fprintf(fid, '"LSSGagMza"               - Rotor torque and low-speed shaft 0- and 90-bending moments at the main bearing\n');
% fprintf(fid, '"YawBrFxp"                - Fore-aft shear, side-to-side shear, and vertical forces at the top of the tower (not rotating with nacelle yaw)\n');
% fprintf(fid, '"YawBrFyp"                - Fore-aft shear, side-to-side shear, and vertical forces at the top of the tower (not rotating with nacelle yaw)\n');
% fprintf(fid, '"YawBrFzp"                - Fore-aft shear, side-to-side shear, and vertical forces at the top of the tower (not rotating with nacelle yaw)\n');
% fprintf(fid, '"YawBrMxp"                - Side-to-side bending, fore-aft bending, and yaw moments at the top of the tower (not rotating with nacelle yaw)\n');
% fprintf(fid, '"YawBrMyp"                - Side-to-side bending, fore-aft bending, and yaw moments at the top of the tower (not rotating with nacelle yaw)\n');
% fprintf(fid, '"YawBrMzp"                - Side-to-side bending, fore-aft bending, and yaw moments at the top of the tower (not rotating with nacelle yaw)\n');
% fprintf(fid, '"TwrBsFxt"                - Fore-aft shear, side-to-side shear, and vertical forces at the base of the tower (mudline)\n');
% fprintf(fid, '"TwrBsFyt"                - Fore-aft shear, side-to-side shear, and vertical forces at the base of the tower (mudline)\n');
% fprintf(fid, '"TwrBsFzt"                - Fore-aft shear, side-to-side shear, and vertical forces at the base of the tower (mudline)\n');
% fprintf(fid, '"TwrBsMxt"                - Side-to-side bending, fore-aft bending, and yaw moments at the base of the tower (mudline)\n');
% fprintf(fid, '"TwrBsMyt"                - Side-to-side bending, fore-aft bending, and yaw moments at the base of the tower (mudline)\n');
% fprintf(fid, '"TwrBsMzt"                - Side-to-side bending, fore-aft bending, and yaw moments at the base of the tower (mudline)\n');
fprintf(fid, 'END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n');
fprintf(fid, '---------------------------------------------------------------------------------------\n');
fclose(fid);

%% Blade input file
fid = fopen([pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'ElastoDyn_blade.dat'], 'wt');
fprintf(fid, '------- ELASTODYN V1.00.* INDIVIDUAL BLADE INPUT FILE --------------------------\n');
fprintf(fid, 'Created %s.\n', datestr(now));
fprintf(fid, '---------------------- BLADE PARAMETERS ----------------------------------------\n');
fprintf(fid, '  %i        NBlInpSt    - Number of blade input stations (-)\n', length(Blade.Radius));
fprintf(fid, '   0.477465   BldFlDmp(1) - Blade flap mode #1 structural damping in percent of critical (%%)\n');
fprintf(fid, '   0.477465   BldFlDmp(2) - Blade flap mode #2 structural damping in percent of critical (%%)\n');
fprintf(fid, '   0.477465   BldEdDmp(1) - Blade edge mode #1 structural damping in percent of critical (%%)\n');
fprintf(fid, '---------------------- BLADE ADJUSTMENT FACTORS --------------------------------\n');
fprintf(fid, '          1   FlStTunr(1) - Blade flapwise modal stiffness tuner, 1st mode (-)\n');
fprintf(fid, '          1   FlStTunr(2) - Blade flapwise modal stiffness tuner, 2nd mode (-)\n');
fprintf(fid, '   1.057344   AdjBlMs     - Factor to adjust blade mass density (-)  !bjj: value for AD14=1.04536; value for AD15=1.057344 (it would be nice to enter the requested blade mass instead of a factor here)\n');
fprintf(fid, '          1   AdjFlSt     - Factor to adjust blade flap stiffness (-)\n');
fprintf(fid, '          1   AdjEdSt     - Factor to adjust blade edge stiffness (-)\n');
fprintf(fid, '---------------------- DISTRIBUTED BLADE PROPERTIES ----------------------------\n');
fprintf(fid, '    BlFract      PitchAxis      StrcTwst       BMassDen        FlpStff        EdgStff\n');
fprintf(fid, '      (-)           (-)          (deg)          (kg/m)         (Nm^2)         (Nm^2)\n');
for i = 1:length(Blade.Radius)
    fprintf(fid, '%7.7E\t%7.7E\t%7.7E\t%7.7E\t%7.7E\t%7.7E\n', ...
        (Blade.Radius(i)-Blade.Radius(1))/(Blade.Radius(end)-Blade.Radius(1)), ...
        Blade.PitchAxis(i), ...
        Blade.Twist(i), ...
        Blade.Mass(i), ...
        Blade.EIflap(i), ...
        Blade.EIedge(i));
end
fprintf(fid, '---------------------- BLADE MODE SHAPES ---------------------------------------\n');
fprintf(fid, '    %9.4f   BldFl1Sh(2) - Flap mode 1, coeff of x^2\n', Blade.Flap1_coeff(5));
fprintf(fid, '    %9.4f   BldFl1Sh(3) -            , coeff of x^3\n', Blade.Flap1_coeff(4));
fprintf(fid, '    %9.4f   BldFl1Sh(4) -            , coeff of x^4\n', Blade.Flap1_coeff(3));
fprintf(fid, '    %9.4f   BldFl1Sh(5) -            , coeff of x^5\n', Blade.Flap1_coeff(2));
fprintf(fid, '    %9.4f   BldFl1Sh(6) -            , coeff of x^6\n', Blade.Flap1_coeff(1));
fprintf(fid, '    %9.4f   BldFl2Sh(2) - Flap mode 2, coeff of x^2\n', Blade.Flap2_coeff(5));
fprintf(fid, '    %9.4f   BldFl2Sh(3) -            , coeff of x^3\n', Blade.Flap2_coeff(4));
fprintf(fid, '    %9.4f   BldFl2Sh(4) -            , coeff of x^4\n', Blade.Flap2_coeff(3));
fprintf(fid, '    %9.4f   BldFl2Sh(5) -            , coeff of x^5\n', Blade.Flap2_coeff(2));
fprintf(fid, '    %9.4f   BldFl2Sh(6) -            , coeff of x^6\n', Blade.Flap2_coeff(1));
fprintf(fid, '    %9.4f   BldEdgSh(2) - Edge mode 1, coeff of x^2\n', Blade.Edge1_coeff(5));
fprintf(fid, '    %9.4f   BldEdgSh(3) -            , coeff of x^3\n', Blade.Edge1_coeff(4));
fprintf(fid, '    %9.4f   BldEdgSh(4) -            , coeff of x^4\n', Blade.Edge1_coeff(3));
fprintf(fid, '    %9.4f   BldEdgSh(5) -            , coeff of x^5\n', Blade.Edge1_coeff(2));
fprintf(fid, '    %9.4f   BldEdgSh(6) -            , coeff of x^6\n', Blade.Edge1_coeff(1));
fclose(fid);

%% Tower input file
fid = fopen([pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'ElastoDyn_tower.dat'], 'wt');
fprintf(fid, '------- ELASTODYN V1.00.* TOWER INPUT FILE -------------------------------------\n');
fprintf(fid, 'Created %s.\n', datestr(now));
fprintf(fid, '---------------------- TOWER PARAMETERS ----------------------------------------\n');
fprintf(fid, '  %i        NTwInpSt    - Number of input stations to specify tower geometry\n', length(Tower.Height));
fprintf(fid, '          1   TwrFADmp(1) - Tower 1st fore-aft mode structural damping ratio (%%)\n');
fprintf(fid, '          1   TwrFADmp(2) - Tower 2nd fore-aft mode structural damping ratio (%%)\n');
fprintf(fid, '          1   TwrSSDmp(1) - Tower 1st side-to-side mode structural damping ratio (%%)\n');
fprintf(fid, '          1   TwrSSDmp(2) - Tower 2nd side-to-side mode structural damping ratio (%%)\n');
fprintf(fid, '---------------------- TOWER ADJUSTMUNT FACTORS --------------------------------\n');
fprintf(fid, '          1   FAStTunr(1) - Tower fore-aft modal stiffness tuner, 1st mode (-)\n');
fprintf(fid, '          1   FAStTunr(2) - Tower fore-aft modal stiffness tuner, 2nd mode (-)\n');
fprintf(fid, '          1   SSStTunr(1) - Tower side-to-side stiffness tuner, 1st mode (-)\n');
fprintf(fid, '          1   SSStTunr(2) - Tower side-to-side stiffness tuner, 2nd mode (-)\n');
fprintf(fid, '          1   AdjTwMa     - Factor to adjust tower mass density (-)\n');
fprintf(fid, '          1   AdjFASt     - Factor to adjust tower fore-aft stiffness (-)\n');
fprintf(fid, '          1   AdjSSSt     - Factor to adjust tower side-to-side stiffness (-)\n');
fprintf(fid, '---------------------- DISTRIBUTED TOWER PROPERTIES ----------------------------\n');
fprintf(fid, '  HtFract       TMassDen         TwFAStif       TwSSStif\n');
fprintf(fid, '   (-)           (kg/m)           (Nm^2)         (Nm^2)\n');

% Ensure that tower top < hub height (to prevent negative Twr2Sft)
if Tower.Height(end) >= Tower.HubHeight
    Tower.Height(end) = Tower.HubHeight - 0.1;
end
for i = 1:length(Tower.Height)
    fprintf(fid, '%7.7E\t%7.7E\t%7.7E\t%7.7E\n', ...
        Tower.Height(i)/Tower.Height(end), ...
        Tower.Mass(i), ...
        Tower.EI(i), ...
        Tower.EI(i));
end

fprintf(fid, '---------------------- TOWER FORE-AFT MODE SHAPES ------------------------------\n');
fprintf(fid, '    %9.4f   TwFAM1Sh(2) - Mode 1, coefficient of x^2 term\n', Tower.ForeAft1_coeff(5));
fprintf(fid, '    %9.4f   TwFAM1Sh(3) -       , coefficient of x^3 term\n', Tower.ForeAft1_coeff(4));
fprintf(fid, '    %9.4f   TwFAM1Sh(4) -       , coefficient of x^4 term\n', Tower.ForeAft1_coeff(3));
fprintf(fid, '    %9.4f   TwFAM1Sh(5) -       , coefficient of x^5 term\n', Tower.ForeAft1_coeff(2));
fprintf(fid, '    %9.4f   TwFAM1Sh(6) -       , coefficient of x^6 term\n', Tower.ForeAft1_coeff(1));
fprintf(fid, '    %9.4f   TwFAM2Sh(2) - Mode 2, coefficient of x^2 term\n', Tower.ForeAft2_coeff(5));
fprintf(fid, '    %9.4f   TwFAM2Sh(3) -       , coefficient of x^3 term\n', Tower.ForeAft2_coeff(4));
fprintf(fid, '    %9.4f   TwFAM2Sh(4) -       , coefficient of x^4 term\n', Tower.ForeAft2_coeff(3));
fprintf(fid, '    %9.4f   TwFAM2Sh(5) -       , coefficient of x^5 term\n', Tower.ForeAft2_coeff(2));
fprintf(fid, '    %9.4f   TwFAM2Sh(6) -       , coefficient of x^6 term\n', Tower.ForeAft2_coeff(1));
fprintf(fid, '---------------------- TOWER SIDE-TO-SIDE MODE SHAPES --------------------------\n');
fprintf(fid, '    %9.4f   TwSSM1Sh(2) - Mode 1, coefficient of x^2 term\n', Tower.SideSide1_coeff(5));
fprintf(fid, '    %9.4f   TwSSM1Sh(3) -       , coefficient of x^3 term\n', Tower.SideSide1_coeff(4));
fprintf(fid, '    %9.4f   TwSSM1Sh(4) -       , coefficient of x^4 term\n', Tower.SideSide1_coeff(3));
fprintf(fid, '    %9.4f   TwSSM1Sh(5) -       , coefficient of x^5 term\n', Tower.SideSide1_coeff(2));
fprintf(fid, '    %9.4f   TwSSM1Sh(6) -       , coefficient of x^6 term\n', Tower.SideSide1_coeff(1));
fprintf(fid, '    %9.4f   TwSSM2Sh(2) - Mode 2, coefficient of x^2 term\n', Tower.SideSide2_coeff(5));
fprintf(fid, '    %9.4f   TwSSM2Sh(3) -       , coefficient of x^3 term\n', Tower.SideSide2_coeff(4));
fprintf(fid, '    %9.4f   TwSSM2Sh(4) -       , coefficient of x^4 term\n', Tower.SideSide2_coeff(3));
fprintf(fid, '    %9.4f   TwSSM2Sh(5) -       , coefficient of x^5 term\n', Tower.SideSide2_coeff(2));
fprintf(fid, '    %9.4f   TwSSM2Sh(6) -       , coefficient of x^6 term\n', Tower.SideSide2_coeff(1));
fclose(fid);
