function text = servoDyn(Drivetrain, Control, mode, varargin)
% servoDyn writes data in a character array, which can be written to file. The writing of the character array is separated by this function from writing the file in FileIOClass.m, to separate agnostic fileIO from the generation of the content of the file.
%
%   Syntax
%     text = servoDyn(Drivetrain, Control, mode, varargin)
%
%   Input arguments
%     Drivetrain - Structure with drivetrain data
%     Control - Structure with control specifications
%     mode - Mode of the type of simulation (linearisation or simulation,
%     with simulation as default)
%     varargin{1} - Time of event, if relevant (turn generator off;
%     initiate deployment of HSS brake)
%
%   Output
%     text - Character array formatted as a SERVODYN v1.05.* input file
%     that can be written to a plain-text file named ServoDyn.dat.
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
TPitManS = 9999.9;
TimGenOn = 0;
TimGenOf = 9999.9;
THSSBrDp = 9999.9;
BlPitchF = Control.Pitch.Max;
PitManRat = Control.Pitch.Maxrate;
GenTiStr = 'True';

VS_SlPc = 1e-6;
VS_RtTq =  Control.Torque.Demanded;

if strcmpi(mode,'Linearize')
    PCMode = 0;
    VSContrl = 1;
    HSSBrMode = 0;
    
    VS_RtGnSp = Control.Torque.SpeedC;
    VS_Rgn2K = (pi/30)^2 * Control.Torque.OptGain;
elseif strcmpi(mode,'LinearizeWithNRELSettings')
    % Settings recommended by NREL (FAST manual (2005), p.40, bottom of left-hand side column and top of right-hand side column)
    % However, these regularly cause problems at higher wind speeds
    PCMode = 0;
    VSContrl = 1;
    HSSBrMode = 0;
    
    VS_RtGnSp = 1e-6;
    VS_Rgn2K = 1e-6;
else
    PCMode = 4;
    VSContrl = 4;
    HSSBrMode = 1;
    
    VS_RtGnSp = Control.Torque.SpeedC;
    VS_Rgn2K = (pi/30)^2 * Control.Torque.OptGain;
    
    if str2double(mode) == 1       % Power production
    elseif str2double(mode) == 2   % Power production with fault
        TimGenOf = varargin{1};
        THSSBrDp = TimGenOf + Control.Brake.Delay;
    elseif str2double(mode) == 3   % Startup
        GenTiStr = 'False';
    elseif str2double(mode) == 5   % Emergency shutdown
        THSSBrDp = varargin{1};
    elseif str2double(mode) == 6   % Idling
        TimGenOn = 9999.9;
    elseif str2double(mode) == 7   % Parked
        TimGenOn = 9999.9;
        THSSBrDp = 0;
        Control.Brake.Deploytime = 0;
    end
end

text = '';
text = append(text, sprintf('------- SERVODYN v1.05.* INPUT FILE --------------------------------------------\n'));
text = append(text, sprintf('Created %s.\n', datetime));
text = append(text, sprintf('---------------------- SIMULATION CONTROL --------------------------------------\n'));
text = append(text, sprintf('False         Echo         - Echo input data to <RootName>.ech (flag)\n'));
text = append(text, sprintf('     %s     DT           - Communication interval for controllers (s) (or "default")\n', num2str(Control.DT)));
text = append(text, sprintf('---------------------- PITCH CONTROL -------------------------------------------\n'));
text = append(text, sprintf('          %i   PCMode       - Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n', PCMode));
text = append(text, sprintf('          0   TPCOn        - Time to enable active pitch control (s) [unused when PCMode=0]\n'));
text = append(text, sprintf('      %5.1f   TPitManS(1)  - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n', TPitManS));
text = append(text, sprintf('      %5.1f   TPitManS(2)  - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n', TPitManS));
text = append(text, sprintf('      %5.1f   TPitManS(3)  - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n', TPitManS));
text = append(text, sprintf('      %5.1f   PitManRat(1) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)\n', PitManRat));
text = append(text, sprintf('      %5.1f   PitManRat(2) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)\n', PitManRat));
text = append(text, sprintf('      %5.1f   PitManRat(3) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]\n', PitManRat));
text = append(text, sprintf('      %5.1f   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n', BlPitchF));
text = append(text, sprintf('      %5.1f   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n', BlPitchF));
text = append(text, sprintf('      %5.1f   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n', BlPitchF));
text = append(text, sprintf('---------------------- GENERATOR AND TORQUE CONTROL ----------------------------\n'));
text = append(text, sprintf('          %i  VSContrl     - Variable-speed control mode {0: none, 1: simple VS, 3: user-defined from routine UserVSCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n', VSContrl));
text = append(text, sprintf('          1   GenModel     - Generator model {1: simple, 2: Thevenin, 3: user-defined from routine UserGen} (switch) [used only when VSContrl=0]\n'));
text = append(text, sprintf('      %5.1f   GenEff      - Generator efficiency [ignored by the Thevenin and user-defined generator models] (percent)\n', 100*Drivetrain.Generator.Efficiency));
text = append(text, sprintf('%s          GenTiStr     - Method to start the generator {T: timed using TimGenOn, F: generator speed using SpdGenOn} (flag)\n', GenTiStr));
text = append(text, sprintf('True          GenTiStp     - Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)\n'));
text = append(text, sprintf('      %5.4f   SpdGenOn     - Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]\n', Control.Torque.SpeedA));
text = append(text, sprintf('      %5.1f   TimGenOn     - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n', TimGenOn));
text = append(text, sprintf('      %5.1f   TimGenOf     - Time to turn off the generator (s) [used only when GenTiStp=True]\n', TimGenOf));
text = append(text, sprintf('---------------------- SIMPLE VARIABLE-SPEED TORQUE CONTROL --------------------\n'));
text = append(text, sprintf('      %7.6f   VS_RtGnSp    - Rated generator speed for simple variable-speed generator control (HSS side) (rpm) [used only when VSContrl=1]\n', VS_RtGnSp));
text = append(text, sprintf('      %7.6f   VS_RtTq      - Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) [used only when VSContrl=1]\n', VS_RtTq));
text = append(text, sprintf('      %7.6f   VS_Rgn2K     - Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) [used only when VSContrl=1]\n', VS_Rgn2K));
text = append(text, sprintf('      %7.6f   VS_SlPc      - Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%%) [used only when VSContrl=1]\n', VS_SlPc));
text = append(text, sprintf('---------------------- SIMPLE INDUCTION GENERATOR ------------------------------\n'));
text = append(text, sprintf('    9999.9    SIG_SlPc     - Rated generator slip percentage (%%) [used only when VSContrl=0 and GenModel=1]\n'));
text = append(text, sprintf('    9999.9    SIG_SySp     - Synchronous (zero-torque) generator speed (rpm) [used only when VSContrl=0 and GenModel=1]\n'));
text = append(text, sprintf('    9999.9    SIG_RtTq     - Rated torque (N-m) [used only when VSContrl=0 and GenModel=1]\n'));
text = append(text, sprintf('    9999.9    SIG_PORt     - Pull-out ratio (Tpullout/Trated) (-) [used only when VSContrl=0 and GenModel=1]\n'));
text = append(text, sprintf('---------------------- THEVENIN-EQUIVALENT INDUCTION GENERATOR -----------------\n'));
text = append(text, sprintf('     9999.9   TEC_Freq     - Line frequency [50 or 60] (Hz) [used only when VSContrl=0 and GenModel=2]\n'));
text = append(text, sprintf('       9998   TEC_NPol     - Number of poles [even integer > 0] (-) [used only when VSContrl=0 and GenModel=2]\n'));
text = append(text, sprintf('     9999.9   TEC_SRes     - Stator resistance (ohms) [used only when VSContrl=0 and GenModel=2]\n'));
text = append(text, sprintf('     9999.9   TEC_RRes     - Rotor resistance (ohms) [used only when VSContrl=0 and GenModel=2]\n'));
text = append(text, sprintf('     9999.9   TEC_VLL      - Line-to-line RMS voltage (volts) [used only when VSContrl=0 and GenModel=2]\n'));
text = append(text, sprintf('9999.9   TEC_SLR      - Stator leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n'));
text = append(text, sprintf('     9999.9   TEC_RLR      - Rotor leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n'));
text = append(text, sprintf('     9999.9   TEC_MR       - Magnetizing reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n'));
text = append(text, sprintf('---------------------- HIGH-SPEED SHAFT BRAKE ----------------------------------\n'));
text = append(text, sprintf('          %i   HSSBrMode    - HSS brake model {0: none, 1: simple, 3: user-defined from routine UserHSSBr, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n', HSSBrMode));
text = append(text, sprintf(' %5.4f      THSSBrDp     - Time to initiate deployment of the HSS brake (s)\n', THSSBrDp));
text = append(text, sprintf(' %5.4f      HSSBrDT     - Time for HSS-brake to reach full deployment once initiated (sec) [used only when HSSBrMode=1]\n', Control.Brake.Deploytime));
text = append(text, sprintf(' %5.1f      HSSBrTqF    - Fully deployed HSS-brake torque (N-m)\n', Control.Brake.Torque));
text = append(text, sprintf('---------------------- NACELLE-YAW CONTROL -------------------------------------\n'));
text = append(text, sprintf('          0   YCMode       - Yaw control mode {0: none, 3: user-defined from routine UserYawCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'));
text = append(text, sprintf('     9999.9   TYCOn        - Time to enable active yaw control (s) [unused when YCMode=0]\n'));
text = append(text, sprintf('          0   YawNeut      - Neutral yaw position--yaw spring force is zero at this yaw (degrees)\n'));
text = append(text, sprintf('9.02832E+09   YawSpr       - Nacelle-yaw spring constant (N-m/rad)\n'));
text = append(text, sprintf('  1.916E+07   YawDamp      - Nacelle-yaw damping constant (N-m/(rad/s))\n'));
text = append(text, sprintf('     9999.9   TYawManS     - Time to start override yaw maneuver and end standard yaw control (s)\n'));
text = append(text, sprintf('          2   YawManRat    - Yaw maneuver rate (in absolute value) (deg/s)\n'));
text = append(text, sprintf('          0   NacYawF      - Final yaw angle for override yaw maneuvers (degrees)\n'));
text = append(text, sprintf('---------------------- TUNED MASS DAMPER ---------------------------------------\n'));
text = append(text, sprintf('False         CompNTMD     - Compute nacelle tuned mass damper {true/false} (flag)\n'));
text = append(text, sprintf('"unused"    NTMDfile     - Name of the file for nacelle tuned mass damper (quoted string) [unused when CompNTMD is false]\n'));
text = append(text, sprintf('False         CompTTMD     - Compute tower tuned mass damper {true/false} (flag)\n'));
text = append(text, sprintf('"unused"    TTMDfile     - Name of the file for tower tuned mass damper (quoted string) [unused when CompTTMD is false]\n'));
text = append(text, sprintf('---------------------- BLADED INTERFACE ---------------------------------------- [used only with Bladed Interface]\n'));
text = append(text, sprintf('"unused"    DLL_FileName - Name/location of the dynamic library {.dll [Windows] or .so [Linux]} in the Bladed-DLL format (-) [used only with Bladed Interface]\n'));
text = append(text, sprintf('"unused"    DLL_InFile   - Name of input file sent to the DLL (-) [used only with Bladed Interface]\n'));
text = append(text, sprintf('"unused"      DLL_ProcName - Name of procedure in DLL to be called (-) [case sensitive; used only with DLL Interface]\n'));
text = append(text, sprintf('"default"     DLL_DT       - Communication interval for dynamic library (s) (or "default") [used only with Bladed Interface]\n'));
text = append(text, sprintf('false         DLL_Ramp     - Whether a linear ramp should be used between DLL_DT time steps [introduces time shift when true] (flag) [used only with Bladed Interface]\n'));
text = append(text, sprintf('     9999.9   BPCutoff     - Cuttoff frequency for low-pass filter on blade pitch from DLL (Hz) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   NacYaw_North - Reference yaw angle of the nacelle when the upwind end points due North (deg) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   Ptch_Cntrl   - Record 28: Use individual pitch control {0: collective pitch; 1: individual pitch control} (switch) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   Ptch_SetPnt  - Record  5: Below-rated pitch angle set-point (deg) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   Ptch_Min     - Record  6: Minimum pitch angle (deg) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   Ptch_Max     - Record  7: Maximum pitch angle (deg) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   PtchRate_Min - Record  8: Minimum pitch rate (most negative value allowed) (deg/s) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   PtchRate_Max - Record  9: Maximum pitch rate  (deg/s) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   Gain_OM      - Record 16: Optimal mode gain (Nm/(rad/s)^2) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   GenSpd_MinOM - Record 17: Minimum generator speed (rpm) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   GenSpd_MaxOM - Record 18: Optimal mode maximum speed (rpm) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   GenSpd_Dem   - Record 19: Demanded generator speed above rated (rpm) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   GenTrq_Dem   - Record 22: Demanded generator torque above rated (Nm) [used only with Bladed Interface]\n'));
text = append(text, sprintf('          0   GenPwr_Dem   - Record 13: Demanded power (W) [used only with Bladed Interface]\n'));
text = append(text, sprintf('---------------------- BLADED INTERFACE TORQUE-SPEED LOOK-UP TABLE -------------\n'));
text = append(text, sprintf('          0   DLL_NumTrq   - Record 26: No. of points in torque-speed look-up table {0 = none and use the optimal mode parameters; nonzero = ignore the optimal mode PARAMETERs by setting Record 16 to 0.0} (-) [used only with Bladed Interface]\n'));
text = append(text, sprintf(' GenSpd_TLU   GenTrq_TLU\n'));
text = append(text, sprintf(' (rpm)          (Nm)\n'));
text = append(text, sprintf('---------------------- OUTPUT --------------------------------------------------\n'));
text = append(text, sprintf('True          SumPrint     - Print summary data to <RootName>.sum (flag) (currently unused)\n'));
text = append(text, sprintf('          1   OutFile      - Switch to determine where output will be placed: {1: in module output file only; 2: in glue code output file only; 3: both} (currently unused)\n'));
text = append(text, sprintf('True          TabDelim     - Use tab delimiters in text tabular output file? (flag) (currently unused)\n'));
text = append(text, sprintf('"ES10.3E2"    OutFmt       - Format used for text tabular output (except time).  Resulting field should be 10 characters. (quoted string) (currently unused)\n'));
text = append(text, sprintf('          0   TStart       - Time to begin tabular output (s) (currently unused)\n'));
text = append(text, sprintf('              OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n'));
text = append(text, sprintf('"GenPwr"                  - Electrical generator power and torque\n'));
text = append(text, sprintf('"GenTq"                   - Electrical generator power and torque\n'));
text = append(text, sprintf('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n'));
text = append(text, sprintf('---------------------------------------------------------------------------------------\n'));
