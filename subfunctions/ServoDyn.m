function ServoDyn(Drivetrain,Control,mode,varargin)

%% Operating conditions
TPitManS = 9999.9;
TimGenOn = 0;
TimGenOf = 9999.9;
THSSBrDp = 9999.9;
BlPitchF = Control.Pitch.Max;
PitManRat = Control.Pitch.Maxrate;
GenTiStr = 'True';

if strcmpi(mode,'Linearize')
    
    PCMode = 0;
    VSContrl = 1;
    HSSBrMode = 0;
    
else
    
    PCMode = 4;
    VSContrl = 4;
    HSSBrMode = 1;

    if mode == 1       % Power production
    elseif mode == 2   % Power production with fault
        TimGenOf = varargin{1};
        THSSBrDp = TimGenOf + Control.Brake.Delay;
    elseif mode == 3   % Startup
        GenTiStr = 'False';
    elseif mode == 5   % Emergency shutdown
        THSSBrDp = varargin{1};
    elseif mode == 6   % Idling
        TimGenOn = 9999.9;
    elseif mode == 7   % Parked
        TimGenOn = 9999.9;
        THSSBrDp = 0;
        Control.Brake.Deploytime = 0;
    end
    
end

%% ServoDyn input file
fid = fopen([pwd, '\subfunctions\inputfiles\ServoDyn.dat'], 'wt');
fprintf(fid, '------- SERVODYN v1.05.* INPUT FILE --------------------------------------------\n');
fprintf(fid, 'Created %s.\n', datestr(now));
fprintf(fid, '---------------------- SIMULATION CONTROL --------------------------------------\n');
fprintf(fid, 'False         Echo         - Echo input data to <RootName>.ech (flag)\n');
fprintf(fid, '     %s     DT           - Communication interval for controllers (s) (or "default")\n', num2str(Control.DT));
fprintf(fid, '---------------------- PITCH CONTROL -------------------------------------------\n');
fprintf(fid, '          %i   PCMode       - Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n', PCMode);
fprintf(fid, '          0   TPCOn        - Time to enable active pitch control (s) [unused when PCMode=0]\n');
fprintf(fid, '      %5.1f   TPitManS(1)  - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n', TPitManS);
fprintf(fid, '      %5.1f   TPitManS(2)  - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n', TPitManS);
fprintf(fid, '      %5.1f   TPitManS(3)  - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n', TPitManS);
fprintf(fid, '      %5.1f   PitManRat(1) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)\n', PitManRat);
fprintf(fid, '      %5.1f   PitManRat(2) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)\n', PitManRat);
fprintf(fid, '      %5.1f   PitManRat(3) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]\n', PitManRat);
fprintf(fid, '      %5.1f   BlPitchF(1)  - Blade 1 final pitch for pitch maneuvers (degrees)\n', BlPitchF);
fprintf(fid, '      %5.1f   BlPitchF(2)  - Blade 2 final pitch for pitch maneuvers (degrees)\n', BlPitchF);
fprintf(fid, '      %5.1f   BlPitchF(3)  - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n', BlPitchF);
fprintf(fid, '---------------------- GENERATOR AND TORQUE CONTROL ----------------------------\n');
fprintf(fid, '          %i  VSContrl     - Variable-speed control mode {0: none, 1: simple VS, 3: user-defined from routine UserVSCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n', VSContrl);
fprintf(fid, '          1   GenModel     - Generator model {1: simple, 2: Thevenin, 3: user-defined from routine UserGen} (switch) [used only when VSContrl=0]\n');
fprintf(fid, '      %5.1f   GenEff      - Generator efficiency [ignored by the Thevenin and user-defined generator models] (percent)\n', 100*Drivetrain.Generator.Efficiency);
fprintf(fid, '%s          GenTiStr     - Method to start the generator {T: timed using TimGenOn, F: generator speed using SpdGenOn} (flag)\n', GenTiStr);
fprintf(fid, 'True          GenTiStp     - Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)\n');
fprintf(fid, '      %5.4f   SpdGenOn     - Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]\n', Control.Torque.SpeedA);
fprintf(fid, '      %5.1f   TimGenOn     - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n', TimGenOn);
fprintf(fid, '      %5.1f   TimGenOf     - Time to turn off the generator (s) [used only when GenTiStp=True]\n', TimGenOf);
fprintf(fid, '---------------------- SIMPLE VARIABLE-SPEED TORQUE CONTROL --------------------\n');
fprintf(fid, '      %5.4f   VS_RtGnSp    - Rated generator speed for simple variable-speed generator control (HSS side) (rpm) [used only when VSContrl=1]\n', Control.Torque.SpeedC);
fprintf(fid, '      %5.4f   VS_RtTq      - Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) [used only when VSContrl=1]\n', Control.Torque.Demanded);
fprintf(fid, '      %5.4f   VS_Rgn2K     - Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) [used only when VSContrl=1]\n', (pi/30)^2 * Control.Torque.OptGain);
fprintf(fid, ' 9.9999E-06   VS_SlPc      - Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%%) [used only when VSContrl=1]\n');
fprintf(fid, '---------------------- SIMPLE INDUCTION GENERATOR ------------------------------\n');
fprintf(fid, '     9999.9   SIG_SlPc     - Rated generator slip percentage (%%) [used only when VSContrl=0 and GenModel=1]\n');
fprintf(fid, '     9999.9   SIG_SySp     - Synchronous (zero-torque) generator speed (rpm) [used only when VSContrl=0 and GenModel=1]\n');
fprintf(fid, '     9999.9   SIG_RtTq     - Rated torque (N-m) [used only when VSContrl=0 and GenModel=1]\n');
fprintf(fid, '     9999.9   SIG_PORt     - Pull-out ratio (Tpullout/Trated) (-) [used only when VSContrl=0 and GenModel=1]\n');
fprintf(fid, '---------------------- THEVENIN-EQUIVALENT INDUCTION GENERATOR -----------------\n');
fprintf(fid, '     9999.9   TEC_Freq     - Line frequency [50 or 60] (Hz) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '       9998   TEC_NPol     - Number of poles [even integer > 0] (-) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '     9999.9   TEC_SRes     - Stator resistance (ohms) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '     9999.9   TEC_RRes     - Rotor resistance (ohms) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '     9999.9   TEC_VLL      - Line-to-line RMS voltage (volts) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '9999.9   TEC_SLR      - Stator leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '     9999.9   TEC_RLR      - Rotor leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '     9999.9   TEC_MR       - Magnetizing reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '---------------------- HIGH-SPEED SHAFT BRAKE ----------------------------------\n');
fprintf(fid, '          %i   HSSBrMode    - HSS brake model {0: none, 1: simple, 3: user-defined from routine UserHSSBr, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n', HSSBrMode);
fprintf(fid, ' %5.4f      THSSBrDp     - Time to initiate deployment of the HSS brake (s)\n', THSSBrDp);
fprintf(fid, ' %5.4f      HSSBrDT     - Time for HSS-brake to reach full deployment once initiated (sec) [used only when HSSBrMode=1]\n', Control.Brake.Deploytime);
fprintf(fid, ' %5.1f      HSSBrTqF    - Fully deployed HSS-brake torque (N-m)\n', Control.Brake.Torque);
fprintf(fid, '---------------------- NACELLE-YAW CONTROL -------------------------------------\n');
fprintf(fid, '          0   YCMode       - Yaw control mode {0: none, 3: user-defined from routine UserYawCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n');
fprintf(fid, '     9999.9   TYCOn        - Time to enable active yaw control (s) [unused when YCMode=0]\n');
fprintf(fid, '          0   YawNeut      - Neutral yaw position--yaw spring force is zero at this yaw (degrees)\n');
fprintf(fid, '9.02832E+09   YawSpr       - Nacelle-yaw spring constant (N-m/rad)\n');
fprintf(fid, '  1.916E+07   YawDamp      - Nacelle-yaw damping constant (N-m/(rad/s))\n');
fprintf(fid, '     9999.9   TYawManS     - Time to start override yaw maneuver and end standard yaw control (s)\n');
fprintf(fid, '          2   YawManRat    - Yaw maneuver rate (in absolute value) (deg/s)\n');
fprintf(fid, '          0   NacYawF      - Final yaw angle for override yaw maneuvers (degrees)\n');
fprintf(fid, '---------------------- TUNED MASS DAMPER ---------------------------------------\n');
fprintf(fid, 'False         CompNTMD     - Compute nacelle tuned mass damper {true/false} (flag)\n');
fprintf(fid, '"unused"    NTMDfile     - Name of the file for nacelle tuned mass damper (quoted string) [unused when CompNTMD is false]\n');
fprintf(fid, 'False         CompTTMD     - Compute tower tuned mass damper {true/false} (flag)\n');
fprintf(fid, '"unused"    TTMDfile     - Name of the file for tower tuned mass damper (quoted string) [unused when CompTTMD is false]\n');
fprintf(fid, '---------------------- BLADED INTERFACE ---------------------------------------- [used only with Bladed Interface]\n');
fprintf(fid, '"unused"    DLL_FileName - Name/location of the dynamic library {.dll [Windows] or .so [Linux]} in the Bladed-DLL format (-) [used only with Bladed Interface]\n');
fprintf(fid, '"unused"    DLL_InFile   - Name of input file sent to the DLL (-) [used only with Bladed Interface]\n');
fprintf(fid, '"unused"      DLL_ProcName - Name of procedure in DLL to be called (-) [case sensitive; used only with DLL Interface]\n');
fprintf(fid, '"default"     DLL_DT       - Communication interval for dynamic library (s) (or "default") [used only with Bladed Interface]\n');
fprintf(fid, 'false         DLL_Ramp     - Whether a linear ramp should be used between DLL_DT time steps [introduces time shift when true] (flag) [used only with Bladed Interface]\n');
fprintf(fid, '     9999.9   BPCutoff     - Cuttoff frequency for low-pass filter on blade pitch from DLL (Hz) [used only with Bladed Interface]\n');
fprintf(fid, '          0   NacYaw_North - Reference yaw angle of the nacelle when the upwind end points due North (deg) [used only with Bladed Interface]\n');
fprintf(fid, '          0   Ptch_Cntrl   - Record 28: Use individual pitch control {0: collective pitch; 1: individual pitch control} (switch) [used only with Bladed Interface]\n');
fprintf(fid, '          0   Ptch_SetPnt  - Record  5: Below-rated pitch angle set-point (deg) [used only with Bladed Interface]\n');
fprintf(fid, '          0   Ptch_Min     - Record  6: Minimum pitch angle (deg) [used only with Bladed Interface]\n');
fprintf(fid, '          0   Ptch_Max     - Record  7: Maximum pitch angle (deg) [used only with Bladed Interface]\n');
fprintf(fid, '          0   PtchRate_Min - Record  8: Minimum pitch rate (most negative value allowed) (deg/s) [used only with Bladed Interface]\n');
fprintf(fid, '          0   PtchRate_Max - Record  9: Maximum pitch rate  (deg/s) [used only with Bladed Interface]\n');
fprintf(fid, '          0   Gain_OM      - Record 16: Optimal mode gain (Nm/(rad/s)^2) [used only with Bladed Interface]\n');
fprintf(fid, '          0   GenSpd_MinOM - Record 17: Minimum generator speed (rpm) [used only with Bladed Interface]\n');
fprintf(fid, '          0   GenSpd_MaxOM - Record 18: Optimal mode maximum speed (rpm) [used only with Bladed Interface]\n');
fprintf(fid, '          0   GenSpd_Dem   - Record 19: Demanded generator speed above rated (rpm) [used only with Bladed Interface]\n');
fprintf(fid, '          0   GenTrq_Dem   - Record 22: Demanded generator torque above rated (Nm) [used only with Bladed Interface]\n');
fprintf(fid, '          0   GenPwr_Dem   - Record 13: Demanded power (W) [used only with Bladed Interface]\n');
fprintf(fid, '---------------------- BLADED INTERFACE TORQUE-SPEED LOOK-UP TABLE -------------\n');
fprintf(fid, '          0   DLL_NumTrq   - Record 26: No. of points in torque-speed look-up table {0 = none and use the optimal mode parameters; nonzero = ignore the optimal mode PARAMETERs by setting Record 16 to 0.0} (-) [used only with Bladed Interface]\n');
fprintf(fid, ' GenSpd_TLU   GenTrq_TLU\n');
fprintf(fid, ' (rpm)          (Nm)\n');
fprintf(fid, '---------------------- OUTPUT --------------------------------------------------\n');
fprintf(fid, 'True          SumPrint     - Print summary data to <RootName>.sum (flag) (currently unused)\n');
fprintf(fid, '          1   OutFile      - Switch to determine where output will be placed: {1: in module output file only; 2: in glue code output file only; 3: both} (currently unused)\n');
fprintf(fid, 'True          TabDelim     - Use tab delimiters in text tabular output file? (flag) (currently unused)\n');
fprintf(fid, '"ES10.3E2"    OutFmt       - Format used for text tabular output (except time).  Resulting field should be 10 characters. (quoted string) (currently unused)\n');
fprintf(fid, '          0   TStart       - Time to begin tabular output (s) (currently unused)\n');
fprintf(fid, '              OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n');
fprintf(fid, '"GenPwr"                  - Electrical generator power and torque\n');
fprintf(fid, '"GenTq"                   - Electrical generator power and torque\n');
fprintf(fid, 'END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n');
fprintf(fid, '---------------------------------------------------------------------------------------\n');
fclose(fid);
