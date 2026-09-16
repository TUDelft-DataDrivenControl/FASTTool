function text = FASTinput(DT, TMax, varargin)
% FASTinput writes data in a character array, which can be written to file.
% The writing of the character array is separated by this function from
% writing the file in FileIOClass.m, to separate agnostic fileIO from the
% generation of the content of the file.
%
%   Syntax
%     text = FASTinput(DT, TMax, varargin)
%
%   Input arguments
%     DT - Time step for the simulation
%     TMax - Duration of the simulation. For linearisation this is the
%     simulation time until the moment for which linearisation is done
%     varargin{1} - Mode of the type of simulation (linearisation or
%     simulation, with simulation as default)
%     varargin{2} - Array with times at which to linearise
%
%   Output
%     text - Character array formatted as a FAST v8.16.* input file that
%     can be written to a plain-text file named FAST.fst.
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

linearize = 'False';
NLinTimes = 1;
LinTimes = '60';
if length(varargin) >= 1
    if contains(varargin{1},'Linearize')
        linearize = 'True';
    end
    T = varargin{2};
    NLinTimes = length(T);
    LinTimes = num2str(T(1));
    for i = 2:NLinTimes
        LinTimes = [LinTimes, ', ', num2str(T(i))];
    end
end

% FAST input file
text = '';
text = append(text, sprintf('------- FAST v8.16.* INPUT FILE ------------------------------------------------\n'));
text = append(text, sprintf('Created %s.\n', datetime));
text = append(text, sprintf('---------------------- SIMULATION CONTROL --------------------------------------\n'));
text = append(text, sprintf('False         Echo            - Echo input data to <RootName>.ech (flag)\n'));
text = append(text, sprintf('"FATAL"       AbortLevel      - Error level when simulation should abort (string) {"WARNING", "SEVERE", "FATAL"}\n'));
text = append(text, sprintf('         %s   TMax            - Total run time (s)\n', num2str(TMax)));
text = append(text, sprintf('      %s   DT              - Recommended module time step (s)\n', num2str(DT)));
text = append(text, sprintf('          2   InterpOrder     - Interpolation order for input/output time history (-) {1=linear, 2=quadratic}\n'));
text = append(text, sprintf('          1   NumCrctn        - Number of correction iterations (-) {0=explicit calculation, i.e., no corrections}\n'));
text = append(text, sprintf('    99999.9   DT_UJac         - Time between calls to get Jacobians (s)\n'));
text = append(text, sprintf('      1E+06   UJacSclFact     - Scaling factor used in Jacobians (-)\n'));
text = append(text, sprintf('---------------------- FEATURE SWITCHES AND FLAGS ------------------------------\n'));
text = append(text, sprintf('          1   CompElast       - Compute structural dynamics (switch) {1=ElastoDyn; 2=ElastoDyn + BeamDyn for blades}\n'));
text = append(text, sprintf('          1   CompInflow      - Compute inflow wind velocities (switch) {0=still air; 1=InflowWind; 2=external from OpenFOAM}\n'));
text = append(text, sprintf('          2   CompAero        - Compute aerodynamic loads (switch) {0=None; 1=AeroDyn v14; 2=AeroDyn v15}\n'));
text = append(text, sprintf('          1   CompServo       - Compute control and electrical-drive dynamics (switch) {0=None; 1=ServoDyn}\n'));
text = append(text, sprintf('          0   CompHydro       - Compute hydrodynamic loads (switch) {0=None; 1=HydroDyn}\n'));
text = append(text, sprintf('          0   CompSub         - Compute sub-structural dynamics (switch) {0=None; 1=SubDyn}\n'));
text = append(text, sprintf('          0   CompMooring     - Compute mooring system (switch) {0=None; 1=MAP++; 2=FEAMooring; 3=MoorDyn; 4=OrcaFlex}\n'));
text = append(text, sprintf('          0   CompIce         - Compute ice loads (switch) {0=None; 1=IceFloe; 2=IceDyn}\n'));
text = append(text, sprintf('---------------------- INPUT FILES ---------------------------------------------\n'));
text = append(text, sprintf('"ElastoDyn.dat"   EDFile      - Name of file containing ElastoDyn input parameters (quoted string)\n'));
text = append(text, sprintf('"unused"      BDBldFile(1)    - Name of file containing BeamDyn input parameters for blade 1 (quoted string)\n'));
text = append(text, sprintf('"unused"      BDBldFile(2)    - Name of file containing BeamDyn input parameters for blade 2 (quoted string)\n'));
text = append(text, sprintf('"unused"      BDBldFile(3)    - Name of file containing BeamDyn input parameters for blade 3 (quoted string)\n'));
text = append(text, sprintf('"InflowWind.dat"  InflowFile  - Name of file containing inflow wind input parameters (quoted string)\n'));
text = append(text, sprintf('"AeroDyn.dat"     AeroFile    - Name of file containing aerodynamic input parameters (quoted string)\n'));
text = append(text, sprintf('"ServoDyn.dat"    ServoFile   - Name of file containing control and electrical-drive input parameters (quoted string)\n'));
text = append(text, sprintf('"unused"      HydroFile       - Name of file containing hydrodynamic input parameters (quoted string)\n'));
text = append(text, sprintf('"unused"      SubFile         - Name of file containing sub-structural input parameters (quoted string)\n'));
text = append(text, sprintf('"unused"      MooringFile     - Name of file containing mooring system input parameters (quoted string)\n'));
text = append(text, sprintf('"unused"      IceFile         - Name of file containing ice input parameters (quoted string)\n'));
text = append(text, sprintf('---------------------- OUTPUT --------------------------------------------------\n'));
text = append(text, sprintf('True          SumPrint        - Print summary data to "<RootName>.sum" (flag)\n'));
text = append(text, sprintf('          1   SttsTime        - Amount of time between screen status messages (s)\n'));
text = append(text, sprintf('      99999   ChkptTime       - Amount of time between creating checkpoint files for potential restart (s)\n'));
text = append(text, sprintf('  "default"   DT_Out          - Time step for tabular output (s) (or "default")\n'));
text = append(text, sprintf('          0   TStart          - Time to begin tabular output (s)\n'));
text = append(text, sprintf('          1   OutFileFmt      - Format for tabular (time-marching) output file (switch) {1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both}\n'));
text = append(text, sprintf('True          TabDelim        - Use tab delimiters in text tabular output file? (flag) {uses spaces if false}\n'));
text = append(text, sprintf('"ES10.3E2"    OutFmt          - Format used for text tabular output, excluding the time channel.  Resulting field should be 10 characters. (quoted string)\n'));
text = append(text, sprintf('---------------------- LINEARIZATION -------------------------------------------\n'));
text = append(text, sprintf('%s         Linearize       - Linearization analysis (flag)\n', linearize));
text = append(text, sprintf('          %i   NLinTimes       - Number of times to linearize (-) [>=1] [unused if Linearize=False]\n', NLinTimes));
text = append(text, sprintf('%s            LinTimes        - List of times at which to linearize (s) [1 to NLinTimes] [unused if Linearize=False]\n', LinTimes));
text = append(text, sprintf('          1   LinInputs       - Inputs included in linearization (switch) {0=none; 1=standard; 2=all module inputs (debug)} [unused if Linearize=False]\n'));
text = append(text, sprintf('          1   LinOutputs      - Outputs included in linearization (switch) {0=none; 1=from OutList(s)); 2=all module outputs (debug)} [unused if Linearize=False]\n'));
text = append(text, sprintf('False         LinOutJac       - Include full Jacobians in linearization output (for debug) (flag) [unused if Linearize=False; used only if LinInputs=LinOutputs=2]\n'));
text = append(text, sprintf('False         LinOutMod       - Write module-level linearization output files in addition to output for full system? (flag) [unused if Linearize=False]\n'));
text = append(text, sprintf('---------------------- VISUALIZATION ------------------------------------------\n'));
text = append(text, sprintf('          0   WrVTK           - VTK visualization data output: (switch) {0=none; 1=initialization data only; 2=animation}\n'));
text = append(text, sprintf('          2   VTK_type        - Type of VTK visualization data: (switch) {1=surfaces; 2=basic meshes (lines/points)); 3=all meshes (debug)} [unused if WrVTK=0]\n'));
text = append(text, sprintf('false         VTK_fields      - Write mesh fields to VTK data files? (flag) {true/false} [unused if WrVTK=0]\n'));
text = append(text, sprintf('         15   VTK_fps         - Frame rate for VTK output (frames per second){will use closest integer multiple of DT} [used only if WrVTK=2]\n'));
