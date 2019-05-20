function FASTinput(DT, TMax, varargin)

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


%% FAST input file
fid = fopen([pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'FAST.fst'], 'wt');
fprintf(fid, '------- FAST v8.16.* INPUT FILE ------------------------------------------------\n');
fprintf(fid, 'Created %s.\n', datestr(now));
fprintf(fid, '---------------------- SIMULATION CONTROL --------------------------------------\n');
fprintf(fid, 'False         Echo            - Echo input data to <RootName>.ech (flag)\n');
fprintf(fid, '"FATAL"       AbortLevel      - Error level when simulation should abort (string) {"WARNING", "SEVERE", "FATAL"}\n');
fprintf(fid, '         %s   TMax            - Total run time (s)\n', num2str(TMax));
fprintf(fid, '      %s   DT              - Recommended module time step (s)\n', num2str(DT));
fprintf(fid, '          2   InterpOrder     - Interpolation order for input/output time history (-) {1=linear, 2=quadratic}\n');
fprintf(fid, '          1   NumCrctn        - Number of correction iterations (-) {0=explicit calculation, i.e., no corrections}\n');
fprintf(fid, '    99999.9   DT_UJac         - Time between calls to get Jacobians (s)\n');
fprintf(fid, '      1E+06   UJacSclFact     - Scaling factor used in Jacobians (-)\n');
fprintf(fid, '---------------------- FEATURE SWITCHES AND FLAGS ------------------------------\n');
fprintf(fid, '          1   CompElast       - Compute structural dynamics (switch) {1=ElastoDyn; 2=ElastoDyn + BeamDyn for blades}\n');
fprintf(fid, '          1   CompInflow      - Compute inflow wind velocities (switch) {0=still air; 1=InflowWind; 2=external from OpenFOAM}\n');
fprintf(fid, '          2   CompAero        - Compute aerodynamic loads (switch) {0=None; 1=AeroDyn v14; 2=AeroDyn v15}\n');
fprintf(fid, '          1   CompServo       - Compute control and electrical-drive dynamics (switch) {0=None; 1=ServoDyn}\n');
fprintf(fid, '          0   CompHydro       - Compute hydrodynamic loads (switch) {0=None; 1=HydroDyn}\n');
fprintf(fid, '          0   CompSub         - Compute sub-structural dynamics (switch) {0=None; 1=SubDyn}\n');
fprintf(fid, '          0   CompMooring     - Compute mooring system (switch) {0=None; 1=MAP++; 2=FEAMooring; 3=MoorDyn; 4=OrcaFlex}\n');
fprintf(fid, '          0   CompIce         - Compute ice loads (switch) {0=None; 1=IceFloe; 2=IceDyn}\n');
fprintf(fid, '---------------------- INPUT FILES ---------------------------------------------\n');
fprintf(fid, '"ElastoDyn.dat"   EDFile      - Name of file containing ElastoDyn input parameters (quoted string)\n');
fprintf(fid, '"unused"      BDBldFile(1)    - Name of file containing BeamDyn input parameters for blade 1 (quoted string)\n');
fprintf(fid, '"unused"      BDBldFile(2)    - Name of file containing BeamDyn input parameters for blade 2 (quoted string)\n');
fprintf(fid, '"unused"      BDBldFile(3)    - Name of file containing BeamDyn input parameters for blade 3 (quoted string)\n');
fprintf(fid, '"InflowWind.dat"  InflowFile  - Name of file containing inflow wind input parameters (quoted string)\n');
fprintf(fid, '"AeroDyn.dat"     AeroFile    - Name of file containing aerodynamic input parameters (quoted string)\n');
fprintf(fid, '"ServoDyn.dat"    ServoFile   - Name of file containing control and electrical-drive input parameters (quoted string)\n');
fprintf(fid, '"unused"      HydroFile       - Name of file containing hydrodynamic input parameters (quoted string)\n');
fprintf(fid, '"unused"      SubFile         - Name of file containing sub-structural input parameters (quoted string)\n');
fprintf(fid, '"unused"      MooringFile     - Name of file containing mooring system input parameters (quoted string)\n');
fprintf(fid, '"unused"      IceFile         - Name of file containing ice input parameters (quoted string)\n');
fprintf(fid, '---------------------- OUTPUT --------------------------------------------------\n');
fprintf(fid, 'True          SumPrint        - Print summary data to "<RootName>.sum" (flag)\n');
fprintf(fid, '          1   SttsTime        - Amount of time between screen status messages (s)\n');
fprintf(fid, '      99999   ChkptTime       - Amount of time between creating checkpoint files for potential restart (s)\n');
fprintf(fid, '  "default"   DT_Out          - Time step for tabular output (s) (or "default")\n');
fprintf(fid, '          0   TStart          - Time to begin tabular output (s)\n');
fprintf(fid, '          1   OutFileFmt      - Format for tabular (time-marching) output file (switch) {1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both}\n');
fprintf(fid, 'True          TabDelim        - Use tab delimiters in text tabular output file? (flag) {uses spaces if false}\n');
fprintf(fid, '"ES10.3E2"    OutFmt          - Format used for text tabular output, excluding the time channel.  Resulting field should be 10 characters. (quoted string)\n');
fprintf(fid, '---------------------- LINEARIZATION -------------------------------------------\n');
fprintf(fid, '%s         Linearize       - Linearization analysis (flag)\n', linearize);
fprintf(fid, '          %i   NLinTimes       - Number of times to linearize (-) [>=1] [unused if Linearize=False]\n', NLinTimes);
fprintf(fid, '%s            LinTimes        - List of times at which to linearize (s) [1 to NLinTimes] [unused if Linearize=False]\n', LinTimes);
fprintf(fid, '          1   LinInputs       - Inputs included in linearization (switch) {0=none; 1=standard; 2=all module inputs (debug)} [unused if Linearize=False]\n');
fprintf(fid, '          1   LinOutputs      - Outputs included in linearization (switch) {0=none; 1=from OutList(s); 2=all module outputs (debug)} [unused if Linearize=False]\n');
fprintf(fid, 'False         LinOutJac       - Include full Jacobians in linearization output (for debug) (flag) [unused if Linearize=False; used only if LinInputs=LinOutputs=2]\n');
fprintf(fid, 'False         LinOutMod       - Write module-level linearization output files in addition to output for full system? (flag) [unused if Linearize=False]\n');
fprintf(fid, '---------------------- VISUALIZATION ------------------------------------------\n');
fprintf(fid, '          0   WrVTK           - VTK visualization data output: (switch) {0=none; 1=initialization data only; 2=animation}\n');
fprintf(fid, '          2   VTK_type        - Type of VTK visualization data: (switch) {1=surfaces; 2=basic meshes (lines/points); 3=all meshes (debug)} [unused if WrVTK=0]\n');
fprintf(fid, 'false         VTK_fields      - Write mesh fields to VTK data files? (flag) {true/false} [unused if WrVTK=0]\n');
fprintf(fid, '         15   VTK_fps         - Frame rate for VTK output (frames per second){will use closest integer multiple of DT} [used only if WrVTK=2]\n');
fclose(fid);
